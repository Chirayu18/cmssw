// Check the Phase-2 digi morphing in FindClus against a whole-module reference implementation on the host.
//
// The device processes each module one strip of columns at a time. The reference builds the full module bitmap,
// applies the same morphological closing (3x3 dilation, then erosion with a plus-shaped kernel), and clusters the
// valid and recovered pixels with 8-neighbour connectivity. The two must assign the valid pixels to the same
// clusters, in every module.

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <vector>

#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include <alpaka/alpaka.hpp>

#include "DataFormats/SiPixelClusterSoA/interface/alpaka/SiPixelClustersSoACollection.h"
#include "DataFormats/SiPixelDigiSoA/interface/SiPixelDigisHost.h"
#include "DataFormats/SiPixelDigiSoA/interface/alpaka/SiPixelDigisSoACollection.h"
#include "FWCore/Utilities/interface/stringize.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "RecoLocalTracker/SiPixelClusterizer/plugins/alpaka/PixelClustering.h"

using namespace ALPAKA_ACCELERATOR_NAMESPACE;
using TrackerTraits = pixelTopology::Phase2;

namespace {

  struct Pixel {
    uint16_t x;  // row
    uint16_t y;  // column
  };

  struct Module {
    uint32_t rawId;
    uint16_t rows;
    uint16_t cols;
    bool morphing;
    std::vector<Pixel> pixels;
  };

  // Reference: cluster index of each pixel of the module, after the (optional) whole-module digi morphing
  std::vector<int> referenceClusters(Module const& module) {
    const int rows = module.rows;
    const int cols = module.cols;
    auto at = [rows](int x, int y) { return y * rows + x; };
    std::vector<uint8_t> image(rows * cols, 0);
    for (auto const& p : module.pixels)
      image[at(p.x, p.y)] = 1;

    std::vector<uint8_t> fake(rows * cols, 0);
    if (module.morphing) {
      // dilate: empty pixels with at least one of the 8 neighbours set
      std::vector<uint8_t> dilated(rows * cols, 0);
      for (int y = 0; y < cols; ++y)
        for (int x = 0; x < rows; ++x) {
          if (image[at(x, y)])
            continue;
          for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx) {
              int nx = x + dx, ny = y + dy;
              if (nx >= 0 and nx < rows and ny >= 0 and ny < cols and image[at(nx, ny)])
                dilated[at(x, y)] = 1;
            }
        }
      // erode: keep the dilated pixels whose 4 neighbours are not empty;
      // neighbours beyond the first or last row count as not empty, beyond the first or last column as empty
      auto filled = [&](int x, int y) {
        if (x < 0 or x >= rows)
          return true;
        if (y < 0 or y >= cols)
          return false;
        return image[at(x, y)] or dilated[at(x, y)];
      };
      for (int y = 0; y < cols; ++y)
        for (int x = 0; x < rows; ++x)
          if (dilated[at(x, y)] and filled(x - 1, y) and filled(x + 1, y) and filled(x, y - 1) and filled(x, y + 1))
            fake[at(x, y)] = 1;
    }

    // union-find over the valid and recovered pixels, with 8-neighbour connectivity
    std::vector<int> parent(rows * cols);
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](int i) {
      while (parent[i] != i)
        i = parent[i] = parent[parent[i]];
      return i;
    };
    auto used = [&](int x, int y) {
      return x >= 0 and x < rows and y >= 0 and y < cols and (image[at(x, y)] or fake[at(x, y)]);
    };
    for (int y = 0; y < cols; ++y)
      for (int x = 0; x < rows; ++x) {
        if (not used(x, y))
          continue;
        for (int dy = -1; dy <= 1; ++dy)
          for (int dx = -1; dx <= 1; ++dx)
            if (used(x + dx, y + dy))
              parent[find(at(x, y))] = find(at(x + dx, y + dy));
      }

    std::vector<int> clusters;
    clusters.reserve(module.pixels.size());
    for (auto const& p : module.pixels)
      clusters.push_back(find(at(p.x, p.y)));
    return clusters;
  }

  // Two labellings describe the same partition if the mapping between them is one-to-one
  bool samePartition(std::vector<int> const& a, std::vector<int> const& b) {
    std::map<int, int> ab, ba;
    for (size_t i = 0; i < a.size(); ++i) {
      auto [ia, newa] = ab.emplace(a[i], b[i]);
      auto [ib, newb] = ba.emplace(b[i], a[i]);
      if (ia->second != b[i] or ib->second != a[i])
        return false;
    }
    return true;
  }

  // Broken tracks: straight segments with missing pixels, which the morphing should bridge, plus edge cases
  Module makeModule(std::mt19937& rng, uint32_t rawId, uint16_t rows, uint16_t cols, bool morphing) {
    Module module{rawId, rows, cols, morphing, {}};
    std::vector<uint8_t> occupied(rows * cols, 0);
    auto add = [&](int x, int y) {
      if (x < 0 or x >= rows or y < 0 or y >= cols or occupied[y * rows + x])
        return;
      occupied[y * rows + x] = 1;
      module.pixels.push_back({static_cast<uint16_t>(x), static_cast<uint16_t>(y)});
    };

    std::uniform_int_distribution<int> row(0, rows - 1), col(0, cols - 1), length(2, 40), width(1, 3);
    std::bernoulli_distribution hole(0.15);
    // about 30 broken tracks in a dual module, keeping well below maxPixInModule
    const int tracks = std::max(2, rows * cols / 10000);
    for (int track = 0; track < tracks; ++track) {
      // segments along the columns (long clusters, as at high |eta| in the barrel) and along the rows
      bool alongColumns = track % 3 != 0;
      int x0 = row(rng), y0 = col(rng), n = length(rng), w = width(rng);
      for (int k = 0; k < n; ++k)
        for (int j = 0; j < w; ++j)
          if (not hole(rng))
            alongColumns ? add(x0 + j, y0 + k) : add(x0 + k, y0 + j);
    }
    // single isolated pixels
    for (int i = 0; i < 100; ++i)
      add(row(rng), col(rng));
    // pixels on the module edges, and across the strip boundaries
    for (int y : {0, 1, 2, cols - 3, cols - 2, cols - 1}) {
      add(0, y);
      add(2, y);
      add(rows - 1, y);
      add(rows - 3, y);
    }
    for (int y = 28; y < cols; y += 32)
      for (int k = 0; k < 9; k += 2) {
        add(rows / 2, y + k);
        add(rows / 2 + 1, y + k);
        add(rows / 3 + (k % 4), y + k);
      }
    // a hole surrounded on all sides, at a word boundary in x
    for (int dx : {-1, 0, 1})
      for (int dy : {-1, 0, 1})
        if (dx != 0 or dy != 0)
          add(64 + dx, 100 + dy);

    std::shuffle(module.pixels.begin(), module.pixels.end(), rng);
    return module;
  }

}  // namespace

TEST_CASE("Phase-2 digi morphing on " EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE), "[PixelClustering]") {
  auto const& devices = cms::alpakatools::devices<Platform>();
  if (devices.empty()) {
    FAIL("No devices available for the " EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE) " backend, "
         "the test will be skipped.\n");
  }

  // dual and quad modules, a module with a number of rows that is not a multiple of 32, and a module
  // without morphing; the raw ids are not in the same order as the modules
  std::mt19937 rng(42);
  std::vector<Module> modules;
  modules.push_back(makeModule(rng, 3005, 672, 432, true));
  modules.push_back(makeModule(rng, 3001, 1344, 432, true));
  modules.push_back(makeModule(rng, 3003, 100, 200, true));
  modules.push_back(makeModule(rng, 3002, 672, 432, false));
  modules.push_back(makeModule(rng, 3004, 1344, 432, true));
  modules.push_back(makeModule(rng, 3000, 672, 70, true));

  // sorted raw ids and sizes of the morphing modules
  std::vector<std::pair<uint32_t, uint32_t>> morphing;
  for (auto const& m : modules)
    if (m.morphing)
      morphing.emplace_back(m.rawId, ::pixelClustering::phase2Morphing::packSize(m.rows, m.cols));
  std::sort(morphing.begin(), morphing.end());
  std::vector<uint32_t> morphingIds, morphingSizes;
  for (auto const& [id, size] : morphing) {
    morphingIds.push_back(id);
    morphingSizes.push_back(size);
  }
  const uint32_t nMorphing = morphingIds.size();
  const uint32_t maxFakes = ::pixelClustering::phase2Morphing::defaultMaxFakesInModule;

  uint32_t nDigis = 0;
  for (auto const& m : modules)
    nDigis += m.pixels.size();

  for (auto const& device : devices) {
    std::cout << "Phase-2 digi morphing on " << alpaka::getName(device) << '\n';
    Queue queue(device);

    // digis, grouped by module
    SiPixelDigisHost digis_h(queue, nDigis);
    uint32_t i = 0;
    for (uint16_t moduleId = 0; moduleId < modules.size(); ++moduleId) {
      for (auto const& p : modules[moduleId].pixels) {
        auto digi = digis_h.view()[i++];
        digi.clus() = 0;
        digi.pdigi() = 0;
        digi.rawIdArr() = modules[moduleId].rawId;
        digi.adc() = 100;
        digi.xx() = p.x;
        digi.yy() = p.y;
        digi.moduleId() = moduleId;
        digi.rawADC() = 0;
      }
    }
    SiPixelDigisSoACollection digis_d(queue, nDigis);
    alpaka::memcpy(queue, digis_d.buffer(), digis_h.buffer());

    SiPixelClustersSoACollection clusters_d(queue, TrackerTraits::numberOfModules);
    clusters_d.zeroInitialise(queue);

    auto morphingIds_d = cms::alpakatools::make_device_buffer<uint32_t[]>(queue, nMorphing);
    auto morphingSizes_d = cms::alpakatools::make_device_buffer<uint32_t[]>(queue, nMorphing);
    alpaka::memcpy(queue, morphingIds_d, cms::alpakatools::make_host_view(morphingIds.data(), nMorphing));
    alpaka::memcpy(queue, morphingSizes_d, cms::alpakatools::make_host_view(morphingSizes.data(), nMorphing));
    SiPixelDigisSoACollection fakes_d(queue, nMorphing * maxFakes);

    {
      auto attributes = alpaka::getFunctionAttributes<Acc1D>(device,
                                                             ALPAKA_ACCELERATOR_NAMESPACE::pixelClustering::FindClus<TrackerTraits>{},
                                                             digis_d.view(),
                                                             fakes_d.view(),
                                                             true,
                                                             morphingIds_d.data(),
                                                             static_cast<uint32_t const*>(morphingSizes_d.data()),
                                                             nMorphing,
                                                             maxFakes,
                                                             clusters_d.view(),
                                                             nDigis);
      std::cout << "  FindClus: " << attributes.numRegs << " registers, at most " << attributes.maxThreadsPerBlock
                << " threads per block, " << ALPAKA_ACCELERATOR_NAMESPACE::pixelClustering::FindClus<TrackerTraits>::maxElementsPerBlockMorph
                << " needed\n";
    }
    const uint32_t threads = 256;
    alpaka::exec<Acc1D>(queue,
                        cms::alpakatools::make_workdiv<Acc1D>(cms::alpakatools::divide_up_by(nDigis, threads), threads),
                        ALPAKA_ACCELERATOR_NAMESPACE::pixelClustering::CountModules<TrackerTraits>{},
                        digis_d.view(),
                        clusters_d.view(),
                        nDigis);
    alpaka::exec<Acc1D>(queue,
                        cms::alpakatools::make_workdiv<Acc1D>(
                            modules.size(), ALPAKA_ACCELERATOR_NAMESPACE::pixelClustering::FindClus<TrackerTraits>::maxElementsPerBlockMorph),
                        ALPAKA_ACCELERATOR_NAMESPACE::pixelClustering::FindClus<TrackerTraits>{},
                        digis_d.view(),
                        fakes_d.view(),
                        true,
                        morphingIds_d.data(),
                        morphingSizes_d.data(),
                        nMorphing,
                        maxFakes,
                        clusters_d.view(),
                        nDigis);

    SiPixelDigisHost result(queue, nDigis);
    alpaka::memcpy(queue, result.buffer(), digis_d.buffer());
    alpaka::wait(queue);

    i = 0;
    for (auto const& m : modules) {
      std::vector<int> device(m.pixels.size());
      for (auto& c : device)
        c = result.view()[i++].clus();
      auto reference = referenceClusters(m);
      Module plain = m;
      plain.morphing = false;
      auto unmorphed = referenceClusters(plain);
      auto count = [](std::vector<int> const& v) { return std::set<int>(v.begin(), v.end()).size(); };
      std::cout << "  module " << m.rawId << " (" << m.rows << " x " << m.cols << ", morphing "
                << (m.morphing ? "on" : "off") << "): " << m.pixels.size() << " pixels, " << count(unmorphed)
                << " clusters before morphing, " << count(reference) << " expected, " << count(device)
                << " found\n";
      CHECK(samePartition(device, reference));
    }
  }
}
