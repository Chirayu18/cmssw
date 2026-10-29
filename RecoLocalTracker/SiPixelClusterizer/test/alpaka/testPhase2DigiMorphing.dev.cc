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
  std::vector<int> referenceClusters(Module const& module, int* nFakes = nullptr) {
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

    if (nFakes)
      *nFakes = std::count(fake.begin(), fake.end(), 1);

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

  // The labelling `fine` is a refinement of `coarse`: pixels in the same cluster in `fine` are in the same in `coarse`
  bool refines(std::vector<int> const& fine, std::vector<int> const& coarse) {
    std::map<int, int> map;
    for (size_t i = 0; i < fine.size(); ++i) {
      auto [it, inserted] = map.emplace(fine[i], coarse[i]);
      if (it->second != coarse[i])
        return false;
    }
    return true;
  }

  size_t countClusters(std::vector<int> const& v) { return std::set<int>(v.begin(), v.end()).size(); }

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

  // Result of FindClus for every digi, in the same order as the pixels of the modules
  struct DeviceResult {
    std::vector<std::vector<int>> clus;
    std::vector<std::vector<uint32_t>> moduleId;
  };

  // Run CountModules and FindClus on the modules. The morphing modules are those with `morphing` set, plus the raw
  // ids in `extraMorphing`, which have no digis.
  DeviceResult runFindClus(Device const& device,
                           std::vector<Module> const& modules,
                           uint32_t maxFakes,
                           std::vector<std::pair<uint32_t, uint32_t>> extraMorphing = {},
                           bool printAttributes = false) {
    std::vector<std::pair<uint32_t, uint32_t>> morphing = extraMorphing;
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

    uint32_t nDigis = 0;
    for (auto const& m : modules)
      nDigis += m.pixels.size();

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

    if (printAttributes) {
      auto attributes =
          alpaka::getFunctionAttributes<Acc1D>(device,
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
                << " threads per block, "
                << ALPAKA_ACCELERATOR_NAMESPACE::pixelClustering::FindClus<TrackerTraits>::maxElementsPerBlockMorph
                << " needed\n";
    }
    const uint32_t threads = 256;
    alpaka::exec<Acc1D>(queue,
                        cms::alpakatools::make_workdiv<Acc1D>(cms::alpakatools::divide_up_by(nDigis, threads), threads),
                        ALPAKA_ACCELERATOR_NAMESPACE::pixelClustering::CountModules<TrackerTraits>{},
                        digis_d.view(),
                        clusters_d.view(),
                        nDigis);
    alpaka::exec<Acc1D>(
        queue,
        cms::alpakatools::make_workdiv<Acc1D>(
            modules.size(),
            ALPAKA_ACCELERATOR_NAMESPACE::pixelClustering::FindClus<TrackerTraits>::maxElementsPerBlockMorph),
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

    DeviceResult out;
    i = 0;
    for (auto const& m : modules) {
      out.clus.emplace_back();
      out.moduleId.emplace_back();
      for (size_t k = 0; k < m.pixels.size(); ++k, ++i) {
        out.clus.back().push_back(result.view()[i].clus());
        out.moduleId.back().push_back(result.view()[i].moduleId());
      }
    }
    return out;
  }

  // A module with the given pixels, in random order
  Module exactModule(std::mt19937& rng, uint32_t rawId, uint16_t rows, uint16_t cols, std::vector<Pixel> pixels) {
    std::set<std::pair<int, int>> seen;
    Module module{rawId, rows, cols, true, {}};
    for (auto const& p : pixels) {
      REQUIRE(p.x < rows);
      REQUIRE(p.y < cols);
      if (seen.emplace(p.x, p.y).second)
        module.pixels.push_back(p);
    }
    std::shuffle(module.pixels.begin(), module.pixels.end(), rng);
    return module;
  }

  // A track along the columns at row x, from column y0 to y1 included, with the columns in `gaps` missing
  void addAlongColumns(std::vector<Pixel>& pixels, int x, int y0, int y1, std::set<int> const& gaps) {
    for (int y = y0; y <= y1; ++y)
      if (not gaps.contains(y))
        pixels.push_back({static_cast<uint16_t>(x), static_cast<uint16_t>(y)});
  }

  // A track along the rows at column y, from row x0 to x1 included, with the rows in `gaps` missing
  void addAlongRows(std::vector<Pixel>& pixels, int y, int x0, int x1, std::set<int> const& gaps) {
    for (int x = x0; x <= x1; ++x)
      if (not gaps.contains(x))
        pixels.push_back({static_cast<uint16_t>(x), static_cast<uint16_t>(y)});
  }

  // The 8 pixels around (x, y), leaving (x, y) empty
  void addRing(std::vector<Pixel>& pixels, int x, int y) {
    for (int dx : {-1, 0, 1})
      for (int dy : {-1, 0, 1})
        if (dx != 0 or dy != 0)
          pixels.push_back({static_cast<uint16_t>(x + dx), static_cast<uint16_t>(y + dy)});
  }

  // Compare the device clusters of every module with the reference, and check that the morphing did something
  void checkAgainstReference(std::vector<Module> const& modules, DeviceResult const& result) {
    for (size_t m = 0; m < modules.size(); ++m) {
      auto const& module = modules[m];
      int nFakes = 0;
      auto reference = referenceClusters(module, &nFakes);
      Module plain = module;
      plain.morphing = false;
      auto unmorphed = referenceClusters(plain);
      std::cout << "  module " << module.rawId << " (" << module.rows << " x " << module.cols << ", morphing "
                << (module.morphing ? "on" : "off") << "): " << module.pixels.size() << " pixels, " << nFakes
                << " recovered, " << countClusters(unmorphed) << " clusters before morphing, "
                << countClusters(reference) << " expected, " << countClusters(result.clus[m]) << " found\n";
      CHECK(samePartition(result.clus[m], reference));
    }
  }

}  // namespace

TEST_CASE("Phase-2 digi morphing on " EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE), "[PixelClustering]") {
  auto const& devices = cms::alpakatools::devices<Platform>();
  if (devices.empty()) {
    FAIL("No devices available for the " EDM_STRINGIZE(ALPAKA_ACCELERATOR_NAMESPACE) " backend, "
         "the test will be skipped.\n");
  }
  const uint32_t maxFakes = ::pixelClustering::phase2Morphing::defaultMaxFakesInModule;

  SECTION("random broken tracks") {
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

    for (auto const& device : devices) {
      std::cout << "Phase-2 digi morphing, random broken tracks, on " << alpaka::getName(device) << '\n';
      checkAgainstReference(modules, runFindClus(device, modules, maxFakes, {}, true));
    }
  }

  SECTION("strip, word and module boundaries") {
    std::mt19937 rng(7);
    std::vector<Module> modules;

    // TBPX L1 size, 216 = 4 x 44 + 40 columns: one missing pixel on every strip boundary, next to the first and
    // last columns, and a missing pixel surrounded by 8 pixels on the strip boundaries
    {
      std::vector<Pixel> pixels;
      int x = 4;
      for (int gap : {1, 43, 44, 45, 87, 88, 131, 132, 175, 176, 214}) {
        addAlongColumns(pixels, x, std::max(gap - 3, 0), std::min(gap + 3, 215), {gap});
        x += 6;
      }
      // two missing pixels across a strip boundary
      addAlongColumns(pixels, x, 40, 48, {43, 44});
      x += 6;
      for (int y : {43, 44, 87, 88, 175, 176, 214}) {
        addRing(pixels, x + 1, y);
        x += 6;
      }
      // missing pixels on the word boundaries, along the rows
      for (int gap : {1, 31, 32, 63, 64, 670}) {
        addAlongRows(pixels, 20 + gap % 7, std::max(gap - 3, 0), std::min(gap + 3, 671), {gap});
        addAlongRows(pixels, 66 + gap % 5, std::max(gap - 3, 0), std::min(gap + 3, 671), {gap});
      }
      // a track along the last column and one along the first, with a missing pixel
      addAlongRows(pixels, 215, 200, 210, {205});
      addAlongRows(pixels, 0, 200, 210, {205});
      modules.push_back(exactModule(rng, 4002, 672, 216, pixels));
    }
    // quad module with the big pixels: 1354 x 434 = 9 x 44 + 38 columns, and 1354 = 42 x 32 + 10 rows
    {
      std::vector<Pixel> pixels;
      for (int gap : {1343, 1344, 1351, 1352})
        addAlongRows(pixels, 100 + gap % 13, gap - 3, std::min(gap + 3, 1353), {gap});
      addRing(pixels, 1352, 300);
      addRing(pixels, 1343, 301);
      int x = 600;
      for (int gap : {395, 396, 397, 431, 432}) {
        addAlongColumns(pixels, x, gap - 3, std::min(gap + 3, 433), {gap});
        x += 6;
      }
      addRing(pixels, 1352, 432);
      addRing(pixels, 1, 1);
      modules.push_back(exactModule(rng, 4001, 1354, 434, pixels));
    }
    // one pixel, and two pixels in opposite corners
    modules.push_back(exactModule(rng, 4004, 672, 434, {{300, 200}}));
    modules.push_back(exactModule(rng, 4003, 1354, 434, {{0, 0}, {1353, 433}}));

    for (auto const& device : devices) {
      std::cout << "Phase-2 digi morphing, boundaries, on " << alpaka::getName(device) << '\n';
      // 4000 is a morphing module without digis
      auto result =
          runFindClus(device, modules, maxFakes, {{4000, ::pixelClustering::phase2Morphing::packSize(672, 216)}});
      checkAgainstReference(modules, result);
      // the test is only useful if the morphing changes the clusters
      for (size_t m = 0; m < 2; ++m) {
        Module plain = modules[m];
        plain.morphing = false;
        CHECK(countClusters(referenceClusters(plain)) > countClusters(referenceClusters(modules[m])));
      }
    }
  }

  SECTION("more recovered pixels than maxFakesInModule") {
    // With fewer slots than recovered pixels, the slots go to the pixels found first, which is not deterministic.
    // The clusters must then lie between those without morphing and those with the full morphing.
    std::mt19937 rng(11);
    std::vector<Module> modules;
    modules.push_back(makeModule(rng, 5001, 672, 432, true));
    modules.push_back(makeModule(rng, 5002, 1344, 432, true));
    const uint32_t fewFakes = 5;

    for (auto const& device : devices) {
      std::cout << "Phase-2 digi morphing, " << fewFakes << " recovered pixels per module, on "
                << alpaka::getName(device) << '\n';
      auto first = runFindClus(device, modules, fewFakes);
      auto second = runFindClus(device, modules, fewFakes);
      for (size_t m = 0; m < modules.size(); ++m) {
        int nFakes = 0;
        auto reference = referenceClusters(modules[m], &nFakes);
        Module plain = modules[m];
        plain.morphing = false;
        auto unmorphed = referenceClusters(plain);
        REQUIRE(nFakes > static_cast<int>(fewFakes));
        std::cout << "  module " << modules[m].rawId << ": " << nFakes << " recoverable, " << countClusters(unmorphed)
                  << " clusters before morphing, " << countClusters(reference) << " with the full morphing, "
                  << countClusters(first.clus[m]) << " and " << countClusters(second.clus[m]) << " found in two runs, "
                  << (samePartition(first.clus[m], second.clus[m]) ? "same" : "different") << " clusters\n";
        CHECK(refines(unmorphed, first.clus[m]));
        CHECK(refines(first.clus[m], reference));
        CHECK(refines(unmorphed, second.clus[m]));
        CHECK(refines(second.clus[m], reference));
      }
    }
  }

  SECTION("more valid pixels than the histogram keeps for them") {
    // In the morphing modules maxFakesInModule of the maxPixInModule histogram entries are kept for the recovered
    // pixels: the valid pixels past the first maxPixInModule - maxFakesInModule are dropped, without morphing they
    // are all kept.
    std::mt19937 rng(3);
    const uint32_t keep = TrackerTraits::maxPixInModule - maxFakes;
    const uint32_t nPixels = keep + 500;
    // isolated pixels 4 rows and columns apart, which the morphing does not connect, and for one in 20 a second
    // pixel 2 columns away, connected by one recovered pixel: random pixels at this occupancy would give more
    // recovered pixels than maxFakesInModule
    std::uniform_int_distribution<int> row(0, 1343 / 4), col(0, 431 / 4 - 1);
    std::vector<Pixel> pixels;
    std::set<std::pair<int, int>> seen;
    while (pixels.size() < nPixels) {
      int x = 4 * row(rng), y = 4 * col(rng);
      if (not seen.emplace(x, y).second)
        continue;
      pixels.push_back({static_cast<uint16_t>(x), static_cast<uint16_t>(y)});
      if (seen.size() % 20 == 0 and pixels.size() < nPixels)
        pixels.push_back({static_cast<uint16_t>(x), static_cast<uint16_t>(y + 2)});
    }
    std::vector<Module> modules;
    modules.push_back(Module{6001, 1344, 432, true, pixels});
    modules.push_back(Module{6002, 1344, 432, false, pixels});

    Module kept = modules[0];
    kept.pixels.resize(keep);
    int nFakes = 0;
    auto reference = referenceClusters(kept, &nFakes);
    REQUIRE(nFakes <= static_cast<int>(maxFakes));

    for (auto const& device : devices) {
      std::cout << "Phase-2 digi morphing, " << nPixels << " valid pixels, on " << alpaka::getName(device) << '\n';
      auto result = runFindClus(device, modules, maxFakes);
      std::vector<int> keptClus(result.clus[0].begin(), result.clus[0].begin() + keep);
      CHECK(samePartition(keptClus, reference));
      uint32_t dropped = 0;
      for (uint32_t k = keep; k < nPixels; ++k)
        dropped += (result.clus[0][k] == ::pixelClustering::invalidClusterId and
                    result.moduleId[0][k] == ::pixelClustering::invalidModuleId);
      CHECK(dropped == nPixels - keep);
      std::cout << "  morphing module: " << keep << " pixels kept, " << dropped << " dropped, " << nFakes
                << " recovered\n";
      Module plain = modules[1];
      CHECK(samePartition(result.clus[1], referenceClusters(plain)));
    }
  }
}
