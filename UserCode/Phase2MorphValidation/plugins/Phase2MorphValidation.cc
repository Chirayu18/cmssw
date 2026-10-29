// Validation of the Phase-2 digi morphing in the Alpaka pixel clusterizer (not for the PR).
//
// Phase2MorphValidation compares, module by module:
//  - the pixels recovered on the device with those added by the CPU SiPixelDigiMorphing;
//  - the device clusters (morphing on) with an 8-neighbour clustering of the CPU morphed digis,
//    and the device clusters (morphing off) with an 8-neighbour clustering of the input digis.
// Phase2ClusDump writes the device cluster index and charge of every digi, to compare two releases.

#include <algorithm>
#include <cstdint>
#include <map>
#include <numeric>
#include <set>
#include <unordered_map>
#include <vector>

#include "TTree.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDigiSoA/interface/SiPixelDigisHost.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

namespace {

  using Key = uint32_t;
  inline Key key(int x, int y) { return (static_cast<uint32_t>(x) << 16) | static_cast<uint32_t>(y); }

  // 8-neighbour connected components of the given pixels; returns the component of each pixel
  std::vector<int> connected(std::vector<std::pair<int, int>> const& pixels) {
    std::unordered_map<Key, int> index;
    for (size_t i = 0; i < pixels.size(); ++i)
      index.emplace(key(pixels[i].first, pixels[i].second), i);
    std::vector<int> parent(pixels.size());
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](int i) {
      while (parent[i] != i)
        i = parent[i] = parent[parent[i]];
      return i;
    };
    for (size_t i = 0; i < pixels.size(); ++i) {
      auto [x, y] = pixels[i];
      for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy) {
          if (x + dx < 0 or y + dy < 0)
            continue;
          auto it = index.find(key(x + dx, y + dy));
          if (it != index.end())
            parent[find(i)] = find(it->second);
        }
    }
    std::vector<int> label(pixels.size());
    for (size_t i = 0; i < pixels.size(); ++i)
      label[i] = find(i);
    return label;
  }

  // number of pixels that would have to change cluster to turn partition a into b
  int mismatches(std::vector<int> const& a, std::vector<int> const& b) {
    std::map<int, std::map<int, int>> overlap;
    for (size_t i = 0; i < a.size(); ++i)
      ++overlap[a[i]][b[i]];
    int matched = 0;
    for (auto const& [la, m] : overlap) {
      int best = 0;
      for (auto const& [lb, n] : m)
        best = std::max(best, n);
      matched += best;
    }
    // count both directions, so that both splits and merges are seen
    std::map<int, std::map<int, int>> reverse;
    for (size_t i = 0; i < a.size(); ++i)
      ++reverse[b[i]][a[i]];
    int matchedReverse = 0;
    for (auto const& [lb, m] : reverse) {
      int best = 0;
      for (auto const& [la, n] : m)
        best = std::max(best, n);
      matchedReverse += best;
    }
    return static_cast<int>(a.size()) - std::min(matched, matchedReverse);
  }

  // Exact whole-module digi morphing: 3x3 dilation, then erosion with a plus-shaped kernel.
  // rowEdgeFilled: a neighbour beyond the first or last row counts as filled (as in the Alpaka clusterizer),
  // otherwise as empty (as in the CPU SiPixelDigiMorphing). Beyond the first or last column it is always empty.
  std::set<Key> exactMorphing(std::vector<std::pair<int, int>> const& real, int rows, int cols, bool rowEdgeFilled) {
    std::vector<uint8_t> image(rows * cols, 0), dilated(rows * cols, 0);
    auto at = [rows](int x, int y) { return y * rows + x; };
    for (auto const& [x, y] : real)
      image[at(x, y)] = 1;
    for (auto const& [x0, y0] : real)
      for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy) {
          int x = x0 + dx, y = y0 + dy;
          if (x >= 0 and x < rows and y >= 0 and y < cols and not image[at(x, y)])
            dilated[at(x, y)] = 1;
        }
    auto filled = [&](int x, int y) {
      if (x < 0 or x >= rows)
        return rowEdgeFilled;
      if (y < 0 or y >= cols)
        return false;
      return image[at(x, y)] or dilated[at(x, y)];
    };
    std::set<Key> fakes;
    for (int y = 0; y < cols; ++y)
      for (int x = 0; x < rows; ++x)
        if (dilated[at(x, y)] and filled(x - 1, y) and filled(x + 1, y) and filled(x, y - 1) and filled(x, y + 1))
          fakes.insert(key(x, y));
    return fakes;
  }

  int distinct(std::vector<int> const& v) { return std::set<int>(v.begin(), v.end()).size(); }

}  // namespace

class Phase2MorphValidation : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit Phase2MorphValidation(edm::ParameterSet const& config)
      : digisToken_(consumes(config.getParameter<edm::InputTag>("digis"))),
        cpuMorphedToken_(consumes(config.getParameter<edm::InputTag>("cpuMorphed"))),
        gpuOnToken_(consumes(config.getParameter<edm::InputTag>("gpuOn"))),
        gpuOffToken_(consumes(config.getParameter<edm::InputTag>("gpuOff"))),
        gpuRecoveredToken_(consumes(config.getParameter<edm::InputTag>("gpuRecovered"))),
        geomToken_(esConsumes()),
        topoToken_(esConsumes()),
        maxDisplays_(config.getParameter<int>("maxDisplays")) {
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    modules_ = fs->make<TTree>("modules", "one entry per module and event");
    modules_->Branch("event", &event_);
    modules_->Branch("detId", &detId_);
    modules_->Branch("barrel", &barrel_);
    modules_->Branch("layer", &layer_);
    modules_->Branch("rows", &rows_);
    modules_->Branch("cols", &cols_);
    modules_->Branch("nReal", &nReal_);
    modules_->Branch("nFakeCPU", &nFakeCPU_);
    modules_->Branch("nFakeGPU", &nFakeGPU_);
    modules_->Branch("nFakeBoth", &nFakeBoth_);
    modules_->Branch("nFakeRef", &nFakeRef_);
    modules_->Branch("nFakeRefZero", &nFakeRefZero_);
    modules_->Branch("nDiffGPURef", &nDiffGPURef_);
    modules_->Branch("nDiffCPURefZero", &nDiffCPURefZero_);
    modules_->Branch("nClusOffCPU", &nClusOffCPU_);
    modules_->Branch("nClusOffGPU", &nClusOffGPU_);
    modules_->Branch("nClusOnCPU", &nClusOnCPU_);
    modules_->Branch("nClusOnGPU", &nClusOnGPU_);
    modules_->Branch("nMismatchOff", &nMismatchOff_);
    modules_->Branch("nMismatchOn", &nMismatchOn_);
    modules_->Branch("nInvalidGPU", &nInvalidGPU_);

    clusters_ = fs->make<TTree>("clusters", "device clusters in the barrel layer 1, with and without morphing");
    clusters_->Branch("event", &event_);
    clusters_->Branch("detId", &detId_);
    clusters_->Branch("morphing", &clusMorphing_);
    clusters_->Branch("npix", &clusNpix_);
    clusters_->Branch("sizeX", &clusSizeX_);
    clusters_->Branch("sizeY", &clusSizeY_);

    fakeDiff_ = fs->make<TTree>("fakeDiff", "recovered pixels found by only one of CPU and GPU");
    fakeDiff_->Branch("event", &event_);
    fakeDiff_->Branch("detId", &detId_);
    fakeDiff_->Branch("rows", &rows_);
    fakeDiff_->Branch("cols", &cols_);
    fakeDiff_->Branch("x", &diffX_);
    fakeDiff_->Branch("y", &diffY_);
    fakeDiff_->Branch("onlyGPU", &diffOnlyGPU_);
    fakeDiff_->Branch("inRef", &diffInRef_);
    fakeDiff_->Branch("inRefZero", &diffInRefZero_);

    displays_ = fs->make<TTree>("displays", "pixel maps of modules where the morphing merged clusters");
    displays_->Branch("event", &event_);
    displays_->Branch("detId", &detId_);
    displays_->Branch("rows", &rows_);
    displays_->Branch("cols", &cols_);
    displays_->Branch("x", &dispX_);
    displays_->Branch("y", &dispY_);
    displays_->Branch("off", &dispOff_);
    displays_->Branch("on", &dispOn_);
    displays_->Branch("gpuFakeX", &dispGpuFakeX_);
    displays_->Branch("gpuFakeY", &dispGpuFakeY_);
    displays_->Branch("cpuFakeX", &dispCpuFakeX_);
    displays_->Branch("cpuFakeY", &dispCpuFakeY_);
  }

  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override {
    auto const& digis = iEvent.get(digisToken_);
    auto const& cpuMorphed = iEvent.get(cpuMorphedToken_);
    auto const gpuOn = iEvent.get(gpuOnToken_).const_view();
    auto const gpuOff = iEvent.get(gpuOffToken_).const_view();
    auto const recovered = iEvent.get(gpuRecoveredToken_).const_view();
    auto const& geom = iSetup.getData(geomToken_);
    auto const& topo = iSetup.getData(topoToken_);
    event_ = iEvent.id().event();

    // pixels recovered on the device, by module index
    std::unordered_map<int, std::vector<std::pair<int, int>>> gpuFakes;
    for (int i = 0; i < recovered.metadata().size(); ++i)
      if (recovered[i].moduleId() != std::numeric_limits<uint16_t>::max())
        gpuFakes[recovered[i].moduleId()].emplace_back(recovered[i].xx(), recovered[i].yy());

    // the device digis follow the order of the input DetSetVector
    int offset = 0;
    for (auto const& det : digis) {
      const int n = det.size();
      const int first = offset;
      offset += n;
      if (n == 0)
        continue;
      DetId id(det.detId());
      auto const* unit = dynamic_cast<PixelGeomDetUnit const*>(geom.idToDetUnit(id));
      detId_ = id.rawId();
      barrel_ = id.subdetId() == PixelSubdetector::PixelBarrel;
      layer_ = topo.layer(id);
      rows_ = unit->specificTopology().nrows();
      cols_ = unit->specificTopology().ncolumns();
      const int moduleIndex = unit->index();

      std::vector<std::pair<int, int>> real;
      std::vector<int> off, on;
      nInvalidGPU_ = 0;
      for (int k = 0; k < n; ++k) {
        int i = first + k;
        if (gpuOn[i].moduleId() != moduleIndex or gpuOn[i].xx() != det.data[k].row() or
            gpuOn[i].yy() != det.data[k].column())
          throw cms::Exception("Validation") << "device digis are not in the order of the input digis";
        real.emplace_back(det.data[k].row(), det.data[k].column());
        off.push_back(gpuOff[i].clus());
        on.push_back(gpuOn[i].clus());
        nInvalidGPU_ += (gpuOff[i].clus() < 0) + (gpuOn[i].clus() < 0);
      }
      nReal_ = n;

      // pixels added by the CPU morphing
      std::set<Key> realKeys;
      for (auto const& [x, y] : real)
        realKeys.insert(key(x, y));
      std::vector<std::pair<int, int>> cpuFakes;
      auto morphed = cpuMorphed.find(det.detId());
      if (morphed != cpuMorphed.end())
        for (auto const& p : *morphed)
          if (not realKeys.count(key(p.row(), p.column())))
            cpuFakes.emplace_back(p.row(), p.column());
      auto const& deviceFakes = gpuFakes[moduleIndex];
      std::set<Key> cpuKeys;
      for (auto const& [x, y] : cpuFakes)
        cpuKeys.insert(key(x, y));
      nFakeCPU_ = cpuFakes.size();
      nFakeGPU_ = deviceFakes.size();
      std::set<Key> ref, refZero;
      if (barrel_ and layer_ == 1) {
        ref = exactMorphing(real, rows_, cols_, true);
        refZero = exactMorphing(real, rows_, cols_, false);
      }
      nFakeRef_ = ref.size();
      nFakeRefZero_ = refZero.size();
      nFakeBoth_ = 0;
      std::set<Key> gpuKeys;
      for (auto const& [x, y] : deviceFakes) {
        gpuKeys.insert(key(x, y));
        bool both = cpuKeys.count(key(x, y));
        nFakeBoth_ += both;
        if (not both) {
          diffX_ = x;
          diffY_ = y;
          diffOnlyGPU_ = 1;
          diffInRef_ = ref.count(key(x, y));
          diffInRefZero_ = refZero.count(key(x, y));
          fakeDiff_->Fill();
        }
      }
      for (auto const& [x, y] : cpuFakes)
        if (not gpuKeys.count(key(x, y))) {
          diffX_ = x;
          diffY_ = y;
          diffOnlyGPU_ = 0;
          diffInRef_ = ref.count(key(x, y));
          diffInRefZero_ = refZero.count(key(x, y));
          fakeDiff_->Fill();
        }
      // symmetric differences: GPU vs exact reference, CPU vs exact zero-padded reference
      auto symdiff = [](std::set<Key> const& a, std::set<Key> const& b) {
        int n = 0;
        for (auto k : a)
          n += not b.count(k);
        for (auto k : b)
          n += not a.count(k);
        return n;
      };
      nDiffGPURef_ = symdiff(gpuKeys, ref);
      nDiffCPURefZero_ = symdiff(cpuKeys, refZero);

      // reference clusterings: the input digis, and the input plus the CPU recovered digis
      auto refOff = connected(real);
      auto all = real;
      all.insert(all.end(), cpuFakes.begin(), cpuFakes.end());
      auto refOn = connected(all);
      refOn.resize(n);

      nClusOffCPU_ = distinct(refOff);
      nClusOffGPU_ = distinct(off);
      nClusOnCPU_ = distinct(refOn);
      nClusOnGPU_ = distinct(on);
      nMismatchOff_ = mismatches(off, refOff);
      nMismatchOn_ = mismatches(on, refOn);
      modules_->Fill();

      if (barrel_ and layer_ == 1) {
        for (int morphing : {0, 1}) {
          auto const& labels = morphing ? on : off;
          std::map<int, std::vector<int>> members;
          for (int k = 0; k < n; ++k)
            members[labels[k]].push_back(k);
          for (auto const& [label, pix] : members) {
            int x0 = 1 << 30, x1 = -1, y0 = 1 << 30, y1 = -1;
            for (int k : pix) {
              x0 = std::min(x0, real[k].first);
              x1 = std::max(x1, real[k].first);
              y0 = std::min(y0, real[k].second);
              y1 = std::max(y1, real[k].second);
            }
            clusMorphing_ = morphing;
            clusNpix_ = pix.size();
            clusSizeX_ = x1 - x0 + 1;
            clusSizeY_ = y1 - y0 + 1;
            clusters_->Fill();
          }
        }
        if (displays_->GetEntries() < maxDisplays_ and nClusOnGPU_ + 3 <= nClusOffGPU_) {
          dispX_.clear();
          dispY_.clear();
          for (auto const& [x, y] : real) {
            dispX_.push_back(x);
            dispY_.push_back(y);
          }
          dispOff_ = off;
          dispOn_ = on;
          dispGpuFakeX_.clear();
          dispGpuFakeY_.clear();
          for (auto const& [x, y] : deviceFakes) {
            dispGpuFakeX_.push_back(x);
            dispGpuFakeY_.push_back(y);
          }
          dispCpuFakeX_.clear();
          dispCpuFakeY_.clear();
          for (auto const& [x, y] : cpuFakes) {
            dispCpuFakeX_.push_back(x);
            dispCpuFakeY_.push_back(y);
          }
          displays_->Fill();
        }
      }
    }
  }

private:
  const edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> digisToken_;
  const edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> cpuMorphedToken_;
  const edm::EDGetTokenT<SiPixelDigisHost> gpuOnToken_;
  const edm::EDGetTokenT<SiPixelDigisHost> gpuOffToken_;
  const edm::EDGetTokenT<SiPixelDigisHost> gpuRecoveredToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const int maxDisplays_;

  TTree* modules_;
  TTree* clusters_;
  TTree* displays_;
  TTree* fakeDiff_;
  int diffX_, diffY_, diffOnlyGPU_, diffInRef_, diffInRefZero_;
  int nFakeRef_, nFakeRefZero_, nDiffGPURef_, nDiffCPURefZero_;
  uint32_t event_, detId_;
  bool barrel_;
  int layer_, rows_, cols_, nReal_, nFakeCPU_, nFakeGPU_, nFakeBoth_;
  int nClusOffCPU_, nClusOffGPU_, nClusOnCPU_, nClusOnGPU_, nMismatchOff_, nMismatchOn_, nInvalidGPU_;
  int clusMorphing_, clusNpix_, clusSizeX_, clusSizeY_;
  std::vector<int> dispX_, dispY_, dispOff_, dispOn_, dispGpuFakeX_, dispGpuFakeY_, dispCpuFakeX_, dispCpuFakeY_;
};

class Phase2ClusDump : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit Phase2ClusDump(edm::ParameterSet const& config)
      : digisToken_(consumes(config.getParameter<edm::InputTag>("gpuDigis"))) {
    usesResource("TFileService");
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("dump", "device cluster index and charge of every digi");
    tree_->Branch("event", &event_);
    tree_->Branch("clus", &clus_);
    tree_->Branch("adc", &adc_);
    tree_->Branch("moduleId", &moduleId_);
    tree_->Branch("x", &x_);
    tree_->Branch("y", &y_);
  }

  void analyze(edm::Event const& iEvent, edm::EventSetup const&) override {
    auto const view = iEvent.get(digisToken_).const_view();
    event_ = iEvent.id().event();
    clus_.clear();
    adc_.clear();
    moduleId_.clear();
    x_.clear();
    y_.clear();
    for (int i = 0; i < view.metadata().size(); ++i) {
      x_.push_back(view[i].xx());
      y_.push_back(view[i].yy());
      clus_.push_back(view[i].clus());
      adc_.push_back(view[i].adc());
      moduleId_.push_back(view[i].moduleId());
    }
    tree_->Fill();
  }

private:
  const edm::EDGetTokenT<SiPixelDigisHost> digisToken_;
  TTree* tree_;
  uint32_t event_;
  std::vector<int> clus_, adc_, moduleId_, x_, y_;
};

DEFINE_FWK_MODULE(Phase2MorphValidation);
DEFINE_FWK_MODULE(Phase2ClusDump);
