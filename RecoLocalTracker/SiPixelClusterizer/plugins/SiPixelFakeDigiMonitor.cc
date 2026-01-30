#include "RecoLocalTracker/SiPixelClusterizer/plugins/SiPixelFakeDigiMonitor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/PixelBarrelName.h"
#include "DataFormats/TrackerCommon/interface/PixelEndcapName.h"

SiPixelFakeDigiMonitor::SiPixelFakeDigiMonitor(const edm::ParameterSet& iConfig)
    : clustersSoASrc_(iConfig.getParameter<edm::InputTag>("clustersSoASrc")),
      digisSoASrc_(iConfig.getParameter<edm::InputTag>("digisSoASrc")),
      digiSrc_(iConfig.getParameter<edm::InputTag>("digiSrc")),
      topFolderName_(iConfig.getParameter<std::string>("topFolderName")),
      clustersSoAToken_(consumes<SiPixelClustersHost>(clustersSoASrc_)),
      digisSoAToken_(consumes<SiPixelDigisHost>(digisSoASrc_)),
      digiToken_(consumes<edm::DetSetVector<PixelDigi>>(digiSrc_)),
      trackerTopoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>()),
      trackerGeomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()),
      fakeDigisThisLS_(0),
      duplicatePixelsThisLS_(0),
      totalDigisThisLS_(0),
      eventsThisLS_(0) {
  edm::LogInfo("SiPixelFakeDigiMonitor") << "Monitoring fake digis from clusters SoA: " << clustersSoASrc_;
}

SiPixelFakeDigiMonitor::~SiPixelFakeDigiMonitor() {}

void SiPixelFakeDigiMonitor::bookHistograms(DQMStore::IBooker& iBooker,
                                            edm::Run const& iRun,
                                            edm::EventSetup const& iSetup) {
  
  // Constants for per-layer occupancy binning
  const int nModulesPerLadder = 8;  // Modules along z per ladder
  const int nLaddersPerLayer[4] = {12, 28, 44, 64};  // Ladders per layer

  iBooker.setCurrentFolder(topFolderName_ + "/FakeDigis");

  // Event-level histograms
  meFakeDigisPerEvent_ = iBooker.book1D("FakeDigisPerEvent",
                                        "Number of Fake Digis per Event;Fake Digis;Events",
                                        500, 0., 8000.);

  // Module-level histograms
  meFakeDigisPerModule_ = iBooker.book1D("FakeDigisPerModule",
                                         "Number of Fake Digis per Module;Fake Digis;Modules (RawID)",
                                         2000, 0., 8000.);

  meFakeDigisPerModuleBarrel_ = iBooker.book1D("FakeDigisPerModuleBarrel",
                                               "Fake Digis per Module (Barrel);Fake Digis;Modules",
                                               100, 0., 8000.);

  meModulesWithFakeDigis_ = iBooker.book1D("ModulesWithFakeDigis",
                                           "Number of Modules with Fake Digis;Modules;Events",
                                           100, 0., 2000.);

  meNFakeDigisBarrel_ = iBooker.book1D("NFakeDigisBarrel",
                                       "Number of Fake Digis (Barrel);Fake Digis;Events",
                                       500, 0., 8000.);

  // Per-layer histograms
  iBooker.setCurrentFolder(topFolderName_ + "/FakeDigis/Barrel");
  for (int i = 0; i < 4; i++) {
    std::string hname = "NFakeDigisLayer" + std::to_string(i + 1);
    std::string htitle = "Number of Fake Digis (Layer " + std::to_string(i + 1) + ");Fake Digis;Events";
    meNFakeDigisLayer_[i] = iBooker.book1D(hname, htitle, 200, 0., 8000.);
  }

  // Per-disk histograms
  iBooker.setCurrentFolder(topFolderName_ + "/FakeDigis/Endcap");
  for (int i = 0; i < 3; i++) {
    std::string hname = "NFakeDigisDisk" + std::to_string(i + 1);
    std::string htitle = "Number of Fake Digis (Disk " + std::to_string(i + 1) + " +/-);Fake Digis;Events";
    meNFakeDigisDisk_[i] = iBooker.book1D(hname, htitle, 200, 0., 8000.);
  }

  // Per-layer occupancy maps: Module (x) vs Ladder (y)
  iBooker.setCurrentFolder(topFolderName_ + "/FakeDigis/Occupancy");
  for (int i = 0; i < 4; i++) {
    int nLadders = nLaddersPerLayer[i];
    int halfLadders = nLadders / 2;
    std::string hname = "FakeDigiOccupancyBarrelLayer" + std::to_string(i + 1);
    std::string htitle = "Fake Digi Occupancy (Barrel Layer " + std::to_string(i + 1) +
                         ");Module;Ladder";
    meFakeDigiOccupancyBarrelLayer_[i] = iBooker.book2D(hname, htitle,
                                                        nModulesPerLadder + 1, -4.5, 4.5,
                                                        nLadders + 1, -(halfLadders + 0.5), halfLadders + 0.5);
    // Set COLZ draw option, keep 0-bins white, and add bin labels for clear module borders
    meFakeDigiOccupancyBarrelLayer_[i]->setOption("colz");
    meFakeDigiOccupancyBarrelLayer_[i]->getTH2F()->SetMinimum(0.001);
    for (int j = 1; j <= nModulesPerLadder + 1; j++) {
      meFakeDigiOccupancyBarrelLayer_[i]->setBinLabel(j, std::to_string(j - 5), 1);
    }
  }

  // HLT monitoring
  iBooker.setCurrentFolder(topFolderName_ + "/FakeDigis/HLT");
  meFakeDigisPerEventVsLS_ = iBooker.bookProfile("FakeDigisPerEventVsLS", "Fake Digis per Event vs LS;Lumisection;Fake Digis/Event", 2000, 0., 2000., 0., 8000.);
  meFakeDigisFractionVsLS_ = iBooker.bookProfile("FakeDigisFractionVsLS", "Fake Digis Fraction vs LS;Lumisection;Fake Fraction", 2000, 0., 2000., 0., 1.0);

  // ============ Duplicate Pixels Histograms ============
  iBooker.setCurrentFolder(topFolderName_ + "/DuplicatePixels");
  meDuplicatePixelsPerEvent_ = iBooker.book1D("DuplicatePixelsPerEvent", "Number of Duplicate Pixels per Event;Duplicate Pixels;Events", 500, 0., 500.);
  meDuplicatePixelsPerModule_ = iBooker.book1D("DuplicatePixelsPerModule", "Number of Duplicate Pixels per Module;Duplicate Pixels;Modules", 200, 0., 500.);
  meDuplicatePixelsPerModuleBarrel_ = iBooker.book1D("DuplicatePixelsPerModuleBarrel", "Duplicate Pixels per Module (Barrel);Duplicate Pixels;Modules", 200, 0., 500.);
  meDuplicatePixelsPerModuleEndcap_ = iBooker.book1D("DuplicatePixelsPerModuleEndcap", "Duplicate Pixels per Module (Endcap);Duplicate Pixels;Modules", 200, 0., 500.);
  meModulesWithDuplicatePixels_ = iBooker.book1D("ModulesWithDuplicatePixels", "Number of Modules with Duplicate Pixels;Modules;Events", 100, 0., 100.);
  meNDuplicatePixelsBarrel_ = iBooker.book1D("NDuplicatePixelsBarrel", "Number of Duplicate Pixels (Barrel);Duplicate Pixels;Events", 500, 0., 500.);
  meNDuplicatePixelsEndcap_ = iBooker.book1D("NDuplicatePixelsEndcap", "Number of Duplicate Pixels (Endcap);Duplicate Pixels;Events", 500, 0., 500.);

  iBooker.setCurrentFolder(topFolderName_ + "/DuplicatePixels/Barrel");
  for (int i = 0; i < 4; i++) {
    meNDuplicatePixelsLayer_[i] = iBooker.book1D("NDuplicatePixelsLayer" + std::to_string(i + 1), "Number of Duplicate Pixels (Layer " + std::to_string(i + 1) + ");Duplicate Pixels;Events", 200, 0., 500.);
  }

  iBooker.setCurrentFolder(topFolderName_ + "/DuplicatePixels/Endcap");
  for (int i = 0; i < 3; i++) {
    meNDuplicatePixelsDisk_[i] = iBooker.book1D("NDuplicatePixelsDisk" + std::to_string(i + 1), "Number of Duplicate Pixels (Disk " + std::to_string(i + 1) + " +/-);Duplicate Pixels;Events", 200, 0., 500.);
  }

  // Per-layer occupancy maps: Module (x) vs Ladder (y)
  iBooker.setCurrentFolder(topFolderName_ + "/DuplicatePixels/Occupancy");
  for (int i = 0; i < 4; i++) {
    int nLadders = nLaddersPerLayer[i];
    int halfLadders = nLadders / 2;
    std::string hname = "DuplicatePixelOccupancyBarrelLayer" + std::to_string(i + 1);
    std::string htitle = "Duplicate Pixel Occupancy (Barrel Layer " + std::to_string(i + 1) +
                         ");Module;Ladder";
    meDuplicatePixelOccupancyBarrelLayer_[i] = iBooker.book2D(hname, htitle,
                                                              nModulesPerLadder + 1, -4.5, 4.5,
                                                              nLadders + 1, -(halfLadders + 0.5), halfLadders + 0.5);
    // Set COLZ draw option, keep 0-bins white, and add bin labels for clear module borders
    meDuplicatePixelOccupancyBarrelLayer_[i]->setOption("colz");
    meDuplicatePixelOccupancyBarrelLayer_[i]->getTH2F()->SetMinimum(0.001);
    for (int j = 1; j <= nModulesPerLadder + 1; j++) {
      meDuplicatePixelOccupancyBarrelLayer_[i]->setBinLabel(j, std::to_string(j - 5), 1);
    }
  }

  iBooker.setCurrentFolder(topFolderName_ + "/DuplicatePixels/HLT");
  meDuplicatePixelsPerEventVsLS_ = iBooker.bookProfile("DuplicatePixelsPerEventVsLS", "Duplicate Pixels per Event vs LS;Lumisection;Duplicate Pixels/Event", 2000, 0., 2000., 0., 500.);
  meDuplicatePixelsFractionVsLS_ = iBooker.bookProfile("DuplicatePixelsFractionVsLS", "Duplicate Pixels Fraction vs LS;Lumisection;Duplicate Fraction", 2000, 0., 2000., 0., 1.0);

  edm::LogInfo("SiPixelFakeDigiMonitor") << "Histograms booked with per-layer Module vs Ladder occupancy maps.";
}

void SiPixelFakeDigiMonitor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const TrackerTopology& tTopo = iSetup.getData(trackerTopoToken_);

  edm::Handle<SiPixelClustersHost> clustersSoA;
  iEvent.getByToken(clustersSoAToken_, clustersSoA);

  if (!clustersSoA.isValid()) return;

  edm::Handle<edm::DetSetVector<PixelDigi>> digis;
  iEvent.getByToken(digiToken_, digis);

  int totalDigis = 0;
  if (digis.isValid()) {
    for (const auto& detSet : *digis) totalDigis += detSet.size();
  }

  auto const& clusters_view = clustersSoA->const_view();
  const uint32_t nModules = clusters_view[0].moduleStart();

  int totalFakeDigis = 0;
  int fakeDigisBarrel = 0;
  std::vector<int> fakeDigisPerLayer(4, 0);
  std::vector<int> fakeDigisPerDisk(3, 0);
  int modulesWithFakeDigis = 0;

  int totalDuplicatePixels = 0;
  int duplicatePixelsBarrel = 0;
  int duplicatePixelsEndcap = 0;
  std::vector<int> duplicatePixelsPerLayer(4, 0);
  std::vector<int> duplicatePixelsPerDisk(3, 0);
  int modulesWithDuplicatePixels = 0;

  for (uint32_t i = 1; i <= nModules; ++i) {
    uint32_t nFakeDigis = clusters_view[i].nFakeDigis();
    uint32_t nDuplicatePixels = clusters_view[i].nDuplicatePixels();

    if (nFakeDigis == 0 && nDuplicatePixels == 0) continue;

    uint32_t rawId = clusters_view[i].rawId();
    DetId detIdObj(rawId);
    if (detIdObj.det() != DetId::Tracker) continue;

    int subdet = detIdObj.subdetId();

    if (nFakeDigis > 0) {
      totalFakeDigis += nFakeDigis;
      modulesWithFakeDigis++;
      meFakeDigisPerModule_->Fill(nFakeDigis);
      if (subdet == static_cast<int>(PixelSubdetector::PixelBarrel)) {
        meFakeDigisPerModuleBarrel_->Fill(nFakeDigis);
        fakeDigisBarrel += nFakeDigis;
        int layer = tTopo.pxbLayer(detIdObj);
        if (layer >= 1 && layer <= 4) {
          fakeDigisPerLayer[layer - 1] += nFakeDigis;
          unsigned int module = tTopo.pxbModule(detIdObj);
          unsigned int ladder = tTopo.pxbLadder(detIdObj);
          const unsigned int nLaddersPerLayer[4] = {12, 28, 44, 64};
          unsigned int halfLadders = nLaddersPerLayer[layer - 1] / 2;
          int signedLadder = (ladder <= halfLadders)
                                 ? static_cast<int>(ladder) - static_cast<int>(halfLadders) - 1
                                 : static_cast<int>(ladder) - static_cast<int>(halfLadders);
          int signedModule = (module <= 4u) ? static_cast<int>(module) - 5 : static_cast<int>(module) - 4;
          meFakeDigiOccupancyBarrelLayer_[layer - 1]->Fill(signedModule, signedLadder);
        }
      } else if (subdet == static_cast<int>(PixelSubdetector::PixelEndcap)) {
        int disk = tTopo.pxfDisk(detIdObj);
        if (disk >= 1 && disk <= 3) fakeDigisPerDisk[disk - 1] += nFakeDigis;
      }
    }

    if (nDuplicatePixels > 0) {
      totalDuplicatePixels += nDuplicatePixels;
      modulesWithDuplicatePixels++;
      meDuplicatePixelsPerModule_->Fill(nDuplicatePixels);
      if (subdet == static_cast<int>(PixelSubdetector::PixelBarrel)) {
        meDuplicatePixelsPerModuleBarrel_->Fill(nDuplicatePixels);
        duplicatePixelsBarrel += nDuplicatePixels;
        int layer = tTopo.pxbLayer(detIdObj);
        if (layer >= 1 && layer <= 4) {
          duplicatePixelsPerLayer[layer - 1] += nDuplicatePixels;
          unsigned int module = tTopo.pxbModule(detIdObj);
          unsigned int ladder = tTopo.pxbLadder(detIdObj);
          const unsigned int nLaddersPerLayer[4] = {12, 28, 44, 64};
          unsigned int halfLadders = nLaddersPerLayer[layer - 1] / 2;
          int signedLadder = (ladder <= halfLadders)
                                 ? static_cast<int>(ladder) - static_cast<int>(halfLadders) - 1
                                 : static_cast<int>(ladder) - static_cast<int>(halfLadders);
          int signedModule = (module <= 4u) ? static_cast<int>(module) - 5 : static_cast<int>(module) - 4;
          meDuplicatePixelOccupancyBarrelLayer_[layer - 1]->Fill(signedModule, signedLadder);
        }
      } else if (subdet == static_cast<int>(PixelSubdetector::PixelEndcap)) {
        meDuplicatePixelsPerModuleEndcap_->Fill(nDuplicatePixels);
        duplicatePixelsEndcap += nDuplicatePixels;
        int disk = tTopo.pxfDisk(detIdObj);
        if (disk >= 1 && disk <= 3) duplicatePixelsPerDisk[disk - 1] += nDuplicatePixels;
      }
    }
  }

  // Event-level fills
  meFakeDigisPerEvent_->Fill(totalFakeDigis);
  meModulesWithFakeDigis_->Fill(modulesWithFakeDigis);
  meNFakeDigisBarrel_->Fill(fakeDigisBarrel);
  meDuplicatePixelsPerEvent_->Fill(totalDuplicatePixels);
  meModulesWithDuplicatePixels_->Fill(modulesWithDuplicatePixels);
  meNDuplicatePixelsBarrel_->Fill(duplicatePixelsBarrel);
  meNDuplicatePixelsEndcap_->Fill(duplicatePixelsEndcap);

  // LS Accumulation
  fakeDigisThisLS_ += totalFakeDigis;
  duplicatePixelsThisLS_ += totalDuplicatePixels;
  totalDigisThisLS_ += totalDigis;
  eventsThisLS_++;

  for (int i = 0; i < 4; i++) {
    if (fakeDigisPerLayer[i] > 0) meNFakeDigisLayer_[i]->Fill(fakeDigisPerLayer[i]);
    if (duplicatePixelsPerLayer[i] > 0) meNDuplicatePixelsLayer_[i]->Fill(duplicatePixelsPerLayer[i]);
  }
  for (int i = 0; i < 3; i++) {
    if (fakeDigisPerDisk[i] > 0) meNFakeDigisDisk_[i]->Fill(fakeDigisPerDisk[i]);
    if (duplicatePixelsPerDisk[i] > 0) meNDuplicatePixelsDisk_[i]->Fill(duplicatePixelsPerDisk[i]);
  }
}

DEFINE_FWK_MODULE(SiPixelFakeDigiMonitor);
