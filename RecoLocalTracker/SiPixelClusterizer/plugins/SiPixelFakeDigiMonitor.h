#ifndef RecoLocalTracker_SiPixelClusterizer_SiPixelFakeDigiMonitor_h
#define RecoLocalTracker_SiPixelClusterizer_SiPixelFakeDigiMonitor_h

/** \class SiPixelFakeDigiMonitor
 * DQM Monitoring for fake/recovered digis from pixel clustering
 * 
 * This module monitors the number and properties of fake digis
 * that are recovered during the pixel clustering process through
 * digi morphing. These are pixels that are filled in to maintain
 * cluster contiguity.
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelClusterSoA/interface/SiPixelClustersHost.h"
#include "DataFormats/SiPixelDigiSoA/interface/SiPixelDigisHost.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

class SiPixelFakeDigiMonitor : public DQMEDAnalyzer {
public:
  explicit SiPixelFakeDigiMonitor(const edm::ParameterSet&);
  ~SiPixelFakeDigiMonitor() override;

protected:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  // Parameters
  edm::InputTag clustersSoASrc_;
  edm::InputTag digisSoASrc_;
  edm::InputTag digiSrc_;
  std::string topFolderName_;
  
  // Tokens
  edm::EDGetTokenT<SiPixelClustersHost> clustersSoAToken_;
  edm::EDGetTokenT<SiPixelDigisHost> digisSoAToken_;
  edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> digiToken_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopoToken_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeomToken_;

  // Histograms - Event level
  MonitorElement* meFakeDigisPerEvent_;
  MonitorElement* meFakeDigisFraction_;
  MonitorElement* meFakeDigisPerModule_;
  MonitorElement* meFakeDigisPerModuleBarrel_;
  MonitorElement* meFakeDigisPerModuleEndcap_;
  
  // HLT monitoring - per lumisection
  MonitorElement* meFakeDigisPerEventVsLS_;
  MonitorElement* meFakeDigisFractionVsLS_;
  MonitorElement* meFakeDigisRate_;
  
  MonitorElement* meFakeDigisPerCluster_;
  MonitorElement* meFakeDigisPerClusterBarrel_;
  MonitorElement* meFakeDigisPerClusterEndcap_;
  
  MonitorElement* meNFakeDigisBarrel_;
  MonitorElement* meNFakeDigisEndcap_;
  MonitorElement* meNFakeDigisLayer_[4];  // Barrel layers
  MonitorElement* meNFakeDigisDisk_[3];   // Forward disks (both +/-)
  
  MonitorElement* meFakeDigisPerDetIdLayer_[4];  // Profile per layer: DetId vs fake digis
  MonitorElement* meFakeDigisPerDetIdDisk_[3];   // Profile per disk: DetId vs fake digis
  
  MonitorElement* meFakeDigiDistribution_;
  MonitorElement* meFakeDigiOccupancyBarrelLayer_[4];  // Per-layer: Module vs Ladder
  MonitorElement* meFakeDigiOccupancyEndcap_;
  MonitorElement* meFakeDigisPerDetId_;  // Profile: avg fake digis per module DetId
  
  // Statistics
  MonitorElement* meTotalDigisPerEvent_;
  MonitorElement* meModulesWithFakeDigis_;

  // Duplicate Pixels Histograms - Event level
  MonitorElement* meDuplicatePixelsPerEvent_;
  MonitorElement* meDuplicatePixelsPerModule_;
  MonitorElement* meDuplicatePixelsPerModuleBarrel_;
  MonitorElement* meDuplicatePixelsPerModuleEndcap_;
  MonitorElement* meModulesWithDuplicatePixels_;

  // Duplicate Pixels - Barrel/Endcap summary
  MonitorElement* meNDuplicatePixelsBarrel_;
  MonitorElement* meNDuplicatePixelsEndcap_;
  MonitorElement* meNDuplicatePixelsLayer_[4];  // Barrel layers
  MonitorElement* meNDuplicatePixelsDisk_[3];   // Forward disks (both +/-)

  // Duplicate Pixels - Occupancy
  MonitorElement* meDuplicatePixelOccupancyBarrelLayer_[4];  // Per-layer: Module vs Ladder

  // Duplicate Pixels - HLT monitoring
  MonitorElement* meDuplicatePixelsPerEventVsLS_;
  MonitorElement* meDuplicatePixelsFractionVsLS_;
  
  // Lumisection counters for HLT monitoring
  std::atomic<unsigned long long> fakeDigisThisLS_;
  std::atomic<unsigned long long> duplicatePixelsThisLS_;
  std::atomic<unsigned long long> totalDigisThisLS_;
  std::atomic<unsigned long long> eventsThisLS_;
  
  // Static moduleId -> DetId mapping (built once, reused for all events)
  std::unordered_map<uint32_t, uint32_t> moduleIdToDetId_;
};

#endif
