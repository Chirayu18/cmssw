#include <cstdint>
#include <memory>
#include <vector>
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"  // Include pixel cluster data format
#include "DataFormats/SiPixelClusterSoA/interface/SiPixelClustersSoA.h"  // Include SoA format for pixel clusters
#include "DataFormats/SiPixelClusterSoA/interface/alpaka/SiPixelClustersSoACollection.h" // Include Alpaka SoA format for pixel clusters
#include "DataFormats/Common/interface/DetSetVectorNew.h"  // For handling collections of pixel clusters
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"  // Alpaka event handling
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"  // Alpaka configuration
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"  // Define a stream module
#include "FWCore/ParameterSet/interface/ParameterSet.h"  // Handle configuration parameters
							 //
namespace ALPAKA_ACCELERATOR_NAMESPACE {

class SiPixelClusterSoAConverter : public stream::EDProducer<> {

  
public:
  explicit SiPixelClusterSoAConverter(const edm::ParameterSet& iConfig);
  ~SiPixelClusterSoAConverter() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(device::Event& iEvent, const device::EventSetup& es) override;

  const edm::EDGetTokenT<DetSetVector<SiPixelCluster>> clusterToken_; 
  const edm::EDPutToken<SiPixelClustersSoACollection> clusterSoAToken_;
};


SiPixelClusterSoAConverter::SiPixelClusterSoAConverter(const edm::ParameterSet& iConfig)
    : clusterToken_(consumes<edm::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("clusterSource"))),
      clusterSoAToken_(produces<edm::SiPixelClustersSoACollection>()) {}

void SiPixelClusterSoAConverter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("clusterSource", edm::InputTag("siPixelClusters"));
  descriptions.addWithDefaultLabel(desc);

}

void SiPixelClusterSoAConverter::produce(device::Event& iEvent, device::EventSetup const& iSetup) {
  const auto& pixelClusters = iEvent.get(clusterToken_);
  auto queue = iEvent.queue();
  int totalClusters = 0;

  for (const auto& detSet : pixelClusters) {
    totalClusters += detSet.size();
  }
  cout<< "Total Clusters: "<< totalClusters<<endl;

  SiPixelClustersSoACollection clustersSoA(queue, totalClusters);

  size_t clusterIndex = 0;

  for (const auto& detSet : pixelClusters) {
    for (const auto& cluster : detSet) {
      
      // Copy the cluster properties to the SoA collection
      clustersSoA.view()[clusterIndex].x() = cluster.x();
      clustersSoA.view()[clusterIndex].y() = cluster.y();
      clustersSoA.view()[clusterIndex].charge() = cluster.charge();
      clustersSoA.view()[clusterIndex].sizeX() = cluster.sizeX();
      clustersSoA.view()[clusterIndex].sizeY() = cluster.sizeY();

      // Increment cluster index
      ++clusterIndex;
    }
  }

  iEvent.emplace(clusterSoAToken_, std::move(clustersSoA)); 
}

}

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(SiPixelClusterSoAConverter);
