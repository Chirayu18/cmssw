import FWCore.ParameterSet.Config as cms

# Import the fake digi monitor configured for HLT
from RecoLocalTracker.SiPixelClusterizer.SiPixelFakeDigiMonitor_cfi import SiPixelFakeDigiMonitorHLT

# Configure for HLT environment with correct input tags
# Uses SoA clusters which contain the nFakeDigis column
hltPixelFakeDigiMonitor = SiPixelFakeDigiMonitorHLT.clone(
    topFolderName = cms.string("HLT/Pixel/FakeDigis"),
    clustersSoASrc = cms.InputTag("hltSiPixelClustersAlpaka"),  # SoA clusters from HLT
    digisSoASrc = cms.InputTag("hltSiPixelDigisSoAAlpaka"),  # SoA digis for moduleId->DetId mapping
)

# Sequence for HLT DQM
hltPixelFakeDigiMonitoringSequence = cms.Sequence(
    hltPixelFakeDigiMonitor
)
