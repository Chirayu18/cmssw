import FWCore.ParameterSet.Config as cms

# Configuration for HLT monitoring of fake digis
from RecoLocalTracker.SiPixelClusterizer.SiPixelFakeDigiMonitor_cfi import SiPixelFakeDigiMonitorHLT

# Sequence for HLT DQM
siPixelFakeDigiMonitoringHLT = cms.Sequence(
    SiPixelFakeDigiMonitorHLT
)

# For integration into HLT DQM sequences
# Add to HLT DQM configuration:
# process.load("RecoLocalTracker.SiPixelClusterizer.SiPixelFakeDigiMonitor_HLT_cff")
# process.HLTPixelMonitoring += process.siPixelFakeDigiMonitoringHLT
