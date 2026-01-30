import FWCore.ParameterSet.Config as cms

from DQMOffline.Trigger.SiPixel_OfflineMonitoring_cff import *
from DQM.HLTEvF.HLTPixelFakeDigiMonitoring_cff import *

pixelOnlineMonitorHLTsequence = cms.Sequence(
    sipixelMonitorHLTsequence +
    hltPixelFakeDigiMonitoringSequence
)
