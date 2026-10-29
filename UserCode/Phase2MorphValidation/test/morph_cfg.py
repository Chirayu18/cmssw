import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.register("mode", "validate", VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "validate: GPU morphing on/off vs CPU SiPixelDigiMorphing; dump: stock clusterizer, dump digis")
options.register("cpuRows", 672, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                 "rows of the morphed modules, for the CPU SiPixelDigiMorphing")
options.register("cpuCols", 432, VarParsing.multiplicity.singleton, VarParsing.varType.int,
                 "columns of the morphed modules, for the CPU SiPixelDigiMorphing")
options.setDefault("maxEvents", 50)
options.setDefault("outputFile", "morph.root")
options.parseArguments()

from Configuration.Eras.Era_Phase2C26I13M9_cff import Phase2C26I13M9
process = cms.Process("MORPH", Phase2C26I13M9)

process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryExtendedRun4D128Reco_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Accelerators_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T35", "")

process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.SiPixelPhase2DigiToCluster = dict()
process.MessageLogger.cerr.INFO.limit = 0
process.MessageLogger.cerr.SiPixelPhase2DigiToCluster = cms.untracked.PSet(limit = cms.untracked.int32(10))

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(options.inputFiles))
process.maxEvents.input = options.maxEvents
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile))

# same settings as siPixelClustersPreSplitting_cff for Phase-2
from RecoLocalTracker.SiPixelClusterizer.siPixelPhase2DigiToCluster_cfi import siPixelPhase2DigiToCluster
from SimTracker.SiPhase2Digitizer.phase2TrackerDigitizer_cfi import PixelDigitizerAlgorithmCommon
stock = siPixelPhase2DigiToCluster.clone(
    Phase2ReadoutMode = PixelDigitizerAlgorithmCommon.Phase2ReadoutMode.value(),
    Phase2DigiBaseline = int(PixelDigitizerAlgorithmCommon.ThresholdInElectrons_Barrel.value()),
    Phase2KinkADC = 8,
    ElectronPerADCGain = PixelDigitizerAlgorithmCommon.ElectronPerAdc.value(),
)

if options.mode == "dump":
    process.gpuClusters = stock.clone()
    process.dump = cms.EDAnalyzer("Phase2ClusDump", gpuDigis = cms.InputTag("gpuClusters"))
    process.path = cms.Path(process.gpuClusters + process.dump)

else:
    # the charge cut would remove clusters after the clustering: switch it off, to compare the connectivity only
    region = ["1,1-99,1-99"]
    process.gpuMorphOff = stock.clone(
        clusterThreshold_layer1 = 0,
        clusterThreshold_otherLayers = 0,
    )
    process.gpuMorphOn = process.gpuMorphOff.clone(
        DoDigiMorphing = True,
        barrelRegions = region,
        endcapRegions = [],
        ValidationRecoveredDigis = True,
    )
    from RecoLocalTracker.SiPixelDigiReProducers.siPixelDigiMorphing_cfi import siPixelDigiMorphing
    process.cpuMorphed = siPixelDigiMorphing.clone(
        src = "simSiPixelDigis:Pixel",
        barrelRegions = region,
        endcapRegions = [],
        nrows = options.cpuRows,
        ncols = options.cpuCols,
        nrocs = 8,
    )
    process.validation = cms.EDAnalyzer("Phase2MorphValidation",
        digis = cms.InputTag("simSiPixelDigis", "Pixel"),
        cpuMorphed = cms.InputTag("cpuMorphed"),
        gpuOn = cms.InputTag("gpuMorphOn"),
        gpuOff = cms.InputTag("gpuMorphOff"),
        gpuRecovered = cms.InputTag("gpuMorphOn", "recovered"),
        maxDisplays = cms.int32(40),
    )
    process.path = cms.Path(process.gpuMorphOff + process.gpuMorphOn + process.cpuMorphed + process.validation)
