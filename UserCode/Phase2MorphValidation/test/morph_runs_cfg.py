# Repeated runs of the Phase-2 clusterizer on the same events: reproducibility, overflow of the recovered pixels,
# and timing (run under nsys). master=True runs the stock module of a release without the digi morphing parameters.
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.register("master", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "release without morphing")
options.register("morphing", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "DoDigiMorphing")
options.register("maxFakes", 1000, VarParsing.multiplicity.singleton, VarParsing.varType.int, "MaxFakesInModule")
options.register("dump", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "dump digis and recovered pixels")
options.register("threads", 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "threads and streams")
options.register("backend", "", VarParsing.multiplicity.singleton, VarParsing.varType.string, "force a backend, e.g. serial_sync")
options.register("jsonFile", "", VarParsing.multiplicity.singleton, VarParsing.varType.string, "FastTimerService json")
options.setDefault("maxEvents", 100)
options.setDefault("outputFile", "runs.root")
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
process.MessageLogger.cerr.FwkReport.reportEvery = 50
process.options.numberOfThreads = options.threads
process.options.numberOfStreams = options.threads

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(options.inputFiles))
process.maxEvents.input = options.maxEvents
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile))

from RecoLocalTracker.SiPixelClusterizer.siPixelPhase2DigiToCluster_cfi import siPixelPhase2DigiToCluster
from SimTracker.SiPhase2Digitizer.phase2TrackerDigitizer_cfi import PixelDigitizerAlgorithmCommon
process.gpuClusters = siPixelPhase2DigiToCluster.clone(
    Phase2ReadoutMode = PixelDigitizerAlgorithmCommon.Phase2ReadoutMode.value(),
    Phase2DigiBaseline = int(PixelDigitizerAlgorithmCommon.ThresholdInElectrons_Barrel.value()),
    Phase2KinkADC = 8,
    ElectronPerADCGain = PixelDigitizerAlgorithmCommon.ElectronPerAdc.value(),
)
if not options.master:
    process.gpuClusters.DoDigiMorphing = options.morphing
    process.gpuClusters.barrelRegions = ["1,1-99,1-99"]
    process.gpuClusters.endcapRegions = []
    process.gpuClusters.MaxFakesInModule = options.maxFakes
    process.gpuClusters.ValidationRecoveredDigis = options.morphing and options.dump
if options.backend:
    process.gpuClusters.alpaka = cms.untracked.PSet(backend = cms.untracked.string(options.backend))

process.path = cms.Path(process.gpuClusters)
if options.dump:
    process.dump = cms.EDAnalyzer("Phase2ClusDump", gpuDigis = cms.InputTag("gpuClusters"))
    process.path += process.dump
    if options.morphing and not options.master:
        process.dumpRecovered = cms.EDAnalyzer("Phase2ClusDump", gpuDigis = cms.InputTag("gpuClusters", "recovered"))
        process.path += process.dumpRecovered

# per-module timing; run with CUDA_LAUNCH_BLOCKING=1 so that the time of the producer includes its kernels
options_json = options.jsonFile
if options_json:
    process.load("HLTrigger.Timer.FastTimerService_cfi")
    process.FastTimerService.enableDQM = False
    process.FastTimerService.printEventSummary = False
    process.FastTimerService.printRunSummary = False
    process.FastTimerService.printJobSummary = True
    process.FastTimerService.writeJSONSummary = True
    process.FastTimerService.jsonFileName = options_json
