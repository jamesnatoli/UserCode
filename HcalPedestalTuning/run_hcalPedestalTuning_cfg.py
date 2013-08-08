import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")

process.load('Configuration.StandardSequences.Services_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.default.limit = 10
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometryDB_cff')

process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
   dump = cms.untracked.vstring(''),
   file = cms.untracked.string('')
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_P_V42_AN4::All"
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')
process.prefer("GlobalTag")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("HcalTBSource",
    fileNames = cms.untracked.vstring(
        'file:USC_213237.root'
    )
)

process.options = cms.untracked.PSet( 
    wantSummary = cms.untracked.bool(False) ) ## default is false

# HCAL digis
process.hcalDigis = cms.EDProducer("HcalRawToDigi",
    # Flag to enable unpacking of ZDC channels (default = false)
    UnpackZDC = cms.untracked.bool(False),
    # Optional filter to remove any digi with "data valid" off, "error" on, 
    # or capids not rotating
    FilterDataQuality = cms.bool(True),
    # Do not complain about missing FEDs
    ExceptionEmptyData = cms.untracked.bool(False),
    InputLabel = cms.InputTag("source"),
    # Use the defaults for FED numbers
    # Do not complain about missing FEDs
    ComplainEmptyData = cms.untracked.bool(False),
    # Flag to enable unpacking of calibration channels (default = false)
    UnpackCalib = cms.untracked.bool(True),
    # At most ten samples can be put into a digi, if there are more
    # than ten, firstSample and lastSample select which samples
    # will be copied to the digi
    firstSample = cms.int32(0),
    lastSample = cms.int32(9)
)

process.pedTuner = cms.EDAnalyzer('HcalPedestalTuning',
    digisName          = cms.string("hcalDigis"),
    mapIOV             = cms.int32(3),
    tunePedestals      = cms.bool(True),
    pathToOldXMLBricks = cms.untracked.string('Aug01_2013'),
    xmlBricksOutputDir = cms.untracked.string('test'),
    uniformDACSettings = cms.bool(False),
    initialDACSetting  = cms.int32(2),
    SiPMRBXes          = cms.vstring('HO1P10', 'HO2P12'),
    targetPedestalHPD  = cms.double(3.0),
    targetPedestalSiPM = cms.double(10.0),
    dac2adc            = cms.double(0.6),
    tagName            = cms.string('Aug08_2013')
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('alberto_histo.root')
)

#Paths
process.raw2digi_step = cms.Path(process.hcalDigis)
process.tuning_step = cms.Path(process.pedTuner)

# Schedule
process.schedule = cms.Schedule(process.raw2digi_step,process.tuning_step)

