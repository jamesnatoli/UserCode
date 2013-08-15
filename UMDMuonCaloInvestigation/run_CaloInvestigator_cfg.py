import FWCore.ParameterSet.Config as cms

process = cms.Process('CALOINVESTIGATION')

#############################################################################
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')

process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')
#from Configuration.AlCa.autoCond import autoCond
#process.GlobalTag.globaltag = autoCond['startup']
#############################################################################

process.TFileService = cms.Service(
    "TFileService",
    fileName=cms.string('muon-calo_histograms.root')
    )

process.L1MuonCaloInv = cms.EDAnalyzer(
    'MuonCaloInvestigator',
    doGen = cms.untracked.bool(True),
    genSrc = cms.InputTag("genParticles"),
    rpcSrc = cms.InputTag("L1RPCbTFTrackConverter"),
    dttfSrc = cms.InputTag("L1DTTFTrackConverter"),
    hcalSrc = cms.InputTag("L1TMuonTriggerPrimitives","HCAL"),
    stdmuSrc = cms.InputTag("standAloneMuons"),
    glbmuSrc = cms.InputTag("globalMuons"),
    dRtruthToRpc  = cms.untracked.double(99.0), #0.2),
    dRrpcToDttf   = cms.untracked.double(99.0), #0.2),
    dRdttfToHcal  = cms.untracked.double(99.0), #0.2),
    dRhcalToStdMu = cms.untracked.double(99.0), #0.2),
    dRdttfToStdMu = cms.untracked.double(99.0)  #0.2)
)

infile = 'file:L1TMuon.root'

process.source = cms.Source(
    'PoolSource',
    fileNames = cms.untracked.vstring(infile)
    )

process.L1TMuonSequence = cms.Sequence(process.L1MuonCaloInv)

process.L1TMuonPath = cms.Path(process.L1TMuonSequence)

process.schedule = cms.Schedule(process.L1TMuonPath)
