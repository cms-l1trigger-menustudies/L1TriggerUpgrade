import FWCore.ParameterSet.Config as cms

from UserCode.L1TriggerUpgrade.l1Stage2NtupleFromRAW import *

from UserCode.L1TriggerUpgrade.MCSetup import *
mcSetup(process, False, True)

# job options
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

readFiles.extend( [
#' /store/mc/UpgFall13d/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/PU140bx25_POSTLS261_V3-v1/20000/008E2E98-0A39-E311-833F-0025905938D4.root'
#'/store/group/comm_trigger/L1TrackTrigger/620_SLHC10/Extended2023TTI/Electrons/PU140/m1_SingleElectron_E2023TTI_PU140.root'
    '/store/group/comm_trigger/L1TrackTrigger/620_SLHC12/Extended2023TTI/Neutrinos/PU140/NeutrinoGun_E2023TTI_PU140_1.root'
#    '/store/mc/Summer12/Neutrino_Pt2to20_gun/GEN-SIM-DIGI-RAW/UpgradeL1TDR-PU35_POSTLS161_V12-v2/00000/168A4B24-C53C-E211-A370-003048C68AA6.root'
    ] )

# Debug stuff

#process.output = cms.OutputModule(
#    "PoolOutputModule",
#    outputCommands = cms.untracked.vstring('keep *'),
#    fileName = cms.untracked.string('output.root')
#    )

#process.e = cms.EndPath(process.output)

#process.load("L1TriggerConfig.GctConfigProducers.l1GctConfigDump_cfi")
#process.MessageLogger.cout.placeholder = cms.untracked.bool(False)
#process.MessageLogger.cout.threshold = cms.untracked.string('DEBUG')
#process.MessageLogger.debugModules = cms.untracked.vstring('l1GctConfigDump')
#process.p += ( process.l1GctConfigDump )
