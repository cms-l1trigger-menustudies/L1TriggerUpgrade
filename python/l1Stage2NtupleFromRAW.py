import FWCore.ParameterSet.Config as cms

# import master L1 Ntuple config
from UserCode.L1TriggerUpgrade.l1CurrentNtupleFromRAW import *

print "Adding Stage 2 algorithms"

# add full upgrade objects
process.load('UserCode.L1TriggerUpgrade.L1Upgrade_EGTau_cff')
process.load('UserCode.L1TriggerUpgrade.L1Upgrade_Jet_cff')
process.load('UserCode.L1TriggerUpgrade.l1UpgradeMuonIsolator_cfi')

# load the upgrade tree
process.load('UserCode.L1TriggerUpgrade.l1ExtraUpgradeTreeProducer_cfi')


#Need the L1Tracks and the L1Vertex
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.p += process.L1TkStubsFromPixelDigis 
# --- now we runn the L1Track producer :
process.BeamSpotFromSim =cms.EDProducer("BeamSpotFromSimProducer")

process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi') 
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TTrack_cfi")
process.L1Tracks.geometry = cms.untracked.string('BE5D')
process.p += process.BeamSpotFromSim*process.L1Tracks 

process.load("SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkPrimaryVertexProducer_cfi")
process.p += process.L1TrackPrimaryVertex      

# add upgrade algorithms to the path
process.p += process.SLHCCaloTrigger
process.p += process.l1UpgradeMuonIsolator
                                           

# bug fix for missing HCAL TPs in MC RAW
process.p.insert(1, process.valHcalTriggerPrimitiveDigis)
from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff import HcalTPGCoderULUT
HcalTPGCoderULUT.LUTGenerationMode = cms.bool(True)
process.valRctDigis.hcalDigis             = cms.VInputTag(cms.InputTag('valHcalTriggerPrimitiveDigis'))
process.L1CaloTowerProducer.HCALDigis =  cms.InputTag("valHcalTriggerPrimitiveDigis")

# --- Load the L1TkCaloSequence :

process.load("SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkCaloSequence_cff")

# -- Produce L1TkJets, HT and MHT from the L1Jets :
#process.L1TkCaloL1Jets = cms.Path( process.L1TkCaloSequence )
#process.p += process.L1TkCaloL1Jets
process.p += process.L1TkCaloSequence
# -- Produce the HLT JI Jets and L1TkJets, HT and MHT from  these jets :
#process.L1TkCaloHIJets = cms.Path( process.L1TkCaloSequenceHI )
#process.p += process.L1TkCaloHIJets
process.p += process.L1TkCaloSequenceHI
#
# ---------------------------------------------------------------------------
# TkEtMiss
process.load("SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkEtMissProducer_cfi")
process.p += process.L1TkEtMiss

# L1TkElectrons
process.load("SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkElectronTrackProducer_cfi")

# "Tkphotons" :
process.L1TkPhotons = cms.EDProducer("L1TkEmParticleProducer",
                                             label = cms.string("EGIsoTrk"), # labels the collection of L1TkEmParticleProducer that is produced.
                                                                                     # e.g. EG or IsoEG if all objects are kept, or
                                                                                     # EGIsoTrk or IsoEGIsoTrk if only the EG or IsoEG
                                                                                     # objects that pass a cut RelIso < RelIsoCut are written
                                                                                     # into the new collection.
                                             L1EGammaInputTag = cms.InputTag("SLHCL1ExtraParticlesNewClustering","EGamma"),
    #"SLHCL1ExtraParticles","EGamma"),      # input L1EG collection
                                                                                     # When the standard sequences are used :
                                                                                     #   - for the Run-1 algo, use ("l1extraParticles","NonIsolated")
                                                                                     #     or ("l1extraParticles","Isolated")
                                                                                     #   - for the "old stage-2" algo (2x2 clustering), use
                                                                                     #     ("SLHCL1ExtraParticles","EGamma") or ("SLHCL1ExtraParticles","IsoEGamma")
                                                                                     #   - for the new clustering algorithm of Jean-Baptiste et al,
                                                                                     #     use ("SLHCL1ExtraParticlesNewClustering","IsoEGamma") or
                                                                                     #     ("SLHCL1ExtraParticlesNewClustering","EGamma").
                                             ETmin = cms.double( -1 ),               # Only the L1EG objects that have ET > ETmin in GeV
                                                                                     # are considered. ETmin < 0 means that no cut is applied.
                                             RelativeIsolation = cms.bool( True ),   # default = True. The isolation variable is relative if True,
                                                                                     # else absolute.
                                             IsoCut = cms.double( 0.1 ),             # Cut on the (Trk-based) isolation: only the L1TkEmParticle for which
                                                                                     # the isolation is below RelIsoCut are written into
                                                                                     # the output collection. When RelIsoCut < 0, no cut is applied.
                                                                                     # When RelativeIsolation = False, IsoCut is in GeV.
                                                # Determination of the isolation w.r.t. L1Tracks :
                                             L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
                                             ZMAX = cms.double( 25. ),       # in cm
                                             CHI2MAX = cms.double( 100. ),
                                             PTMINTRA = cms.double( 2. ),    # in GeV
                                             DRmin = cms.double( 0.07),
                                             DRmax = cms.double( 0.30 ),
                                             PrimaryVtxConstrain = cms.bool( False ),  # default = False
                                             DeltaZMax = cms.double( 999. ),    # in cm. Used only when PrimaryVtxConstrain = True
                                             L1VertexInputTag = cms.InputTag("NotUsed"),     # Used only when PrimaryVtxConstrain = True
                                     )

process.p += process.L1TkPhotons

process.p += process.L1TkElectrons

# Isolated (w.r.t. L1Tracks) electrons :
process.L1TkIsoElectrons = process.L1TkElectrons.clone()
process.L1TkIsoElectrons.IsoCut = cms.double(0.1)
process.p += process.L1TkIsoElectrons

# on top of new stage-2 :
process.L1TkElectronsNewclus = process.L1TkElectrons.clone()
process.L1TkElectronsNewclus.L1EGammaInputTag = cms.InputTag("SLHCL1ExtraParticlesNewClustering","EGamma")
process.p += process.L1TkElectronsNewclus

process.L1TkIsoElectronsNewclus = process.L1TkIsoElectrons.clone()
process.L1TkIsoElectronsNewclus.L1EGammaInputTag = cms.InputTag("SLHCL1ExtraParticlesNewClustering","EGamma")
process.p += process.L1TkIsoElectronsNewclus

#L1TkMuons
 #--- Now run the L1TkMuonParticleProducer
process.load("SLHCUpgradeSimulations.L1TrackTriggerObjects.L1TkMuonProducer_cfi")
process.p += process.L1TkMuons

# --- Now run the L1TkTauParticle producer(s)

process.L1TkTauFromCalo = cms.EDProducer("L1TkTauFromCaloProducer",
                                                 L1TausInputTag = cms.InputTag("SLHCL1ExtraParticles","Taus"),
                                                 L1TrackInputTag = cms.InputTag("L1Tracks","Level1TkTracks"),
                                                 ZMAX = cms.double( 25. ),       # in cm
                                                 CHI2MAX = cms.double( 100. ),
                                                 PTMINTRA = cms.double( 2. ),    # in GeV
                                                 DRmax = cms.double( 0.5 ),
                                                 nStubsmin = cms.int32( 5 )        # minimum number of stubs
                                         )
process.p += process.L1TkTauFromCalo 


process.p += process.l1ExtraUpgradeTreeProducer    # upgrade candidates

print "Done"
