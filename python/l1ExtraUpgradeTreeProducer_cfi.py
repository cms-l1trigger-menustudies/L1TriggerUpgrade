import FWCore.ParameterSet.Config as cms

l1ExtraUpgradeTreeProducer = cms.EDAnalyzer("L1ExtraUpgradeTreeProducer",
   egLabel = cms.untracked.InputTag("SLHCL1ExtraParticlesNewClustering","EGamma"),
   isoEGLabel = cms.untracked.InputTag("SLHCL1ExtraParticlesNewClustering","IsoEGamma"),

   tkEGLabel = cms.untracked.InputTag("L1TkElectronsNewclus","EG"),
   tkIsoEGLabel = cms.untracked.InputTag("L1TkIsoElectronsNewclus","EG"),
   tkEMLabel = cms.untracked.InputTag("L1TkPhotons","EGIsoTrk"),
    
   tauLabel = cms.untracked.InputTag("SLHCL1ExtraParticles:Taus"),
   isoTauLabel = cms.untracked.InputTag("SLHCL1ExtraParticles:IsoTaus"),
   jetLabel = cms.untracked.InputTag("L1TowerJetPUSubtractedProducer:PUSubCen8x8"),

   tkMuonLabel = cms.untracked.InputTag("L1TkMuons",""),                                            
   tkJetLabel = cms.untracked.InputTag("L1TkJets","Central"),                                            
   
   fwdJetLabel = cms.untracked.InputTag("l1extraParticles","Forward"),
   #L1CalibFilterTowerJetProducer:Fwd8x8"),
   muonLabel = cms.untracked.InputTag("l1UpgradeMuonIsolator"),
   metLabel = cms.untracked.InputTag("rawSLHCL1ExtraParticles:MET"),
   tkMetLabel = cms.untracked.InputTag("L1TkEtMiss","MET"),
                                            
   mhtLabel = cms.untracked.InputTag("L1CalibFilterTowerJetProducer:TowerMHT"),
   maxL1Extra = cms.uint32(20)
)
