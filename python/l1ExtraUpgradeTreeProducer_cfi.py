import FWCore.ParameterSet.Config as cms

l1ExtraUpgradeTreeProducer = cms.EDAnalyzer("L1ExtraUpgradeTreeProducer",
   egLabel = cms.untracked.InputTag("SLHCL1ExtraParticlesNewClustering","EGamma"),
   isoEGLabel = cms.untracked.InputTag("SLHCL1ExtraParticlesNewClustering","IsoEGamma"),

   tkEGLabel = cms.untracked.InputTag("L1TkElectronsNewclus","EG"),
   tkEG2Label = cms.untracked.InputTag("L1TkElectronsLoose","EG"),

   tkIsoEGLabel = cms.untracked.InputTag("L1TkIsoElectronsNewclus","EG"),
   tkEMLabel = cms.untracked.InputTag("L1TkPhotonsNewclus","EG"),
    
   tauLabel = cms.untracked.InputTag("SLHCL1ExtraParticles:Taus"),
   isoTauLabel = cms.untracked.InputTag("SLHCL1ExtraParticles:IsoTaus"),
   tkTauLabel = cms.untracked.InputTag("L1TkTauFromL1Track",""),                                            

#   jetLabel = cms.untracked.InputTag("L1TowerJetPUSubtractedProducer:PUSubCen8x8"),
   jetLabel = cms.untracked.InputTag("L1CalibFilterTowerJetProducer","CalibratedTowerJets"),   # Calib L1jets                                            

   tkMuonLabel = cms.untracked.InputTag("L1TkMuonsMerge",""),                                            

   tkJetLabel = cms.untracked.InputTag("L1TkJets","Central"),                                            
#   tkJetLabel = cms.untracked.InputTag("L1TkJetsHI","Central"),  ## Use the HI HLT version of L1Jets                                          

   fwdJetLabel = cms.untracked.InputTag("l1extraParticles","Forward"),

   muonLabel = cms.untracked.InputTag("l1extraParticles"), #l1UpgradeMuonIsolator"),
   metLabel = cms.untracked.InputTag("L1EnergySumProducer","MET"),
   #metLabel = cms.untracked.InputTag("rawSLHCL1ExtraParticles:MET"),
   tkMetLabel = cms.untracked.InputTag("L1TkEtMiss","MET"),

   mhtLabel = cms.untracked.InputTag("L1CalibFilterTowerJetProducer:TowerMHT"),
   tkMhtLabel = cms.untracked.InputTag("L1TkHTMissVtx",""),
#   tkMhtLabel = cms.untracked.InputTag("L1TkHTMissCalo",""),                                            
# tkMhtLabel = cms.untracked.InputTag("L1TkHTMissCaloHI",""), # HI version                                            
   maxL1Extra = cms.uint32(20)
)
