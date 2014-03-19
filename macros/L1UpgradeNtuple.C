#define l1ntuple_cxx
#include "L1UpgradeNtuple.h"

double L1UpgradeNtuple::deltaPhi(double phi1, double phi2) {
  double result = phi1 - phi2;
  if(fabs(result) > 9999) return result;
  while (result > TMath::Pi()) result -= 2*TMath::Pi();
  while (result <= -TMath::Pi()) result += 2*TMath::Pi();
  return result;
}

double L1UpgradeNtuple::deltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return sqrt(deta*deta + dphi*dphi);
}

void L1UpgradeNtuple::Test()
{ 

  if (fChain == 0)  return;
 
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  unsigned int nevents =0;

  std::cout << nentries << " events to process"<<std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    //fChain->GetEvent(jentry);
  
    nevents++;
    if (nevents<9)  //eight first events
      { 
    std::cout << "--------------------- Event "<<jentry<<" ---------------------"<<std::endl;

    //event_
    std::cout << "L1Tree         : event_->run = "<<event_->run<<std::endl;

    //gct_
    std::cout << "L1Tree         : gct_->IsoEmEta.Size() = "<<gct_->IsoEmEta.size()<<std::endl;
    if (gct_->IsoEmEta.size()!=0)
      std::cout << "L1Tree         : gct->IsoEmEta[0] = "<<gct_->IsoEmEta[0]<<std::endl; 

    //gmt_
    std::cout << "L1Tree         : gmt_->Ndt        = "<<gmt_->Ndt<<std::endl;
    if (gmt_->Ndt!=0)
      std::cout << "L1Tree         : gmt->Bxdt[0] = "<<gmt_->Bxdt[0]<<std::endl;

    //gt_
    std::cout << "L1Tree         : gt_->tw1.size = "<<gt_->tw1.size()<<std::endl;
    if (gt_->tw1.size()!=0)
      std::cout << "L1Tree         : gt->tw1[0] = "<<gt_->tw1[0]<<std::endl;

    //rct_
    std::cout << "L1Tree         : rct->RegSize    = "<<rct_->RegSize<<std::endl;
    if (rct_->RegSize!=0)
      std::cout << "L1Tree         : rct->RegEta[0] = "<<rct_->RegEta[0]<<std::endl; 

    //dttf_
    std::cout << "L1Tree         : dttf->trSize     = "<<dttf_->trSize<<std::endl;
    if (dttf_->trSize!=0)
      std::cout << "L1Tree         : dttf->trBx[0] = "<<dttf_->trBx[0]<<std::endl;

    //csctf_
    std::cout << "L1Tree         : csctf_->lctSize    = "<<csctf_->lctSize<<std::endl;
    if (csctf_->lctBx.size()!=0)
      std::cout << "L1Tree         : csctf->lctBx[0] = "<<csctf_->lctBx[0]<<std::endl;

    //recoMuon
    if (domuonreco) std::cout << "L1MuonRecoTree : nb muons   = " << recoMuon_->nMuons << std::endl;

    //recoMet
    if (doreco)     std::cout << "L1RecoTree     : met        = " << recoMet_->met     << std::endl;
  
    //recoJet_
    if (doreco)     std::cout << "L1RecoTree     : nb jets    = " << recoJet_->nJets   << std::endl;

    //recoBasicCluster_
    if (doreco)     std::cout << "L1RecoTree     : nb BasicCluster   = " << recoBasicCluster_->nClusters   << std::endl;

    //recoSuperCluster_
    if (doreco)     std::cout << "L1RecoTree     : nb SuperCluster   = " << recoSuperCluster_->nClusters   << std::endl;

    //recoTrack_
    if (doreco)     std::cout << "L1RecoTree     : nTrk       = " << recoTrack_->nTrk     << std::endl;

    //recoVertex
    if (doreco)     std::cout << "L1RecoTree     : nVtx       = " << recoVertex_->nVtx     << std::endl;    

    //l1extra_
    if (dol1extra)  std::cout << "L1ExtraTree    : et         = " << l1extra_->et[0]      << std::endl;

    // l1menu
    if (dol1menu)   std::cout << "L1MenuTree     : Is PrescaleIndex for AlgoTrigg valid ? " << l1menu_->AlgoTrig_PrescaleFactorIndexValid << std::endl;

    if (dol1menu)   std::cout << "L1MenuTree     : PrescaleIndex for AlgoTrigg  " << l1menu_->AlgoTrig_PrescaleFactorIndex << std::endl;

    if (dol1menu)   std::cout << "L1MenuTree     : Is PrescaleIndex for TechTrigg valid ? " << l1menu_->TechTrig_PrescaleFactorIndexValid << std::endl;

    if (dol1menu)   std::cout << "L1MenuTree     : PrescaleIndex for TechTrigg  " << l1menu_->TechTrig_PrescaleFactorIndex << std::endl;

    //    std::cout << "L1ExtraUpgrade  : nEG = " << l1upgrade_->nEG << std::endl;

    }
  }
   
}




void L1UpgradeNtuple::Test2()
{ 

  if (fChain == 0)  return;
 
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  //unsigned int nevents =0;

  std::cout << nentries << " events to process"<<std::endl;

//   TCanvas * can = new TCanvas;
//   can->Divide(2,1);

  TFile *total=new TFile("total.root","RECREATE");

  TH1D *hGen = new TH1D("gelectrons","",25,0.,50.);
  TH1D *hUpgrade1 = new TH1D("hupgrade1","",25,0.,50.); TH1D *hUpgrade2 = new TH1D("hupgrade2","",25,0.,50.);
  
  TH1D *hGENeta = new TH1D("hgeneta","",30,-3.,3.);
  TH1D *hEGeta = new TH1D("hegeta","",30,-3.,3.); TH1D *hisoEGeta = new TH1D("hisoegeta","",30,-3.,3.);

  TH1F *DRmin1 = new TH1F("drmin1","",100,0.,5.); TH1F *DRmin2 = new TH1F("drmin2","",100,0.,5.);

  TH2F *dPTvsdR1 = new TH2F("dptvsdr1","",100,0.,1.,100,0.,2.);
  TH2F *dPTvsdR2 = new TH2F("dptvsdr2","",100,0.,1.,100,0.,2.);

  TH1D *hJETgen = new TH1D("hjetgen","",100,0.,1000.);
  TH1D *hJet25 = new TH1D("hjet25","",100,0.,1000.);
  TH1D *hJet50 = new TH1D("hjet50","",100,0.,1000.);
  TH1D *hJet75 = new TH1D("hjet75","",100,0.,1000.);

  TH1D *hTaugen = new TH1D("htaugen","",100,0.,1000.);
  TH1D *hTau30 = new TH1D("htau30","",100,0.,1000.);
  TH1D *hTau50 = new TH1D("htau50","",100,0.,1000.);
  TH1D *hTau100 = new TH1D("htau100","",100,0.,1000.);
  TH1D *hTau150 = new TH1D("htau150","",100,0.,1000.);
  TH1D *hTau125 = new TH1D("htau125","",100,0.,1000.);

  TH1D *hisoTaugen = new TH1D("hisotaugen","",100,0.,1000.);
  TH1D *hIsoTau30 = new TH1D("hisotau30","",100,0.,1000.);
  TH1D *hIsoTau50 = new TH1D("hisotau50","",100,0.,1000.);
  TH1D *hIsoTau100 = new TH1D("hisotau100","",100,0.,1000.);
  TH1D *hIsoTau150 = new TH1D("hisotau150","",100,0.,1000.);
  TH1D *hIsoTau125 = new TH1D("hisotau125","",100,0.,1000.);

  TH1D *hJet100 = new TH1D("hjet100","",100,0.,1000.);
  TH1D *hJet150 = new TH1D("hjet150","",100,0.,1000.);
  TH1D *hJet200 = new TH1D("hjet200","",100,0.,1000.);
  TH1D *hJet250 = new TH1D("hjet250","",100,0.,1000.);
  TH1D *hJet300 = new TH1D("hjet300","",100,0.,1000.);

  TH1D *hEGgen = new TH1D("heggen","",25,0.,50.);
  TH1D *hEG12 = new TH1D("heg12","",25,0.,50.);
  TH1D *hEG15 = new TH1D("heg15","",25,0.,50.);
  TH1D *hEG20 = new TH1D("heg20","",25,0.,50.);
  TH1D *hEG25 = new TH1D("heg25","",25,0.,50.);
  TH1D *hEG30 = new TH1D("heg30","",25,0.,50.);

  TH1D *hisoEGgen = new TH1D("hisoeggen","",25,0.,50.);
  TH1D *hisoEG12 = new TH1D("hisoeg12","",25,0.,50.);
  TH1D *hisoEG15 = new TH1D("hisoeg15","",25,0.,50.);
  TH1D *hisoEG20 = new TH1D("hisoeg20","",25,0.,50.);
  TH1D *hisoEG25 = new TH1D("hisoeg25","",25,0.,50.);
  TH1D *hisoEG30 = new TH1D("hisoeg30","",25,0.,50.);

  std::map<int,double> taus;
  std::map<int,double> partons;
  std::vector<double> electrons;

  TLorentzVector electron(0.,0.,0.,0.);
  TLorentzVector tau(0.,0.,0.,0.);
  TLorentzVector parton;
  //TLorentzVector egamma, isoegamma;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    nb = fChain->GetEntry(jentry); nbytes += nb;
    
  //   std::cout << "" << std::endl;
//     std::cout << "Event = " << int(jentry) << std::endl;
    
    electrons.clear();
    partons.clear(); taus.clear();

    // genParticles (==Reference collection to efficiency studies)
    // Only pick-up the highest-pT particle - and then loop over the L1 collection to find a match in DR.
    int itau=-1;
    int idp;
    for (int i=0; i<gen_->id.size(); i++) {
      // Exclude initial state partons at 0,1 positions
      //      if (gen_->status.at(i)==3) { //continue; }
      //    std::cout << "PParticle id = " << gen_->id.at(i) << " has status = " << gen_->status.at(i) << std::endl;
      //} else { continue; }
      //if (gen_->status.at(i)==3) { 

      // Fill electrons
      if (abs(gen_->id.at(i))==11) {
	electron.SetPxPyPzE(gen_->px.at(i),gen_->py.at(i),gen_->pz.at(i),gen_->e.at(i));
      }
      // Fill hadronic-taus
      if (gen_->status.at(i)==3) { 
	if  (abs(gen_->id.at(i))==15  && abs(gen_->parent_id.at(i))==24) {
	  tau.SetPxPyPzE(gen_->px.at(i),gen_->py.at(i),gen_->pz.at(i),gen_->e.at(i));
	  taus[i]=tau.Pt(); itau=1;
	  //  std::cout << "PParticle id = " << gen_->id.at(i) << " has status = " << gen_->status.at(i) << " and parent ID: " << gen_->parent_id.at(i) << std::endl;
	}
      }
      if (itau>0) {
// 	if (gen_->status.at(i)==2) { 
// 	  if  (abs(gen_->id.at(i))==15) {
// 	    //	    tau.SetPxPyPzE(gen_->px.at(i),gen_->py.at(i),gen_->pz.at(i),gen_->e.at(i));
// 	    //	    taus[i]=tau.Pt(); itau=1;
// 	    std::cout << "PParticle id = " << gen_->id.at(i) << " has status = " << gen_->status.at(i) << " and parent ID: " << gen_->parent_id.at(i) << std::endl;
// 	  }
// 	}
	if (gen_->status.at(i)==3) { 
	  if  ( (abs(gen_->id.at(i))==11 || abs(gen_->id.at(i))==13) && abs(gen_->parent_id.at(i))==15) {
	    //	  tau.SetPxPyPzE(gen_->px.at(i),gen_->py.at(i),gen_->pz.at(i),gen_->e.at(i));
	    //	  taus[i]=tau.Pt(); itau=1;
	    std::cout << "PParticle id = " << gen_->id.at(i) << " has status = " << gen_->status.at(i) << " and parent ID: " << gen_->parent_id.at(i) << std::endl;
	  }
	}
      }


      // Fill partons
      idp=abs(gen_->id.at(i));
      if (gen_->status.at(i)==3) { //continue; }
	if (idp==1 || idp==2 || idp==3 || idp==4 || idp==5 || idp==21) {
	  parton.SetPxPyPzE(gen_->px.at(i),gen_->py.at(i),gen_->pz.at(i),gen_->e.at(i));
	  if (fabs(parton.Eta())<=3.) { partons[i]=parton.Pt(); }
	}
      }
    }

    // Sort taus in pT
    std::map<int, double>::iterator it = taus.begin();
    double maxpt = it->second;
    int maxipt=-1;

    for (it=taus.begin(); it!=taus.end(); ++it) {
       if (it->second>=maxpt) {
	maxpt=it->second; maxipt=it->first;
      }
       //  std::cout << "Particle at " << it->first << " has pT = " << it->second << std::endl;
    }
    if (maxipt>=0) { 
      tau.SetPxPyPzE(gen_->px.at(maxipt),gen_->py.at(maxipt),gen_->pz.at(maxipt),gen_->e.at(maxipt));
      // std::cout << "Highest-pT tau at "  << maxipt << " has pT = " << maxpt << std::endl;

      int ijt;
      int iTau=-1;
      double dr; 
      double drmax=0.5;
      
       // Loop over Taus
      //  std::cout << "We have : " << l1upgrade_->nIsoTau << " taus" << std::endl;
      for (Int_t i=0; i<l1upgrade_->nTau; i++) {

// 	std::cout << "Counting taus = " << i << std::endl;
// 	std::cout << "IsoTau et = " << l1upgrade_->isoTauEt.at(i) << std::endl;
 	ijt=i;
// 	if (l1upgrade_->nIsoTau>1) {
// 	  if (l1upgrade_->isoTauEt.at(i+1) > l1upgrade_->isoTauEt.at(i)) { ijt=i+1; }
// 	} 
	if (i>0) break;
	dr=deltaR(tau.Eta(),tau.Phi(),l1upgrade_->tauEta.at(ijt),l1upgrade_->tauPhi.at(ijt));

	if (dr<drmax) { drmax=dr; iTau=ijt; }
      } // End loop L1jets

      if (iTau>=0) { // gen-L1 jet match
	hTaugen->Fill(tau.Pt());

	if (l1upgrade_->tauEt.at(iTau)>=30.) hTau30->Fill(tau.Pt());
	if (l1upgrade_->tauEt.at(iTau)>=50.) hTau50->Fill(tau.Pt());
	if (l1upgrade_->tauEt.at(iTau)>=100.) hTau100->Fill(tau.Pt());
	if (l1upgrade_->tauEt.at(iTau)>=150.) hTau150->Fill(tau.Pt());
	if (l1upgrade_->tauEt.at(iTau)>=125.) hTau125->Fill(tau.Pt());

      }
    

      iTau=-1;
      drmax=0.5;
      // Loop over IsoTaus
      //  std::cout << "We have : " << l1upgrade_->nIsoTau << " taus" << std::endl;
      for (Int_t i=0; i<l1upgrade_->nIsoTau; i++) {

// 	std::cout << "Counting taus = " << i << std::endl;
// 	std::cout << "IsoTau et = " << l1upgrade_->isoTauEt.at(i) << std::endl;
 	ijt=i;
// 	if (l1upgrade_->nIsoTau>1) {
// 	  if (l1upgrade_->isoTauEt.at(i+1) > l1upgrade_->isoTauEt.at(i)) { ijt=i+1; }
// 	} 
	if (i>0) break;
	dr=deltaR(tau.Eta(),tau.Phi(),l1upgrade_->isoTauEta.at(ijt),l1upgrade_->isoTauPhi.at(ijt));

	if (dr<drmax) { drmax=dr; iTau=ijt; }
      } // End loop L1jets

      if (iTau>=0) { // gen-L1 jet match
	hisoTaugen->Fill(tau.Pt());

	if (l1upgrade_->isoTauEt.at(iTau)>=30.) hIsoTau30->Fill(tau.Pt());
	if (l1upgrade_->isoTauEt.at(iTau)>=50.) hIsoTau50->Fill(tau.Pt());
	if (l1upgrade_->isoTauEt.at(iTau)>=100.) hIsoTau100->Fill(tau.Pt());
	if (l1upgrade_->isoTauEt.at(iTau)>=150.) hIsoTau150->Fill(tau.Pt());
	if (l1upgrade_->isoTauEt.at(iTau)>=125.) hIsoTau125->Fill(tau.Pt());

      }
      
      
    }
    
    // Sort partons in PT
    std::map<int, double>::iterator tt = partons.begin();
    maxpt = tt->second;
    maxipt=-1;
    for (tt=partons.begin(); tt!=partons.end(); ++tt) {
      if (tt->second>=maxpt) {
	maxpt=tt->second; maxipt=tt->first;
      }
    }
    if (maxipt>0) {      // Found highest-pT genJet
      //      std::cout << "Highest-pT partonJet at "  << maxipt << " has pT = " << maxpt << std::endl;
      parton.SetPxPyPzE(gen_->px.at(maxipt),gen_->py.at(maxipt),gen_->pz.at(maxipt),gen_->e.at(maxipt));

      int ij;
      int iJet=-1;
      double dr; 
      double drmax=0.5;
      // Loop over L1jets
      for (Int_t i=0; i<l1upgrade_->nJets; i++) {

	//	std::cout << "L1Jet et = " << l1upgrade_->jetEt.at(i) << std::endl;
	ij=i;
	if (l1upgrade_->nJets>1) {
	  if (l1upgrade_->jetEt.at(i+1) > l1upgrade_->jetEt.at(i)) { ij=i+1; }
	} 
	if (i>0) break;
	dr=deltaR(parton.Eta(),parton.Phi(),l1upgrade_->jetEta.at(ij),l1upgrade_->jetPhi.at(ij));

	if (dr<drmax) { drmax=dr; iJet=ij; }
      } // End loop L1jets

      if (iJet>=0) { // gen-L1 jet match
	hJETgen->Fill(parton.Pt());

	if (l1upgrade_->jetEt.at(iJet)>=25.) hJet25->Fill(parton.Pt());
	if (l1upgrade_->jetEt.at(iJet)>=50.) hJet50->Fill(parton.Pt());
	if (l1upgrade_->jetEt.at(iJet)>=75.) hJet75->Fill(parton.Pt());
	
	if (l1upgrade_->jetEt.at(iJet)>=100.) hJet100->Fill(parton.Pt());
	if (l1upgrade_->jetEt.at(iJet)>=150.) hJet150->Fill(parton.Pt());
	if (l1upgrade_->jetEt.at(iJet)>=200.) hJet200->Fill(parton.Pt());
	if (l1upgrade_->jetEt.at(iJet)>=250.) hJet250->Fill(parton.Pt());
	if (l1upgrade_->jetEt.at(iJet)>=300.) hJet300->Fill(parton.Pt());
      }

    } // end if found highest-pT ref (parton-)Jet

    //    if (electron.Pt()==0.) continue;
    hGen->Fill(electron.Pt()); hGENeta->Fill(electron.Eta());


    double dptall,drall;

    // loop over Egamma
    int iEG=-1;

    double dR;
    double dRmax=0.5;
    for (Int_t i=0; i<l1upgrade_->nTkEG; i++) {
      // std::cout << "EG(" << i << ") et= " << l1upgrade_->egEt.at(i) << std::endl;
      if (i>0) break;

      //  if (l1upgrade_->egEt.at(i)<5.) continue;
      //  std::cout << "EG(" << i << ") et= " << l1upgrade_->egEt.at(i) << std::endl;

      dptall=fabs(electron.Pt()-l1upgrade_->tkEGEt.at(i))/electron.Pt();
      dR=deltaR(electron.Eta(),electron.Phi(),l1upgrade_->tkEGEta.at(i),l1upgrade_->tkEGPhi.at(i));

      dPTvsdR1->Fill(dR,dptall);

      if (dR<dRmax) { dRmax=dR; iEG=i; }

    }
    if (iEG>=0) {  // if found a match within 0.5
      DRmin1->Fill(dRmax);

     //  hGen->Fill(electron.Pt()); hGENeta->Fill(electron.Eta());
      //      std::cout << "Found a match within DR= " << dRmax << std::endl; 
      hEGgen->Fill(electron.Pt()); 
      
      if (l1upgrade_->tkEGEt.at(iEG)>12.) hEG12->Fill(electron.Pt());
      if (l1upgrade_->tkEGEt.at(iEG)>15.) hEG15->Fill(electron.Pt());
      if (l1upgrade_->tkEGEt.at(iEG)>20.) hEG20->Fill(electron.Pt());
      if (l1upgrade_->tkEGEt.at(iEG)>25.) hEG25->Fill(electron.Pt());
      if (l1upgrade_->tkEGEt.at(iEG)>30.) hEG30->Fill(electron.Pt());

      hUpgrade1->Fill(l1upgrade_->tkEGEt.at(iEG));
      hEGeta->Fill(l1upgrade_->tkEGEta.at(iEG));
    }

    
  // Start with a gen electron & loop over isoEgamma
    iEG=-1;

    //double dR;
    dRmax=0.5;
    for (Int_t i=0; i<l1upgrade_->nTkIsoEG; i++) {

      if (i>0) break;
      // if (l1upgrade_->tkIsoEGEt.at(i)<5.) continue;

      dptall=fabs(electron.Pt()-l1upgrade_->tkIsoEGEt.at(i))/electron.Pt();
      dR=deltaR(electron.Eta(),electron.Phi(),l1upgrade_->tkIsoEGEta.at(i),l1upgrade_->tkIsoEGPhi.at(i));

      dPTvsdR2->Fill(dR,dptall);

      if (dR<dRmax) { dRmax=dR; iEG=i; }

    }
    if (iEG>=0) { DRmin2->Fill(dRmax);
      //      std::cout << "Found a match within DR= " << dRmax << std::endl; 
      hisoEGgen->Fill(electron.Pt());

      if (l1upgrade_->tkIsoEGEt.at(iEG)>12.) hisoEG12->Fill(electron.Pt());
      if (l1upgrade_->tkIsoEGEt.at(iEG)>15.) hisoEG15->Fill(electron.Pt());
      if (l1upgrade_->tkIsoEGEt.at(iEG)>20.) hisoEG20->Fill(electron.Pt());
      if (l1upgrade_->tkIsoEGEt.at(iEG)>25.) hisoEG25->Fill(electron.Pt());
      if (l1upgrade_->tkIsoEGEt.at(iEG)>30.) hisoEG30->Fill(electron.Pt());

      hUpgrade2->Fill(l1upgrade_->tkIsoEGEt.at(iEG));
      hisoEGeta->Fill(l1upgrade_->tkIsoEGEta.at(iEG));
    }

  }
  
  hGen->Sumw2(); hGENeta->Sumw2();
  hUpgrade1->Sumw2(); hUpgrade2->Sumw2();
  hEGeta->Sumw2(); hisoEGeta->Sumw2();

  hEGgen->Sumw2(); hisoEGgen->Sumw2();
  hEG12->Sumw2(); hisoEG12->Sumw2();
  hEG15->Sumw2(); hisoEG15->Sumw2();
  hEG20->Sumw2(); hisoEG20->Sumw2();
  hEG25->Sumw2(); hisoEG25->Sumw2();
  hEG30->Sumw2(); hisoEG30->Sumw2();

  hTaugen->Sumw2();
  hTau30->Sumw2(); hTau50->Sumw2(); hTau100->Sumw2();
  hTau150->Sumw2(); hTau125->Sumw2();
  hTaugen->SetDirectory(total);
  hTau30->SetDirectory(total); hTau50->SetDirectory(total); hTau100->SetDirectory(total);
  hTau150->SetDirectory(total); hTau125->SetDirectory(total);

  hisoTaugen->Sumw2();
  hIsoTau30->Sumw2(); hIsoTau50->Sumw2(); hIsoTau100->Sumw2();
  hIsoTau150->Sumw2(); hIsoTau125->Sumw2();
  hTaugen->SetDirectory(total);
  hIsoTau30->SetDirectory(total); hIsoTau50->SetDirectory(total); hIsoTau100->SetDirectory(total);
  hIsoTau150->SetDirectory(total); hIsoTau125->SetDirectory(total);

  hJETgen->Sumw2();
  hJet25->Sumw2(); hJet50->Sumw2(); hJet75->Sumw2();
  hJet100->Sumw2(); hJet150->Sumw2(); hJet200->Sumw2();
  hJet250->Sumw2(); hJet300->Sumw2();
  hJETgen->SetDirectory(total);
  hJet25->SetDirectory(total); hJet50->SetDirectory(total); hJet75->SetDirectory(total);
  hJet100->SetDirectory(total); hJet150->SetDirectory(total); hJet200->SetDirectory(total);
  hJet250->SetDirectory(total); hJet300->SetDirectory(total);

  hEGgen->SetDirectory(total); hisoEGgen->SetDirectory(total);
  hEG12->SetDirectory(total); hisoEG12->SetDirectory(total);
  hEG15->SetDirectory(total); hisoEG15->SetDirectory(total);
  hEG20->SetDirectory(total); hisoEG20->SetDirectory(total);
  hEG25->SetDirectory(total); hisoEG25->SetDirectory(total);
  hEG30->SetDirectory(total); hisoEG30->SetDirectory(total);

  total->Write(); total->Close();

 //  TGraphAsymmErrors *ratio1 = new TGraphAsymmErrors(hUpgrade1,hGen,"n");
//   TGraphAsymmErrors *ratio2 = new TGraphAsymmErrors(hUpgrade2,hGen,"n");
//  TCanvas * c1 = new TCanvas();
//  c1->Divide(2,1);

  //  gStyle->SetOptFit(0);

//   TGraphAsymmErrors *effisoEG12 = new TGraphAsymmErrors(hisoEG12,hisoEGgen);
//   TGraphAsymmErrors *effisoEG15 = new TGraphAsymmErrors(hisoEG15,hisoEGgen);
//   TGraphAsymmErrors *effisoEG20 = new TGraphAsymmErrors(hisoEG20,hisoEGgen);

  
  // TGraphAsymmErrors *effEG12=drawEff(hEG12,hEGgen,"E_{T}>12","EG E_{T} (GeV)",0,10.,50.,0.,50.,kBlack,20,"",12.);
  // TGraphAsymmErrors *effEG15=drawEff(hEG15,hEGgen,"E_{T}>15","EG E_{T} (GeV)",0,10.,50.,0.,50.,kBlue-6,20,"",15.);
  // TGraphAsymmErrors *effEG20=drawEff(hEG20,hEGgen,"E_{T}>20","EG E_{T} (GeV)",0,10.,50.,0.,50.,kGreen-1,20,"",20.);

  return;

  // EG efficiency
  TH1F *ratio1=new TH1F("ratio1","",25,0.,50.);
  ratio1=(TH1F*)hUpgrade1->Clone();
  ratio1->Divide(hGen); 
  // isoEG efficiency
  TH1F *ratio2=new TH1F("ratio2","",25,0.,50.);
  ratio2=(TH1F*)hUpgrade2->Clone();
  ratio2->Divide(hGen);  

  // EGeta eff
  TH1F *eratio1=new TH1F("eratio1","",25,0.,50.);
  eratio1=(TH1F*)hEGeta->Clone();
  eratio1->Divide(hGENeta); 
  // isoEG efficiency
  TH1F *eratio2=new TH1F("eratio2","",25,0.,50.);
  eratio2=(TH1F*)hisoEGeta->Clone();
  eratio2->Divide(hGENeta);  
  
  TCanvas * can = new TCanvas();
  can->Divide(2,1);

  can->cd(1); // gPad->SetLogy();
  gPad->SetGridx(); gPad->SetGridy();

  TLegend *leg2 = new TLegend(0.17,0.67,0.52,0.90);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.04);
  leg2->SetTextFont(42);

  leg2->AddEntry(hGen,"gen e","L");
  leg2->AddEntry(hUpgrade1,"EG cand","L");
  leg2->AddEntry(hUpgrade2,"isoEG cand","L");

  hGen->Draw("HIST"); 
  hGen->SetLineWidth(2); hGen->SetLineColor(1);
  hGen->GetXaxis()->SetRangeUser(10.,50.);
  hUpgrade1->Draw("HISTSAME"); 
  hUpgrade1->SetLineColor(2); hUpgrade1->SetLineWidth(2);
  hGen->GetYaxis()->SetRangeUser(0.,2000.);
  hGen->SetTitle(";p_{T} (GeV);a.u.");
  hUpgrade2->Draw("HISTSAME"); hUpgrade2->SetLineWidth(2);
  hUpgrade2->SetLineColor(4);

  leg2->Draw("SAME");


  can->cd(2); // csctf_histo->Draw();
  
  gPad->SetGridx(); gPad->SetGridy();

  ratio1->Draw("P");
  ratio1->SetMarkerStyle(20);
  ratio1->SetMarkerColor(2);
  ratio1->GetYaxis()->SetRangeUser(0.,1.01);
  ratio1->SetTitle(";p_{T} (GeV);efficiency");
  
  //ratio2->Divide(hGen); 
  ratio2->Draw("PSAME");
  ratio2->SetMarkerStyle(20);
  ratio2->SetMarkerColor(4);

  
  TCanvas * can0 = new TCanvas();
  can0->Divide(2,1);

  can0->cd(1); // gPad->SetLogy();
  gPad->SetGridx(); gPad->SetGridy();

  TLegend *leg = new TLegend(0.17,0.67,0.52,0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);

  leg->AddEntry(hGENeta,"gen e","L");
  leg->AddEntry(hEGeta,"EG cand","L");
  leg->AddEntry(hisoEGeta,"isoEG cand","L");

  hGENeta->Draw("HIST"); 
  hGENeta->SetLineWidth(2); hGENeta->SetLineColor(1);
  hGENeta->GetXaxis()->SetRangeUser(10.,50.);
  hEGeta->Draw("HISTSAME"); 
  hEGeta->SetLineColor(2); hEGeta->SetLineWidth(2);
  hGENeta->GetYaxis()->SetRangeUser(0.,2000.);
  hGENeta->SetTitle(";#eta;a.u.");
  hisoEGeta->Draw("HISTSAME"); hisoEGeta->SetLineWidth(2);
  hisoEGeta->SetLineColor(4);

  leg->Draw("SAME");


  can0->cd(2); // csctf_histo->Draw();
  
  gPad->SetGridx(); gPad->SetGridy();

  eratio1->Draw("P");
  eratio1->SetMarkerStyle(20);
  eratio1->SetMarkerColor(2);
  eratio1->GetYaxis()->SetRangeUser(0.,1.01);
  eratio1->SetTitle(";p_{T} (GeV);efficiency");
  
  //eratio2->Divide(hGen); 
  eratio2->Draw("PSAME");
  eratio2->SetMarkerStyle(20);
  eratio2->SetMarkerColor(4);


  TCanvas * can2 = new TCanvas();
  gPad->SetLogy();

  DRmin1->Draw(); 
  DRmin1->SetLineColor(2); DRmin1->SetLineWidth(2);
  DRmin2->Draw("SAME");
  DRmin2->SetLineColor(4); DRmin2->SetLineWidth(2);
  DRmin1->GetYaxis()->SetRangeUser(0.1,50000.);

  gStyle->SetPalette(1);

  TCanvas * can3 = new TCanvas();
  can3->Divide(2,1);

  can3->cd(1); dPTvsdR1->Draw("COLZ");
  can3->cd(2); dPTvsdR2->Draw("COLZ");


  //hUpgrade->Divide(hGen);
  // hUpgrade->Draw();
  //  TGraphAsymmErrors *Eff = new TGraphAsymmErrors();
  // Eff->BayesDivide(hUpgrade,hGen);
  // Eff->Draw("AP");
}
