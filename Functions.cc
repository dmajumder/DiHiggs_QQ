#include <fstream>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include <TCanvas.h>
#include <TFile.h>
#include <TArray.h>
#include <TH1D.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TArray.h>
#include <TVector.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TString.h>
#include <TObject.h>
#include <TMath.h>
#include <TVector2.h>
#include <stdio.h>      /* Standard Library of Input and Output */
#include <complex.h>    /* Standard c++ Library of Complex Numbers */
#include "rootHistos.h" //declare the histos
#include "cuts.h" // basic cuts 
using namespace fastjet;
using namespace std;
/////////////////////////////////////////////////////////////////
void hello(){cout<<"\n\n\n HELLO!!!! \n\n"<<endl;}

bool GenLevel4b(vector<PseudoJet> jets, vector<int> btag, vector<int> fattag) {

  bool passGenCut(false) ; 

  //std::vector<fastjet::PseudoJet>const_iterator ijet ;

  int nhjets(0) ; 
  std::vector<fastjet::PseudoJet> hjets ; 

  for (int ijet = 0; ijet < jets.size(); ++ijet) {
    if (jets.at(ijet).pt() < 300 || abs(jets.at(ijet).eta()) > 2.4) continue ;
    if (btag.at(ijet) < 2) continue ; 
    //if (fattag.at(ijet) < 1) continue ;
    std::cout << " jet pt " << jets.at(ijet).pt() << " eta " << jets.at(ijet).eta() << " btag = " << btag.at(ijet) << std::endl ; 
    hjets.push_back(jets.at(ijet)) ; 
    ++nhjets ;  
  }

  std::cout << " nhjets = " << nhjets << std::endl ; 

  if (hjets.size() < 2) return passGenCut ; 

  double deta = abs(hjets.at(0).eta() - hjets.at(1).eta()) ; 
   
  if (deta < 1.2) return passGenCut ; 

  passGenCut = true ; 

  return passGenCut ; 
  
}

/////////////////////////////////////////////////////////////////
bool GenLevelDilep(vector<PseudoJet> jets,  vector<PseudoJet> taus){
    vector<PseudoJet> alltaus; 
    ///////////////////////////////////////////////////////////////////
    // gen level info // had == plus
    bool passGencut = false;
    //int blll=-1, bhhh=-1, ell=-1, ehh=-1, null=-1, nuhh=-1;
    //for(unsigned int nj1=0; nj1< particles.size(); nj1++) if(particles.at(nj1).user_index()==5) bhhh=nj1; else if(particles.at(nj1).user_index()==-5) blll=nj1; 
    //for(unsigned int nj1=0; nj1< leptons.size(); nj1++) {if(leptons.at(nj1).user_index() > 0) ell=nj1; else if(leptons.at(nj1).user_index() < 0) ehh=nj1;}
    //for(unsigned int nj1=0; nj1< neutrinos.size(); nj1++) {if(neutrinos.at(nj1).user_index() > 0) nuhh=nj1; else  null=nj1;}
    //PseudoJet lepTtrue, hadTtrue; int counttruth =-1;
    //if(blll!=-1 && bhhh!=-1 && ell!=-1 && ehh!=-1 && null!=-1 && nuhh!=-1){ 
    //cout<<"here "<<taus.size() <<endl;
    ////////////////////////////////////////////////////////////////////////////
    if(taus.size()>3) 
        if(taus.at(0).pt()>20 && taus.at(1).pt()>20 && taus.at(2).pt()>20 && taus.at(3).pt()>20) {  
            //alltaus = taus.at(0)+taus.at(1)+taus.at(2)+taus.at(3);
            double ptj1=-1, ptj2=-1 , etaj1=20 , etaj2=20;
            if (jets.size() > 0) {ptj1 = jets.at(0).pt(); etaj1 = jets.at(0).eta(); }
            if (jets.size() > 1) {ptj2 = jets.at(1).pt(); etaj2 = jets.at(1).eta(); }
            double hadtop[10]={(taus.at(0)+taus.at(1)+taus.at(2)+taus.at(3)).m(),
                ptj1,etaj1,ptj2,etaj2,0,0,0,0,0};
            for(unsigned i=0;i<10;i++) basicHadtop[i]->Fill(hadtop[i],1); 
            //genbhad->Fill(particles.at(bhhh).pt(),weight);
            //genblep->Fill(particles.at(blll).pt(),weight);
            //genbhadeta->Fill(particles.at(bhhh).eta(),weight);
            //genblepeta->Fill(particles.at(blll).eta(),weight);
            //cout<<"here"<<endl;
            passGencut = true;
        } else {
            double hadtop[10]={-10,-10,-10,-10,-10,0,0,0,0,0};
            
        }// close if cut
    //} // close if have enough particles
    /////////////////////////////////////////////////////////////////////////
    return passGencut;
} // close gen level cuts
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
int recol( vector<PseudoJet> jets,vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, double weight){
    // from the jet collection find the hadronic W
    // construct met from all the rest
    // do a isolation vector
    vector<double> LepIso; int nlep=0; int nlepsurvive=0; //cout<<"njet "<<jets.size()<<endl;
    unsigned int jsize = jets.size();
    for(unsigned int j = 0;j<leptons.size();j++) {
        for(unsigned int i = 0;i<jsize;i++) LepIso.push_back(leptons.at(j).delta_R(jets.at(i)));
        if(leptons.size()>1 && jsize >0) for(unsigned int i = 0;i<leptons.size();i++) if (i!=j) LepIso.push_back(leptons.at(j).delta_R(leptons.at(i)));
        double MinDRLep = TMath::LocMin(LepIso.size(), &LepIso[0]);
        if(1>0 && leptons.size()>1 && jsize >0
      //     && LepIso[MinDRLep] > lepiso  
      //     && leptons.at(j).pt()> ptlepton 
      //     && abs(leptons.at(j).eta())< etal
           ) nlep++; // close basic cuts
    } // close for nlep
    //cout<<nlep<<endl;
    Nlep_passing_kLooseID->Fill(nlep,1);
    return nlep;
}// close recow
/////////////////////////////////////////////////////////////////
int isbtagged(vector<PseudoJet> jets, vector<int> & btag, vector<int> & bmistag, vector<int> & btrue);
/////////////////////////////////////////////////////////////////
int recojets(vector<PseudoJet> particles,vector<PseudoJet> & jets_akt, vector<int> & btag, vector<int> & bmistag, vector<int> & fattag, vector<int> & btrue, double weight){
    JetDefinition akt(antikt_algorithm, RR);
    ClusterSequence cs_akt(particles, akt);
    //vector<PseudoJet> jets_akt;
    Selector jet_selector = SelectorPtMin(jet_ptmin) && SelectorAbsRapMax(rapmax);
    if(shower){jets_akt = sorted_by_pt(jet_selector(cs_akt.inclusive_jets()));} // first we do akt jets from particles
    else{double const ptmin=0.0; Selector jet_selector_parton = SelectorPtMin(ptmin);jets_akt = sorted_by_pt(jet_selector_parton(cs_akt.inclusive_jets()));}//if parton jet pt min zero
    //////////////////// do a jet cut
    vector<PseudoJet> jets_final;     unsigned int njets;
    for (unsigned int i = 0; i < jets_akt.size(); i++) if(jets_akt.at(i).pt()>jet_ptminfinal) jets_final.push_back(jets_akt.at(i));
    //for (unsigned int i = 0; i < particles.size(); i++) if(particles.at(i).pt()<jet_ptminfinal) {njets =0; return njets; } //jets_final.push_back(particles.at(i));
    njets = jets_final.size();
    int nbtag = isbtagged(jets_akt, btag, bmistag,btrue); // check wheather the b(c)jet is b--(mis)tagable
    // fill btags
    int count=0; for(unsigned int i = 0; i < njets; i++) if(btag[i]>0)count++; //btagselected->Fill(count,weight);
    ///////////////////// check tag
    JetDefinition CA10(cambridge_algorithm, Rsb);
    // Filter definition to improve mass resolution
    Filter filter(JetDefinition(cambridge_algorithm, Rfilt), SelectorNHardest(n_subjet));
    PseudoJet tagged_jet;
    for (unsigned int i = 0; i < njets; i++) { // to each akt jet
        // first recluster with some large CA (needed for mass-drop)
        ClusterSequence cs_tmp(jets_akt[i].constituents(), CA10);
        // next get hardest jet
        PseudoJet ca_jet = sorted_by_pt(cs_tmp.inclusive_jets())[0]; // find the cores
        // now run mass drop tagger
        MassDropTagger md_tagger(mu, ycut); // define the cut on mass drop
        // mu: ratio in between mass of cores, symetric splitting
        tagged_jet = md_tagger(ca_jet);
        if(tagged_jet.m()>10) {
            PseudoJet filtered_jet = filter(jets_akt.at(i)); // filter to tag
            if(filtered_jet.m()>Mfat) fattag.push_back(i); //see = 1; else see = 0; // no fat tag 
        } //else see = 0; 
    } // close find mass drop
    //jets = jets_akt;---
    Njets_passing_kLooseID->Fill(njets,weight);
    btagselected->Fill(nbtag,weight); 
    return njets;
} // close cluster jets
////////////////////////////////////////////////////////////////////////////////////////////////
int isbtagged(vector<PseudoJet> jets, vector<int> & btag, vector<int> & bmistag, vector<int> & btrue){ 
    unsigned int jsize = jets.size();
    int nbtag=0;
    for (unsigned int i=0; i<jsize; i++) { // check wheter jet have inside a b's are taggable 
        vector<PseudoJet> constitu=jets.at(i).constituents();
        unsigned int csize = constitu.size();
        int see=0,see2=0,seetruth=0,id;
        for (unsigned int j=0; j<csize; j++) {
            //cout<<"constituents flavour "<<constitu.at(j).user_index()<<endl;
            if(constitu.at(j).user_index() == 5 || constitu.at(j).user_index() == -5) {
                seetruth++; id=constitu.at(j).user_index();
                if(constitu.at(j).pt() > bjetpt && constitu.at(j).eta() < etab) {see++; nbtag++;}// to reco btag 
            }  
            //if( abs(constitu.at(j).user_index()) == 4  // work !!
            //     && constitu.at(j).pt() > bjetpt
            // && constitu.at(j).eta() < etab
            //  ) {see2++;}// bmistag.push_back(1);} bmistag.push_back(0);
        } // close constituents
        //bmistag.push_back(see2);
        //cout<<see<<endl;
        btag.push_back(see); //else btag.push_back(0); // count all tag/jet
        //if(see==0 && see2>0) bmistag.push_back(1); else bmistag.push_back(0); // count only one tag/jet
        btrue.push_back(seetruth); 
        //cout<<"b-quarks mistagged = " <<bmistag[i] <<" b-quark = " <<btag[i] <<endl;
    } // close for each jet
    //int numbb=0; for(unsigned i=0 ; i< btrue.size() ; i++ ) numbb += btrue[i];
    //cout<<" total "<<numbb<<endl;
    return nbtag;
} // close isbtagged
/////////////////////////////////////////////////////////////////////////
// save the histos
int save_hist(int isample){
    const char* Mass;
    Mass = Form("Control_parton_m300_%d.root",isample); 
    cout<<" saved "<< Form("Control_parton_m4000_%d.root",isample)<<endl;
    TFile f1(Mass, "recreate");
    f1.cd();
    Njets_passing_kLooseID->Write();
    Nlep_passing_kLooseID->Write();
    btagselected->Write();
    for(unsigned i=0;i<10;i++) basicHadtop[i]->Write();
    f1.Close();
    Njets_passing_kLooseID->Reset();
    Nlep_passing_kLooseID->Reset();
    btagselected->Reset();

    //  basicLeptons[0]->Reset();
    //  basicLeptons[1]->Reset();
    //  basicLeptons[2]->Reset();
    //  for(unsigned i=0;i<10;i++) basicHadtop[i]->Reset();
    //  for(unsigned i=0;i<13;i++) basicLeptop[i]->Reset();

    return 0;
}
/////////////////////////////////////////////////////////////////////////
int decla(int mass){
    
    delete gDirectory->FindObject("leptop1");
    delete gDirectory->FindObject("hadtop1");
    delete gDirectory->FindObject("njets_passing_kLooseID_ct4");
    delete gDirectory->FindObject("btagselected");
    delete gDirectory->FindObject("E1histpt");
    delete gDirectory->FindObject("E1histeta");
    delete gDirectory->FindObject("MetMass_ct4");
    delete gDirectory->FindObject("H1hist");
    delete gDirectory->FindObject("H1histpt");
    delete gDirectory->FindObject("H1histeta");
    delete gDirectory->FindObject("H1histphi");
    delete gDirectory->FindObject("HW1hist");
    delete gDirectory->FindObject("HW1histpt");
    delete gDirectory->FindObject("HW1histeta");
    delete gDirectory->FindObject("HW1histphi");
    delete gDirectory->FindObject("recotruth");
    delete gDirectory->FindObject("detabb");
    delete gDirectory->FindObject("H1LepThist");
    delete gDirectory->FindObject("H1LepThistpt");
    delete gDirectory->FindObject("H1LepThisteta");
    delete gDirectory->FindObject("H1LepThistphi");
    delete gDirectory->FindObject("HW1LepThist");
    delete gDirectory->FindObject("HW1LepThistpt");
    delete gDirectory->FindObject("HW1LepThisteta");
    delete gDirectory->FindObject("HW1LepThistphi");
    delete gDirectory->FindObject("pnuzerror");
    delete gDirectory->FindObject("recotruthlept");
    delete gDirectory->FindObject("mterror");
    delete gDirectory->FindObject("wmt");
    delete gDirectory->FindObject("tmt");
    delete gDirectory->FindObject("detalb");
    delete gDirectory->FindObject("mraz");
    delete gDirectory->FindObject("mrazt");
    delete gDirectory->FindObject("razratio");
    
    
    const char* label="without btag im reco";
    
    Njets_passing_kLooseID = new TH1D("njets_passing_kLooseID_ct4",  
                                      label, 
                                      13, -0.5, 12.5);
    Njets_passing_kLooseID->GetYaxis()->SetTitle("");
    Njets_passing_kLooseID->GetXaxis()->SetTitle("Njets after showering"); 

    Nlep_passing_kLooseID = new TH1D("nlep_passing_kLooseID_ct4",  
                                      label, 
                                      13, -0.5, 12.5);
    Nlep_passing_kLooseID->GetYaxis()->SetTitle("");
    Nlep_passing_kLooseID->GetXaxis()->SetTitle("Nleptons after showering"); 
    
    btagselected = new TH1D("btagselected",  
                            label, 
                            13, -0.5, 12.5);
    btagselected->GetYaxis()->SetTitle("");
    btagselected->GetXaxis()->SetTitle("b-tagable b's on selected events");
    

    ///////////////////////////////////////////////////////////////////////////
    // for hadronic tops
    TH1D *H1hist = new TH1D("H1hist",  
                            label, 
                            200, 0, 2000);
    H1hist->GetYaxis()->SetTitle("Events/ 2 GeV");
    H1hist->GetXaxis()->SetTitle("mass_{4 #tau} (GeV)");
    basicHadtop.push_back (H1hist); 
    
    TH1D *H1histpt = new TH1D("H1histpt",  
                              label, 
                              100, 0, 1000);
    H1histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
    H1histpt->GetXaxis()->SetTitle("j_1 P_T (GeV)");
    basicHadtop.push_back (H1histpt); 
    
    TH1D *H1histeta = new TH1D("H1histeta",  
                               label, 
                               30, -6, 6);
    H1histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
    H1histeta->GetXaxis()->SetTitle("j_1 #eta");
    basicHadtop.push_back (H1histeta); 
    
    TH1D *H1histphi = new TH1D("H1histphi",  
                               label, 
                               100, 0, 1000);
    H1histphi->GetYaxis()->SetTitle("Events/ 2 GeV");
    H1histphi->GetXaxis()->SetTitle("j_2 P_T  (GeV)");
    basicHadtop.push_back (H1histphi); 
    // had w
    TH1D *HW1hist = new TH1D("HW1hist",  
                             label, 
                             30, -6, 6);
    HW1hist->GetYaxis()->SetTitle("Events/ 2 GeV");
    HW1hist->GetXaxis()->SetTitle("j_2 #eta");
    basicHadtop.push_back (HW1hist); 
    
    TH1D *HW1histpt = new TH1D("HW1histpt",  
                               label, 
                               100, 0, 600);
    HW1histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
    HW1histpt->GetXaxis()->SetTitle("W_{had} P_T (GeV)");
    basicHadtop.push_back (HW1histpt); 
    
    TH1D *HW1histeta = new TH1D("HW1histeta",  
                                label, 
                                30, -6, 6);
    HW1histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
    HW1histeta->GetXaxis()->SetTitle("#eta_{W_{had}} (GeV)");
    basicHadtop.push_back (HW1histeta); 
    
    TH1D *HW1histphi = new TH1D("HW1histphi",  
                                label, 
                                30, 0, 5);
    HW1histphi->GetYaxis()->SetTitle("Events/ 2 GeV");
    HW1histphi->GetXaxis()->SetTitle("#phi_{W_{had}} (GeV)");
    basicHadtop.push_back (HW1histphi); 
    
    TH1D *recotruth = new TH1D("recotruth",  
                               label, 
                               3, -0.5, 2.5);
    recotruth->GetYaxis()->SetTitle("Events/ 2 GeV");
    recotruth->GetXaxis()->SetTitle("reco truth");
    basicHadtop.push_back (recotruth); 
    
    TH1D *detabb = new TH1D("detabb",  
                            label, 
                            30, 0., 5);
    detabb->GetYaxis()->SetTitle("Events/ 2 GeV");
    detabb->GetXaxis()->SetTitle("#Delta#eta bb");
    basicHadtop.push_back (detabb); 
    ///////////////////////////////////////////////////////////////////////////
    // for leptonic tops
    TH1D *H1LepThist = new TH1D("H1LepThist",  
                                label, 
                                70, 0, 1000);
    H1LepThist->GetYaxis()->SetTitle("Events/ 2 GeV");
    H1LepThist->GetXaxis()->SetTitle("mass t_{lep} (GeV)");
    basicLeptop.push_back (H1LepThist); 
    
    TH1D *H1LepThistpt = new TH1D("H1LepThistpt",  
                                  label, 
                                  20, 0, 300);
    H1LepThistpt->GetYaxis()->SetTitle("Events/ 2 GeV");
    H1LepThistpt->GetXaxis()->SetTitle("t_{lep} P_T (GeV)");
    basicLeptop.push_back (H1LepThistpt); 
    
    TH1D *H1LepThisteta = new TH1D("H1LepThisteta",  
                                   label, 
                                   30, -6, 6);
    H1LepThisteta->GetYaxis()->SetTitle("Events/ 2 GeV");
    H1LepThisteta->GetXaxis()->SetTitle("#eta_{t_{lep}} (GeV)");
    basicLeptop.push_back (H1LepThisteta); 
    
    TH1D *H1LepThistphi = new TH1D("H1LepThistphi",  
                                   label, 
                                   30, 0, 5);
    H1LepThistphi->GetYaxis()->SetTitle("Events/ 2 GeV");
    H1LepThistphi->GetXaxis()->SetTitle("#phi_{t_{lep}} (GeV)");
    basicLeptop.push_back (H1LepThistphi); 
    // had w
    TH1D *HW1LepThist = new TH1D("HW1LepThist",  
                                 label, 
                                 100, 80, 81);
    HW1LepThist->GetYaxis()->SetTitle("Events/ 2 GeV");
    HW1LepThist->GetXaxis()->SetTitle("mass W_{lep} (GeV)");
    basicLeptop.push_back (HW1LepThist); 
    
    TH1D *HW1LepThistpt = new TH1D("HW1LepThistpt",  
                                   label, 
                                   20, 0, 300);
    HW1LepThistpt->GetYaxis()->SetTitle("Events/ 2 GeV");
    HW1LepThistpt->GetXaxis()->SetTitle("W_{lep} P_{T} (GeV)"); 
    basicLeptop.push_back (HW1LepThistpt); 
    
    TH1D *HW1LepThisteta = new TH1D("HW1LepThisteta",  
                                    label, 
                                    30, -6, 6);
    HW1LepThisteta->GetYaxis()->SetTitle("Events/ 2 GeV");
    HW1LepThisteta->GetXaxis()->SetTitle("#eta_{W_{lep}} (GeV)");
    basicLeptop.push_back (HW1LepThisteta); 
    
    TH1D *HW1LepThistphi = new TH1D("HW1LepThistphi",  
                                    label, 
                                    30, 0, 5);
    HW1LepThistphi->GetYaxis()->SetTitle("Events/ 2 GeV");
    HW1LepThistphi->GetXaxis()->SetTitle("#phi_{W_{lep}} (GeV)");
    basicLeptop.push_back (HW1LepThistphi); 
    
    TH1D *pnuzerror = new TH1D("pnuzerror",  
                               label, 
                               120, -60, 60);
    pnuzerror->GetYaxis()->SetTitle("Events/ 2 GeV");
    pnuzerror->GetXaxis()->SetTitle("mW(reco) - mW(truth)");
    basicLeptop.push_back (pnuzerror); 
    
    TH1D *recotruthlept = new TH1D("recotruthlept",  
                                   label, 
                                   3, -0.5, 2.5);
    recotruthlept->GetYaxis()->SetTitle("Events/ 2 GeV");
    recotruthlept->GetXaxis()->SetTitle("reco truth lept");
    basicLeptop.push_back (recotruthlept); 
    
    TH1D *mterror = new TH1D("mterror",  
                             label, 
                             120, -60, 60);
    mterror->GetYaxis()->SetTitle("Events/ 2 GeV");
    mterror->GetXaxis()->SetTitle("mt(reco) - mt(truth)");
    basicLeptop.push_back (mterror); 
    
    TH1D *wmt = new TH1D("wmt",  
                         label, 
                         75, 0, 150);
    wmt->GetYaxis()->SetTitle("Events/ 2 GeV");
    wmt->GetXaxis()->SetTitle("W transverse mass");
    basicLeptop.push_back (wmt); 
    
    TH1D *tmt = new TH1D("tmt",  
                         label, 
                         50, 0, 600);
    tmt->GetYaxis()->SetTitle("Events/ 2 GeV");
    tmt->GetXaxis()->SetTitle("Top transverse mass (Wb)");
    basicLeptop.push_back (tmt); 
    
    TH1D *detalb = new TH1D("detalb",  
                            label, 
                            30, 0, 6);
    detalb->GetYaxis()->SetTitle("Events/ 2 GeV");
    detalb->GetXaxis()->SetTitle("#Delta#eta b l");
    basicLeptop.push_back (detalb); 
    
    TH1D *mraz = new TH1D("mraz",  
                          label, 
                          150, 0, 1000);
    mraz->GetYaxis()->SetTitle("Events/ 2 GeV");
    mraz->GetXaxis()->SetTitle("MR");
    basicLeptop.push_back (mraz); 
    
    TH1D *mrazt = new TH1D("mrazt",  
                           label, 
                           150, 0, 300);
    mrazt->GetYaxis()->SetTitle("Events/ 2 GeV");
    mrazt->GetXaxis()->SetTitle("MRt");
    basicLeptop.push_back (mrazt); 
    
    TH1D *razratio = new TH1D("razratio",  
                              label, 
                              100, 0, 1);
    razratio->GetYaxis()->SetTitle("Events/ 2 GeV");
    razratio->GetXaxis()->SetTitle("MRt/MR");
    basicLeptop.push_back (razratio); 
    
    TH1D *tmtal = new TH1D("tmtal",  
                           label, 
                           50, 0, 600);
    tmtal->GetYaxis()->SetTitle("Events/ 2 GeV");
    tmtal->GetXaxis()->SetTitle("top transverse mass (lb - #nu)");
    basicLeptop.push_back (tmtal); 
    
    TH1D *mbl = new TH1D("massbl",  
                         label, 
                         30, 0, 1000);
    mbl->GetYaxis()->SetTitle("Events/ 2 GeV");
    mbl->GetXaxis()->SetTitle("m_{lb}");
    basicLeptop.push_back (mbl); 
    ///////////////////////////////////////////////////////////////////////////////////
    return 0;
}


/*
 */ 
/*  /////////////////////////////////////////////////////////////////////
 // find w -- closest hadronic mass
 vector<int> j1,j2; vector<double> a3; int minM;
 for(unsigned int nj1=0; nj1< jsize; nj1++) 
 for(unsigned int nj2=nj1+1; nj2< jsize; nj2++) 
 if(btag[nj1]==0 && btag[nj2]==0)
 { 
 //std::cout<<nj1<<nj2<<" "<<nj3<<nj4<<std::endl; 
 double invmassA =  (jets.at(nj1)+jets.at(nj2)).m();
 a3.push_back((invmassA-wmass)*(invmassA-wmass)); 
 j1.push_back(nj1); j2.push_back(nj2); 
 // we also what to keep the nj...           
 } // loop on jets  
 minM = TMath::LocMin(a3.size(), &a3[0]);
 PseudoJet hadW = jets[j1[minM]] + jets[j2[minM]];
 //cout<<" the hadronic w's are: "<<j1[minM]<<" "<<j2[minM]<<" "<<btrue[j1[minM]]<<" "<<btrue[j2[minM]]<<endl;
 // find t -- closest jet
 unsigned wjet1=j1[minM], wjet2=j2[minM]; 
 vector<int> j3; vector<double> a4; int minMt;
 for(unsigned int nj1=0; nj1< jsize; nj1++) 
 if(nj1!=wjet1 && nj1!=wjet2) 
 if(btag[nj1]>0) 
 { 
 double invmassA =  (jets.at(nj1)+hadW).m();
 a4.push_back((invmassA-tmass)*(invmassA-tmass)); j3.push_back(nj1); 
 //double dr=  jets.at(nj1).delta_R(hadW);
 //a4.push_back(dr); j3.push_back(nj1);
 //cout<<jsize<<" "<<btag.size()<<" "<<btag[nj1]<<endl;
 } // loop on jets
 //if(a4.size()>0){
 minMt = TMath::LocMin(a4.size(), &a4[0]); bh = j3[minMt];
 //cout<<" the hadronic t's are: "<<j1[minM]<<" "<<j2[minM]<<" "<<j3[minMt]<<endl;
 PseudoJet hadt = jets[bh] + hadW;
 */
/////////////////////////////////////////////////
//TVector2 teste = jets.at(j1[minM]).pt()+jets.at(j2[minM]).pt();
//cout<<" teste "<<
//(jets.at(j1[minM]).px()+jets.at(j2[minM]).px())*(jets.at(j1[minM]).px()+jets.at(j2[minM]).px())+
//(jets.at(j1[minM]).py()+jets.at(j2[minM]).py())*(jets.at(j1[minM]).py()+jets.at(j2[minM]).py())<<" "<<
//(jets.at(j1[minM]).pt()+jets.at(j2[minM]).pt())*(jets.at(j1[minM]).pt()+jets.at(j2[minM]).pt())<<" "<<
// teste*teste<<" "<<endl;
