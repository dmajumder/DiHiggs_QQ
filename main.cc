//////////////////
// to run:
// source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.00/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh
// make HH_VBF
// ./HH_VBF
/////////////////
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include <fstream>
#include <iostream>
#include <vector>
using namespace fastjet;
using namespace std;
#include "Functions.h"
#include "choices.h"
int main() {
    srand( time(NULL) );
    //////////////////////////////////////////////////
    // input
    //string fileparton = "Wt_5.lhe.decayed"; 
     string fileparton = "/Users/Xanda/cernbox/Topone_higgs/Bp_signal/MR_700_on.lhe.shower" ; 
    cout<<"\n\n reading file = "<<fileparton<<endl;
    // open the file
    ifstream in1; in1.open(fileparton.c_str());
    // loop in events
    decla(0);  
    for(unsigned int ievent=0;ievent<10000;ievent++){ // read and process for each event ===> input the number of events
      unsigned int counter=0,countert=0,counterl=0,countern=0, countertau  = 0,  nb = 0; // conter to have sure lates
      double Px, Py , Pz, E; int pID; unsigned int nparticles;
      vector<PseudoJet> particles; //jets 
      vector<PseudoJet> neutrinos;
      vector<PseudoJet> leptons; 
      vector<PseudoJet> tops;   
      vector<PseudoJet> taus;  
      // + gammas .... ....
      /////////////////////////////////////////////////////////////////////
      // read and understand
      string c; in1>>c; // read a line
      in1>>nparticles; // recognize the first line, have only one number that is the number of particles in the event
      for(unsigned int ipart=0;ipart<nparticles;ipart++){ // loop on particles
         in1 >> pID >> Px >> Py >> Pz >> E ;//>> idup;
         if(abs(pID) < 6 || pID==21){particles.push_back(fastjet::PseudoJet(Px,Py,Pz,E));particles.at(counter).set_user_index(pID); if(abs(pID) == 5) nb++; counter++;} 
         else if (abs(pID)==6) {tops.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); countert++;} 
         else if (abs(pID)==11 || abs(pID)==13) {leptons.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); leptons.at(counterl).set_user_index(pID); counterl++;} 
         else if (abs(pID)==12 || abs(pID)==14) {neutrinos.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); neutrinos.at(countern).set_user_index(pID); countern++;} // close if 
         else if (abs(pID)==15) {taus.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); taus.at(countertau).set_user_index(pID); countertau++;} // close if 
      } // close for each particle
      vector<PseudoJet> jets; vector<int> btag, bmistag, fattag, btrue; int bh,bl; double met=0;
      bool lepcuts=false, lepwreco=false; int hadtopreco;
      int njets = recojets(particles, jets,btag,bmistag,fattag,btrue,1); 
      bool ana =  GenLevelDilep(jets, taus );
      bool ana4b =  GenLevel4b(jets, btag, fattag);
      //cout<<" number of isolated leptons: "<< nlep <<endl<<endl; 
    } 
    save_hist(1); 
    in1.close(); // close for file
    return 0;
} //
