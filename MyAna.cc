//////////////////
// to run:
// source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.00/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh
// make HH_VBF
// ./HH_VBF
/////////////////
/*
 To macht gen numbering scheeme
 
 mg5>define bb = b b~
 Defined multiparticle bb = b b~
 mg5>define ww = w+ w-
 Defined multiparticle ww = w+ w-
 mg5>define ll = e+ e- mu+ mu-
 Defined multiparticle ll = e- mu- e+ mu+
 mg5>define nu = vl vl~
 Defined multiparticle nu = ve vm vt ve~ vm~ vt~
 mg5>define j = g u c d s u~ c~ d~ s~ b b~
 Defined multiparticle j = g u c d s u~ c~ d~ s~ b b~
 mg5>generate p p > t t~, (tt > ww bb, w+ > ll nu, w- > j j )
 Command "generate p p > t t~, (tt > ww bb, w+ > ll nu, w- > j j )" interrupted with error:
 InvalidCmd : No particle tt in model
 mg5>define tt = t t~
 Defined multiparticle tt = t t~
 mg5>generate p p > t t~, (tt > ww bb, w+ > ll nu, w- > j j ) on-on: 13 diagrams ?
 generate p p > t w- b~, w- > ll nu, (t > w+ b, w+ > j j) on-off:  21 diagrams
 generate p p > t~ w+ b, w+ > j j , (t~ > w- b~, w- > ll nu) off-on: 21 diagrams
 */
//////////////////////////////////////////////////////////
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "TLorentzVector.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TStyle.h>
using namespace fastjet;
using namespace std;
#include "Functions.h"
#include "choices.h"
int main() {
    srand( time(NULL) );
    hello();
    string fileparton , fileshower , dataparton = ".decayed" , datashower;
    if(showering) datashower = ".shower"; else datashower = ".decayed"; 
    //////////////////////////////////////////////////
    // input
    // we are going to have only four samples, with 180 GeV mt selection, merge them and separate
    string path[1]={"/Users/Xanda/Documents/programs/VLQ_h_light/TT_hh/Graviton_Parton/MGraviton_"};
    string sample[8] = { 
        "260.lhe",
        "350.lhe",
        "450.lhe",
        "800.lhe",
        "1000.lhe",
        "1500.lhe",
        "2000.lhe",
        "2500.lhe"

    };// width
    /////////////////////////////////////////////
    // here we put the options for ploting
    style ();
    ///////////////////////////////////////////////////////// 
    // vector to make a table
    /////////////////////////////////////////////////////////
    double nevents =5000; // /fb
    //////////////////////////////////////////////////////////////////////////////////
    // to each cut deffinition and each one of the four region deffintions I pass by the four files and then save 
    //int cluster = 10;       // cluster number
    //for(unsigned int outlayer=0; outlayer<6; outlayer++) { // six outlayers
    double listweight[12];
    // to efficiencies
    // vector<vector<int>> A(dimension, vector<int>(dimension));
    double total[12];
    double selected[12];
    double selected4b[12];
    for(unsigned int isample=0; isample<8; isample++) {
        total[isample]=0;
        selected[isample]=0;
        selected4b[isample]=0;
    }
    decla(0,0); // sum all samples | take the histo for one cluster
    for(unsigned int isample=0; isample<8; isample++)  { // samples ...    
        
        //////////////////////////////////////////////
        fileshower = path[0] + sample[isample] +  datashower; cout<<"\n\n reading file = "<<fileshower<<endl;
        ifstream in1shower; in1shower.open(fileshower.c_str());
        // counters to categories
        int mhh1=0, mhh2=0 , mhh3=0;
        for(unsigned int ievent=0;ievent<nevents;ievent++){ // read and process for each event 
            double Px, Py , Pz, E; int pID, mother; unsigned int nparticles;
            unsigned int counter=0,countert=0,counterl=0,countern=0, nb = 0;
            vector<PseudoJet> particles; //jets 
            vector<PseudoJet> neutrinos;
            vector<PseudoJet> leptons; 
            vector<PseudoJet> tops; 
            vector<PseudoJet> higgses; // the two first are the gen level ones 
            /////////////////////////////////////////////////////////////////////
            // read and understand
            string cshower; in1shower>>cshower; in1shower>>nparticles; 
            // read showered
            for(unsigned int ipart=0;ipart<nparticles;ipart++){ // loop on particles
                in1shower >> pID >> Px >> Py >> Pz >> E ;//>> idup;
                if(abs(pID) < 6 || pID==21){particles.push_back(fastjet::PseudoJet(Px,Py,Pz,E));particles.at(counter).set_user_index(pID); if(abs(pID) == 5) nb++; counter++;} 
                else if (abs(pID)==6) {tops.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); countert++;} 
                else if (abs(pID)==25) {higgses.push_back(fastjet::PseudoJet(Px,Py,Pz,E));} 
                else if (abs(pID)==11 || abs(pID)==13) {leptons.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); leptons.at(counterl).set_user_index(pID); counterl++;} 
                else if (abs(pID)==12 || abs(pID)==14) {neutrinos.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); neutrinos.at(countern).set_user_index(pID); countern++;} // close if 
            } // close for each particle
            //////////////////////////////////////////////////////////////////
           // for(unsigned int outlayer=0; outlayer<6; outlayer++) for(unsigned int clusterr=0; clusterr<12; clusterr++) {
                //cout<<"event "<< ievent<<" here cluster "<<clusterr<<" here out "<< outlayer<<endl;
            //cout<<higgses.size()<<endl;
                double weight =  GenLevelWeight(higgses,isample,0);
                //for(unsigned int outlayer=0; outlayer<6; outlayer++) for(unsigned int clusterr=2; clusterr<3; clusterr++) {
                //    cout<<"here cluster "<<clusterr<<" here out "<< outlayer<<endl;
                // { // clusters
                // begin analysis
                
                //cout<<"the return "<</Users/Xanda/Documents/codes/git/DiHiggs_QQ/MyAna.ccweight<<endl;
                if( weight > 0){
                    total[isample]+=weight;
                    vector<PseudoJet> jets; vector<int> btag, bmistag, fattag, btrue; int bh,bl; double met=0;
                    bool lepcuts=false, lepwreco=false; int hadtopreco;
                    
                    int njets = recojets(particles, jets,btag,bmistag,fattag,btrue,weight,isample,0); // jet deffinition and basic cuts
                    int hh4b =  HH4b( jets, btag, btrue, weight,isample,0);
                    if(njets >3) selected[isample]+=weight;
                    if(hh4b >  3) selected4b[isample]+=weight;
                    //int nlep =recol(jets,leptons,neutrinos,weightV1); // lepton isolation and basic cuts 
                    //listweight[clusterr]=weight;
                } // close cuts
                
           // }//close for outlayers
            //load();
            //} // close for cluster
            
        } in1shower.close();// in1.close(); // close for each event
        /////////////////////////////////////////////////////////////////////

    } // close to isample
    
    //cout<<"here"<<endl;
    //} // close to clusters
    //save_hist(0); // summ of all samples
    draw();
    //draw_out();
    //for(unsigned int cluster=0; cluster<12; cluster++)cout << "counting: "<<total[cluster][0]<<" "<<selected[cluster][0]<<" "<<selected4b[cluster][0]<<endl;
    //cout<<endl;
    //for(unsigned int cluster=0; cluster<12; cluster++){
    //    cout << "efficiency: ";//<<selected[cluster][0]/total[cluster][0]<<" "<<selected4b[cluster][0]/total[cluster][0];
    //    for(unsigned int outlayer=0; outlayer<7; outlayer++) cout << selected[cluster][outlayer]/total[cluster][outlayer]<<" "<<selected4b[cluster][outlayer]/total[cluster][outlayer]<<" ";
//        cout<<endl;        
    //for(unsigned int cluster=0; cluster<2; cluster++) cout<<cluster+1<<" "<<totalV.at(cluster)<<" "<<selectedV.at(cluster)<<" "<<selected4bV.at(cluster)<<endl;
    //cout<<endl<<" in efficiencies: "<<endl;
    //for(unsigned int cluster=0; cluster<2; cluster++) cout<<cluster+1<<" "<<selectedV.at(cluster)/totalV.at(cluster)<<" "<<selected4bV.at(cluster)/totalV.at(cluster)<<endl;
    // draw all together!!!
    
} // close to main
