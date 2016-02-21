#include <TH1D.h>
#include <TH2D.h>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

vector<TH1D *> basicGen;
vector<TH1D *> basicLeptons;
vector<TH1D *> basicHadtop;
vector<TH1D *> basicLeptop;
//////////
vector<TH1D*>  Njets_passing_kLooseID;
vector<TH1D*>  Nlep_passing_kLooseID;
vector<TH1D*>  btagselected;
vector<TH1D*>  leptop;
vector<TH1D*>  hadtop;
vector<TH1D*>  fullmhh;
vector<TH1D*>  fullpth1;
vector<TH1D*>  fullpth2;
vector<TH1D*>  fullcost;

vector<TH1D*>  genmbb; 
TH1D *genmbb1;
vector<TH1D*>  genpth;
vector<TH1D*>  genpthsub;
vector<TH1D*>  gencost;
vector<TH1D*>  gendeta;
TH2D*  bin1;
TH2D*  bin1re;

TH2D* histoV1_bin1;
TH2D*   histoJHEP_bin1;
vector<TH1D*>   histoJHEP_mhh;
vector<TH1D*>   histoJHEP_pth;
vector<TH1D*>   histoJHEP_cost;
 
TH2D*   histoOutlayers;

TH2D*  JHEP2D;
vector<TH1D*>  JHEPmhh;
vector<TH1D*>  JHEPcost;
vector<TH1D*>  JHEPpt;

vector<TH1D*>  REmhh;
vector<TH1D*>  REcost;
vector<TH1D*>  REpt;


//////////////////////////////////



const char* JHEPbench[12]={
    "260",
    "350",
    "450",
    "500",
    "1000",
    "1500",
    "2000",
    "2500"};

const char* JHEP2Dclu[12]={
    "JHEP2D_cluster0.png","JHEP2D_cluster1.png","JHEP2D_cluster2.png",
    "JHEP2D_cluster3.png","JHEP2D_cluster4.png","JHEP2D_cluster5.png",
    "JHEP2D_cluster6.png","JHEP2D_cluster7.png","JHEP2D_cluster8.png",
    "JHEP2D_cluster9.png","JHEP2D_cluster10.png","JHEP2D_cluster11.png"};
const char* JHEP2Dmhhclu[12]={
    "JHEPmhh_cluster0.png","JHEPmhh_cluster1.png","JHEPmhh_cluster2.png",
    "JHEPmhh_cluster3.png","JHEPmhh_cluster4.png","JHEPmhh_cluster5.png",
    "JHEPmhh_cluster6.png","JHEPmhh_cluster7.png","JHEPmhh_cluster8.png",
    "JHEPmhh_cluster9.png","JHEPmhh_cluster10.png","JHEPmhh_cluster11.png"};
const char* JHEP2Dptclu[12]={
    "JHEPpt_cluster0.png","JHEPpt_cluster1.png","JHEPpt_cluster2.png",
    "JHEPpt_cluster3.png","JHEPpt_cluster4.png","JHEPpt_cluster5.png",
    "JHEPpt_cluster6.png","JHEPpt_cluster7.png","JHEPpt_cluster8.png",
    "JHEPpt_cluster9.png","JHEPpt_cluster10.png","JHEPpt_cluster11.png"};
const char* JHEP2Dcostclu[12]={
    "JHEPcost_cluster0.png","JHEPcost_cluster1.png","JHEPcost_cluster2.png",
    "JHEPcost_cluster3.png","JHEPcost_cluster4.png","JHEPcost_cluster5.png",
    "JHEPcost_cluster6.png","JHEPcost_cluster7.png","JHEPcost_cluster8.png",
    "JHEPcost_cluster9.png","JHEPcost_cluster10.png","JHEPcost_cluster11.png"};


TH1D *genbhad;
TH1D *genblep;
TH1D *genbhadeta;
TH1D *genblepeta;