#include <fstream>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include <TCanvas.h>
#include <TFile.h>
#include <TArray.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TArray.h>
#include <TVector.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TString.h>
#include <TObject.h>
#include <TMath.h>
#include <TStyle.h>
#include "TLorentzVector.h"
#include <TVector2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <stdio.h>   
/* Standard Library of Input and Output */
#include <complex.h>    /* Standard c++ Library of Complex Numbers */
#include "rootHistos.h" //declare the histos
#include "cuts.h" // basic cuts 
using namespace fastjet;
using namespace std;
/////////////////////////////////////////////////////////////////
void hello(){cout<<"\n\n\n HELLO!!!! \n\n"<<endl;}
/////////////////////////////////////////////////////////////////
double GenLevelWeight( vector<PseudoJet>  higgses , vector<PseudoJet>  genjet ,vector<PseudoJet>  topo ,    int cluster , int outlayer ){ 
    double weight;
    //cout<<"cluster "<<cluster<<" outlayer "<<outlayer<<endl;
    //cout<<" here "<<histoJHEP_bin1.size()<<" "<<histoJHEP_bin1[cluster].size()<<endl;
    if(higgses.size()>1){
        // only the first tow are the gen level ones
        PseudoJet higgs0 = higgses.at(0); 
        PseudoJet higgs1 = higgses.at(1);
        //if(higgs1.px() == -higgs0.px()) {
            // compute costheta star         
            TLorentzVector P1boost; P1boost.SetPxPyPzE(higgs1.px(),higgs1.py(),higgs1.pz(),higgs1.e());//  higgs1;
            TLorentzVector P12; P12.SetPxPyPzE((higgs1 +higgs0).px(),(higgs1 +higgs0).py(),(higgs1 +higgs0).pz(),(higgs1 +higgs0).e());//  higgs1;
            P1boost.Boost(-P12.BoostVector());                     
            Double_t costhetast = P1boost.CosTheta();
            //cout<<"here 2 "<<costhetast<<endl;
            Double_t mhh = (higgs0+higgs1).m();
            Double_t pth0 = (higgs0).pt(); 
            Double_t pth1 = (higgs1).pt();// sort by pt
            Double_t pth, pthsub;
            if (pth1>pth0) {pth = pth1; pthsub=pth0;} else {pth = pth0; pthsub=pth1;}
        
            weight = 1;
            //cout<<"here "<<weight<< " "<<mhh  <<endl;
            //genmbb1->Fill(mhh,weight);
            if(topo.size()>0)toponept[cluster]->Fill(topo.at(0).pt(),weight);
            genmbb[cluster]->Fill(mhh,weight); 
            genpth[cluster]->Fill(pth,weight);
            genpthsub[cluster]->Fill(pthsub,weight);
            gencost[cluster]->Fill(higgs0.eta());
            gencost2[cluster]->Fill(abs(costhetast),weight); //higgs1.eta());//
            gendeta[cluster]->Fill(abs(higgs0.eta()-higgs1.eta()));
            //gendeta[cluster]->Fill(higgs0.delta_R(higgs1));
            if(genjet.size()>0) {fullpth1[cluster]->Fill(genjet.at(0).pt());  leptop[cluster]->Fill(genjet.at(0).eta());}
            if(genjet.size()>1) {fullpth2[cluster]->Fill(genjet.at(1).pt()); hadtop[cluster]->Fill(genjet.at(0).eta());}
            //cout<<"here "<<endl;
            //bin1[cluster]->Fill(mhh,costhetast);
            //bin1re[cluster]->Fill(mhh,costhetast,weight);  
            //cout<<"here "<<endl;
            /*
            Int_t binmhhV1 = histoV1_bin1->FindBin(mhh,costhetast); 
            double binV1 = histoV1_bin1->GetBinContent(binmhhV1);
            
            TH2D *dumb = (TH2D*) histoJHEP_bin1.at(cluster)->Clone();
            Int_t binmhhJHEP = dumb->FindBin(mhh,costhetast); 
            double binJHEP = dumb->GetBinContent(binmhhJHEP);// binmhhJHEP;//
            
            //cout<<"here in out "<<cluster<<" bin jhep "<< binmhhJHEP<<endl;
            if (binV1 < 1 || binV1 == NULL) binV1 =0.001;
            weight = ((double)binJHEP)/((double)binV1);//
            if(binJHEP < 1 ) weight = 0;
            //cout<<"here"<<endl;
            //cout<<mhh<<" "<<costhetast<<" V1 "<< binV1 <<" JHEP "<< binJHEP <<" ratio "<<weight<<endl;
            //bool GenLevel =  GenLevelDilep(mhh,costhetast,pth,weight); //cout<<"GenLevel "<<GenLevel<<endl;
            genmbb.at(cluster)->Fill(mhh,weight); 
            genpth.at(cluster)->Fill(pth,weight); 
            gencost.at(cluster)->Fill(abs(costhetast),weight);
            // to sum all 2D to make rweight
//            bin1->Fill(pth,costhetast);
//            bin1re->Fill(pth,costhetast,weight);
            bin1.at(cluster)->Fill(mhh,costhetast);
            bin1re.at(cluster)->Fill(mhh,costhetast,weight);    
             */
        //} else {cout<<"it is not gen level information"<<endl; return -1; }
    } else {cout<<"no two higsses gen level information, stoped in event: "<<endl; return -1; }
    //cout<<weight<<endl;
    return weight;
} // close gen level cuts
/////////////////////////////////////////////////////////////////////////
void draw() {
    //vector<TH2D> JHEP2D;
    TLegend *leg = new TLegend(0.6,0.7,0.95,0.92);
    leg->SetTextSize(0.04146853);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    TLegend *leg2 = new TLegend(0.2,0.7,0.4,0.92);
    leg2->SetTextSize(0.04146853);
    leg2->SetLineColor(1);
    leg2->SetLineStyle(1);
    leg2->SetLineWidth(1);
    leg2->SetFillColor(0);
/*    
    int color[16] = {
        221,//224,
        225,//228,
        2,//208,
        8,
    //
        221,//224,
        225,//228,
        2,//208,
        8,
        1};
    int line[16] = {
        1,1,1,1,//1,1,1,
        2,2,2,2,//2,2,2
        1};
 //*/
    ///*
    int color[19] = {
        8,//224,
        225,//228,
        2,
        //8,
        //
        8,//224,
        225,//228,
        2,
        //8,
        //
        8,//224,
        225,//228,
        2,
        //8,
        1};
    int line[19] = {
        1,1,1,//1,1,1,
        2,2,2,//2,2,2,
        4,4,4,//4,4,4,
        1};
    //*/
    // draw all together
    //TLatex Tl1; Tl1.SetTextFont(43); Tl1.SetTextSize(38); 
    TCanvas* PT_HAT = new TCanvas();
    PT_HAT->cd();
    PT_HAT->Clear();
    for(unsigned int cluster=0; cluster<10; cluster++) { // 10 to hh
        cout<<cluster<<endl;
        stringstream clus;
        clus << cluster;
        leg->SetHeader("m_{Q} (GeV)");
        leg2->SetHeader("#kappa = 1 (ggF)");
        //leg2->SetHeader("Quark initiated");
        //leg2->SetHeader("QQ");
        double normalize0 = 1./(genmbb[cluster]->Integral());
        //double normalize1 = 1./(histoJHEP_mhh[cluster][0]->Integral());
        cout<<"sumV1: "<< 1/normalize0<<" JHEP: "<<1/normalize0<<endl;
        /////////////////////////////////////////////////////////
        PT_HAT->Clear();
        genmbb[cluster]->SetTitle("");
        genmbb[cluster]->GetYaxis()->SetTitle("% of events / bin");
        genmbb[cluster]->SetLineColor(color[cluster]);//cluster+1); //color[cluster]
        genmbb[cluster]->SetLineWidth(3);
        genmbb[cluster]->SetLineStyle(line[cluster]);
        //histoJHEP_mhh[cluster][0]->SetLineColor(8);
        //histoJHEP_mhh[cluster][0]->SetLineWidth(3);
        //histoJHEP_mhh[cluster][0]->Scale(normalize1);
        genmbb[cluster]->Scale(normalize0);
        //////
        // to hh+ jets
        //////
        /*
        if(cluster==8) genmbb[cluster]->Scale(0.4);
        genmbb[cluster]->SetMaximum(0.28);//0.21);//0.61
        //genmbb[cluster][0]->SetMaximum(1.2*genmbb[cluster][0]->GetMaximum());
        if(cluster<4)leg->AddEntry(genmbb[cluster],JHEPbench[cluster],"l");
        //if(cluster==1)leg2->AddEntry(genmbb[cluster],"Top partner","l");
        //if(cluster==1)leg2->AddEntry(genmbb[cluster],"Bottom partner","l");
        if(cluster==1)leg2->AddEntry(genmbb[cluster],"Up partner","l");
        if(cluster==4)leg2->AddEntry(genmbb[cluster],"Down partner","l");
        if(cluster==8)leg2->AddEntry(genmbb[cluster],"SM","l");
        // */
        /////
        // to hh only
        ////
        ///*
        //if(cluster==9) genmbb[cluster]->Scale(0.25);
        genmbb[cluster]->SetMaximum(0.25);//0.21);//0.61
        //genmbb[cluster][0]->SetMaximum(1.2*genmbb[cluster][0]->GetMaximum());
        if(cluster<3)leg->AddEntry(genmbb[cluster],JHEPbench[cluster],"l");
        if(cluster==1)leg2->AddEntry(genmbb[cluster],"Top partner","l");
        //if(cluster==1)leg2->AddEntry(genmbb[cluster],"Bottom partner","l");
        if(cluster==4)leg2->AddEntry(genmbb[cluster],"Up partner","l");
        if(cluster==6)leg2->AddEntry(genmbb[cluster],"Down partner","l");
        if(cluster==9)leg2->AddEntry(genmbb[cluster],"SM","l");
        //*/
         ///
        normalize0 = 1./(genpth[cluster]->Integral());
        genpth[cluster]->SetTitle("");
        genpth[cluster]->GetYaxis()->SetTitle("% of events / bin");
        genpth[cluster]->SetLineColor(color[cluster]);//cluster+1);
        genpth[cluster]->SetLineWidth(3);
        genpth[cluster]->Scale(normalize0);
        genpth[cluster]->SetLineStyle(line[cluster]);
        genpth[cluster]->SetMaximum(0.53);
        //
        normalize0 = 1./(genpthsub[cluster]->Integral());
        genpthsub[cluster]->SetTitle("");
        genpthsub[cluster]->GetYaxis()->SetTitle("% of events / bin");
        genpthsub[cluster]->SetLineColor(color[cluster]);//cluster+1);
        genpthsub[cluster]->SetLineWidth(3);
        genpthsub[cluster]->Scale(normalize0);
        genpthsub[cluster]->SetLineStyle(line[cluster]);
        genpthsub[cluster]->SetMaximum(0.79);
        //
        normalize0 = 1./(gencost[cluster]->Integral());
        gencost[cluster]->SetTitle("");
        gencost[cluster]->GetYaxis()->SetTitle("% of events / bin");
        gencost[cluster]->SetLineColor(color[cluster]);//cluster+1);
        gencost[cluster]->SetLineWidth(3);
        gencost[cluster]->Scale(normalize0);
        gencost[cluster]->SetMinimum(0);
        gencost[cluster]->SetMaximum(0.5);
        gencost[cluster]->SetLineStyle(line[cluster]);
        //
        normalize0 = 1./(gencost2[cluster]->Integral());
        gencost2[cluster]->SetTitle("");
        gencost2[cluster]->GetYaxis()->SetTitle("% of events / bin");
        gencost2[cluster]->SetLineColor(color[cluster]);//cluster+1);
        gencost2[cluster]->SetLineWidth(3);
        gencost2[cluster]->Scale(normalize0);
        gencost2[cluster]->SetMinimum(0);
        gencost2[cluster]->SetMaximum(0.12);//0.5
        gencost2[cluster]->SetLineStyle(line[cluster]);
        //
        normalize0 = 1./(Njets_passing_kLooseID[cluster]->Integral());
        Njets_passing_kLooseID[cluster]->SetTitle("");
        Njets_passing_kLooseID[cluster]->GetYaxis()->SetTitle("% of events / bin");
        Njets_passing_kLooseID[cluster]->SetLineColor(color[cluster]);//cluster+1);
        Njets_passing_kLooseID[cluster]->SetLineWidth(3);
        Njets_passing_kLooseID[cluster]->Scale(normalize0);
        Njets_passing_kLooseID[cluster]->SetLineStyle(line[cluster]);
        //
        normalize0 = 1./(btagselected[cluster]->Integral());
        btagselected[cluster]->SetTitle("");
        btagselected[cluster]->GetYaxis()->SetTitle("% of events / bin");
        btagselected[cluster]->SetLineColor(color[cluster]);//cluster+1);
        btagselected[cluster]->SetLineWidth(3);
        btagselected[cluster]->Scale(normalize0);
        btagselected[cluster]->SetLineStyle(line[cluster]);
        //
        normalize0 = 1./(fullmhh[cluster]->Integral());
        fullmhh[cluster]->SetTitle("");
        fullmhh[cluster]->GetYaxis()->SetTitle("% of events / bin");
        fullmhh[cluster]->SetLineColor(color[cluster]);//cluster+1);
        fullmhh[cluster]->SetLineWidth(3);
        fullmhh[cluster]->Scale(normalize0);
        
        fullmhh[cluster]->SetMaximum(0.4);
        fullmhh[cluster]->SetLineStyle(line[cluster]);
        //
        normalize0 = 1./(gendeta[cluster]->Integral());
        gendeta[cluster]->SetTitle("");
        gendeta[cluster]->GetYaxis()->SetTitle("% of events / bin");
        gendeta[cluster]->SetLineColor(color[cluster]);//cluster+1);
        gendeta[cluster]->SetLineWidth(3);
        gendeta[cluster]->Scale(normalize0);
        gendeta[cluster]->SetMaximum(0.4);
        gendeta[cluster]->SetLineStyle(line[cluster]);
        //
        normalize0 = 1./(fullpth1[cluster]->Integral());
        fullpth1[cluster]->SetTitle("");
        fullpth1[cluster]->GetYaxis()->SetTitle("% of events / bin");
        fullpth1[cluster]->SetLineColor(color[cluster]);//cluster+1);
        fullpth1[cluster]->SetLineWidth(3);
        fullpth1[cluster]->Scale(normalize0);
        fullpth1[cluster]->SetMaximum(0.25);
        fullpth1[cluster]->SetLineStyle(line[cluster]);
        //
        normalize0 = 1./(leptop[cluster]->Integral());
        leptop[cluster]->SetTitle("");
        leptop[cluster]->GetYaxis()->SetTitle("% of events / bin");
        leptop[cluster]->SetLineColor(color[cluster]);//cluster+1);
        leptop[cluster]->SetLineWidth(3);
        leptop[cluster]->Scale(normalize0);
        leptop[cluster]->SetMaximum(0.5);
        leptop[cluster]->SetLineStyle(line[cluster]);
        //
        normalize0 = 1./(fullpth2[cluster]->Integral());
        fullpth2[cluster]->SetTitle("");
        fullpth2[cluster]->GetYaxis()->SetTitle("% of events / bin");
        fullpth2[cluster]->SetLineColor(color[cluster]);//cluster+1);
        fullpth2[cluster]->SetLineWidth(3);
        fullpth2[cluster]->Scale(normalize0);
        fullpth2[cluster]->SetMaximum(0.2);
        fullpth2[cluster]->SetLineStyle(line[cluster]);
        //
        normalize0 = 1./(hadtop[cluster]->Integral());
        hadtop[cluster]->SetTitle("");
        hadtop[cluster]->GetYaxis()->SetTitle("% of events / bin");
        hadtop[cluster]->SetLineColor(color[cluster]);//cluster+1);
        hadtop[cluster]->SetLineWidth(3);
        hadtop[cluster]->Scale(normalize0);
        hadtop[cluster]->SetMaximum(0.35);
        hadtop[cluster]->SetLineStyle(line[cluster]);
        //
        normalize0 = 1./(toponept[cluster]->Integral());
        toponept[cluster]->SetTitle("");
        toponept[cluster]->GetYaxis()->SetTitle("% of events / bin");
        toponept[cluster]->SetLineColor(color[cluster]);//cluster+1);
        toponept[cluster]->SetLineWidth(3);
        toponept[cluster]->Scale(normalize0);
        toponept[cluster]->SetMaximum(0.29);
        toponept[cluster]->SetLineStyle(line[cluster]);
        // fullpth1->Fill(genjet.at(0).pt);  leptop1-
        // [cluster]  fullmhh
        // 
        //if(cluster >0)genmbb[cluster]->Draw("SAME");
        //histoJHEP_mhh[cluster][0]->Draw("same");
        } // close to cluster
        
    TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(38); 
    //Tl.DrawText(.1, .5,   "#hat{#kappa} = 1");
    //
        PT_HAT->SetLogx(0);
        genmbb[0]->Draw();
        leg->Draw("same");
    leg2->Draw("same");
        for(unsigned int cluster=1; cluster<10; cluster++) genmbb[cluster]->Draw("SAME");
        string namefilebenchmhh = "plotKin/clu_mhh.pdf";
        const char * filemhh = namefilebenchmhh.c_str();
        PT_HAT->SaveAs(filemhh);//JHEP2Dmhhclu[cluster]); 
        PT_HAT->Clear();
        genpth[0]->Draw();
        leg->Draw("same");
        
    leg2->Draw("same");
        for(unsigned int cluster=1; cluster<10; cluster++) genpth[cluster]->Draw("SAME");
    //Tl.DrawText(200, .1, "#hat{#kappa} = 1");
        namefilebenchmhh = "plotKin/clu_pth.pdf";
        const char * filepth = namefilebenchmhh.c_str();
        PT_HAT->SaveAs(filepth);//JHEP2Dmhhclu[cluster]); 
        PT_HAT->Clear();
    genpthsub[0]->Draw();
    leg->Draw("same");
    leg2->Draw("same");
    for(unsigned int cluster=1; cluster<10; cluster++) genpthsub[cluster]->Draw("SAME");
    namefilebenchmhh = "plotKin/clu_pthsub.pdf";
    const char * filepthsub = namefilebenchmhh.c_str();
    PT_HAT->SaveAs(filepthsub);//JHEP2Dmhhclu[cluster]); 
    PT_HAT->Clear();
    //
    PT_HAT->SetLogx(0);
    gencost[0]->Draw();
    leg->Draw("same");
    leg2->Draw("same");
    for(unsigned int cluster=1; cluster<10; cluster++) gencost[cluster]->Draw("SAME");
    namefilebenchmhh = "plotKin/clu_cost.pdf";
    const char * filecost = namefilebenchmhh.c_str();
    PT_HAT->SaveAs(filecost);
    //
    PT_HAT->SetLogx(0);
    gencost2[0]->Draw();
    leg->Draw("same");
    leg2->Draw("same");
    for(unsigned int cluster=1; cluster<11; cluster++) gencost2[cluster]->Draw("SAME");
    namefilebenchmhh = "plotKin/clu_cost2.pdf";
    const char * filecost2 = namefilebenchmhh.c_str();
    PT_HAT->SaveAs(filecost2);
    //
    Njets_passing_kLooseID[0]->Draw();
    leg->Draw("same");
    leg2->Draw("same");
    for(unsigned int cluster=1; cluster<10; cluster++) Njets_passing_kLooseID[cluster]->Draw("SAME");
    namefilebenchmhh = "plotKin/clu_njets.pdf";
    const char * filenjets = namefilebenchmhh.c_str();
    PT_HAT->SaveAs(filenjets);
    //
    gendeta[0]->Draw();
    leg->Draw("same");
    leg2->Draw("same");
    for(unsigned int cluster=0; cluster<10; cluster++) gendeta[cluster]->Draw("SAME");
    namefilebenchmhh = "plotKin/clu_Deta.pdf";
    const char * filedr = namefilebenchmhh.c_str();
    PT_HAT->SaveAs(filedr);
    //
    fullpth1[0]->Draw();
    leg->Draw("same");
    leg2->Draw("same");
    for(unsigned int cluster=0; cluster<10; cluster++) fullpth1[cluster]->Draw("SAME");
    namefilebenchmhh = "plotKin/clu_j1pt.pdf";
    const char * filej1pt = namefilebenchmhh.c_str();
    PT_HAT->SaveAs(filej1pt);
    //
    leptop[0]->Draw();
    leg->Draw("same");
    leg2->Draw("same");
    for(unsigned int cluster=0; cluster<10; cluster++) leptop[cluster]->Draw("SAME");
    namefilebenchmhh = "plotKin/clu_j1eta.pdf";
    const char * filej1eta = namefilebenchmhh.c_str();
    PT_HAT->SaveAs(filej1eta);
    //
    fullpth2[0]->Draw();
    leg->Draw("same");
    leg2->Draw("same");
    for(unsigned int cluster=0; cluster<10; cluster++) fullpth2[cluster]->Draw("SAME");
    namefilebenchmhh = "plotKin/clu_j2pt.pdf";
    const char * filej2pt = namefilebenchmhh.c_str();
    PT_HAT->SaveAs(filej2pt);
    //
    hadtop[0]->Draw();
    leg->Draw("same");
    leg2->Draw("same");
    for(unsigned int cluster=0; cluster<10; cluster++) hadtop[cluster]->Draw("SAME");
    namefilebenchmhh = "plotKin/clu_j2eta.pdf";
    const char * filej2eta = namefilebenchmhh.c_str();
    PT_HAT->SaveAs(filej2eta);
    //
    toponept[0]->Draw();
    leg->Draw("same");
    leg2->Draw("same");
    for(unsigned int cluster=0; cluster<10; cluster++) toponept[cluster]->Draw("SAME");
    namefilebenchmhh = "plotKin/clu_toponePt.pdf";
    const char * filetopo = namefilebenchmhh.c_str();
    PT_HAT->SaveAs(filetopo);    
    
    
} // close draw
/////////////////////////////////////////////////////////////////////////
/*
void draw_out() {
    cout<<endl<<"Draw with outlayers "<<endl;
    // the first sample is the benchmark
    //vector<TH2D> JHEP2D;
    TLegend *leg = new TLegend(0.4,0.7,0.85,0.92);
    leg->SetTextSize(0.04146853);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    // draw all together
    TLatex Tl; Tl.SetTextFont(43); Tl.SetTextSize(38); 
    for(unsigned int cluster=0; cluster<8; cluster++) { 
        stringstream clus;
        clus << cluster;
        
        TCanvas* PT_HAT = new TCanvas();
        PT_HAT->cd();
        //
        leg->SetHeader(JHEPbench[cluster]);
        //////////////////////////////////////////////////
         double normalize0 = 1./(genmbb[cluster]->Integral());
        PT_HAT->Clear();
        genmbb[cluster]->SetTitle("");
        genmbb[cluster]->GetYaxis()->SetTitle("% of events / (20 GeV)");
        genmbb[cluster]->SetLineColor(1);
        genmbb[cluster]->SetLineWidth(3);
        genmbb[cluster]->Scale(normalize0);
        genmbb[cluster]->SetMaximum(2.5*genmbb[cluster][0]->GetMaximum());
        genmbb[cluster]->Draw();
        if(cluster ==0){
            leg->AddEntry(genmbb[cluster],"Benchmark","l");
        }
        leg->Draw("same");
        genmbb[cluster]->Draw("same");
        string namefilebenchmhh = "plotKin/clu" + clus.str() + "_mhh_out.pdf";
        const char * filemhh = namefilebenchmhh.c_str();
        PT_HAT->SaveAs(filemhh);
        ////////////////////////////
        normalize0 = 1./(fullmhh[cluster]->Integral());
        PT_HAT->Clear();
        fullmhh[cluster]->SetTitle("");
        fullmhh[cluster]->GetYaxis()->SetTitle("% of events / (20 GeV)");
        fullmhh[cluster]->SetLineColor(1);
        fullmhh[cluster]->SetLineWidth(3);
        fullmhh[cluster]->Scale(normalize0);
        fullmhh[cluster]->SetMaximum(2.8*fullmhh[cluster][0]->GetMaximum());
        fullmhh[cluster]->Draw();
        leg->Draw("same");
        fullmhh[cluster]->Draw("same");
        string namefilebenchmhhfull = "plotKin/clu" + clus.str() + "_mhh_out_full.pdf";
        const char * filemhhfull = namefilebenchmhhfull.c_str();
        PT_HAT->SaveAs(filemhhfull);
        ////////////////////////////
        ////////////////////////////
        // 
        PT_HAT->Clear();
    } // close to cluster
    /////////////////////////////////////////////////
    
} // close draw outlayer 
 */
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
int HH4b(vector<PseudoJet> jets, vector<int> btag, vector<int> btrue, double weight, int cluster, int outlayer){
    ///////////////////////////////////////////////////
    // to use b--tag
    vector<PseudoJet> btaggedjets; 
    bool dobtag = true;
    if (dobtag) {for(unsigned n=0; n<jets.size(); n++) if( btag.at(n)>0) btaggedjets.push_back(jets.at(n)); }
        //if(btrue[n] > 0) btaggedjets.push_back(jets.at(n)); } 
    else btaggedjets = jets; // not to use bt--tag
    ////////////////////////////////
    int njets=btaggedjets.size();
    //cout<<njets<<endl;
    if (njets>3){ // && nbtag >3
        //fullmhh.at(cluster)->Fill((btaggedjets.at(0) + btaggedjets.at(1)+btaggedjets.at(2) + btaggedjets.at(3)).m(),weight);
        fullmhh[cluster]->Fill((btaggedjets.at(0) + btaggedjets.at(1)+btaggedjets.at(2) + btaggedjets.at(3)).m(),weight);
    } 
    return njets;
    
}
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
int isbtagged(vector<PseudoJet> jets, vector<int> & btag, vector<int> & bmistag, vector<int> & btrue);
/////////////////////////////////////////////////////////////////
int recojets(vector<PseudoJet> particles,vector<PseudoJet> & jets_akt, vector<int> & btag, vector<int> & bmistag, vector<int> & fattag, vector<int> & btrue, double weight, int cluster, int outlayer){
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
    //cout <<weight<<endl;
    //Njets_passing_kLooseID.at(cluster)->Fill(njets,weight);
    //btagselected.at(cluster)->Fill(nbtag,weight); 
    Njets_passing_kLooseID[cluster]->Fill(njets,weight);
    btagselected[cluster]->Fill(nbtag,weight); 
    
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
                if(constitu.at(j).pt() > bjetpt && abs( constitu.at(j).eta()) < etab) {see++; nbtag++;}// to reco btag 
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
        if(seetruth>0 && see>0) btrue.push_back(id); else btrue.push_back(0); // to mctruth
        //cout<<"b-quarks mistagged = " <<bmistag[i] <<" b-quark = " <<btag[i] <<endl;
    } // close for each jet
    //int numbb=0; for(unsigned i=0 ; i< btrue.size() ; i++ ) numbb += btrue[i];
    //cout<<" total "<<numbb<<endl;
    return nbtag;
} // close isbtagged
/////////////////////////////////////////////////////////////////////////

// save the histos
int save_hist(int isample){
    cout<<"saved histos "<<Nlep_passing_kLooseID.size()<<endl;
    
    for(unsigned int cluster=0; cluster<19; cluster++){ 
        
        //genmbb[cluster][outlayer]=(TH1D*) genmbb1->Clone();
        //genmbb[cluster][outlayer]->SetDirectory(0);
        
        const char* Mass;
        Mass = Form("higgs_4b_shower_%d_%d.root",cluster,0); 
        cout<<" saved "<< Form("higgs_4b_shower_%d_%d.root",cluster,0)<<endl;
        TFile f1(Mass, "recreate");
        f1.cd();
        
        //Njets_passing_kLooseID[cluster][outlayer]->Write();
        //Nlep_passing_kLooseID[cluster][outlayer]->Write();
        //btagselected[cluster][outlayer]->Write();
        // save to plot
        /*
        REmhh[cluster][outlayer]=genmbb[cluster][outlayer]; 
        REpt[cluster][outlayer]=genpth[cluster][outlayer]; 
        REcost[cluster][outlayer]=gencost[cluster][outlayer]; 
        JHEP2D[cluster][outlayer]=histoJHEP_bin1[cluster][outlayer];
        JHEPmhh[cluster][outlayer]=histoJHEP_mhh[cluster][outlayer]; 
        JHEPpt[cluster][outlayer]=histoJHEP_pth[cluster][outlayer]; 
        JHEPcost[cluster][outlayer]=histoJHEP_cost[cluster][outlayer];
         */
        //
        genmbb[cluster]->Write();
        fullmhh[cluster]->Write();
        gencost[cluster]->Write();
        //bin1re[cluster]->Write(); 
        genpth[cluster]->Write();
        genpthsub[cluster]->Write();
        //
        //fullmhh[cluster][outlayer]->Write();
        //fullpth1[cluster][outlayer]->Write();
        //fullpth2[cluster][outlayer]->Write();
        //fullcost[cluster][outlayer]->Write();
         
        f1.Close();
    }
    /*
    for(unsigned int cluster=0; cluster<12; cluster++) { 
        const char* Mass;
        Mass = Form("higgs_4b_shower_%d.root",cluster); 
        cout<<" saved "<< Form("higgs_4b_shower_%d.root",cluster)<<endl;
        TFile f1(Mass, "recreate");
        f1.cd();
        
        Njets_passing_kLooseID.at(cluster)->Write();
    Nlep_passing_kLooseID.at(cluster)->Write();
    btagselected.at(cluster)->Write();
    // save to plot
    REmhh.push_back(genmbb.at(cluster)); 
    REpt.push_back(genpth.at(cluster)); 
    REcost.push_back(gencost.at(cluster)); 
    JHEP2D.push_back(histoJHEP_bin1.at(cluster));
    JHEPmhh.push_back(histoJHEP_mhh.at(cluster)); 
    JHEPpt.push_back(histoJHEP_pth.at(cluster)); 
    JHEPcost.push_back(histoJHEP_cost.at(cluster));
    //
    genmbb.at(cluster)->Write();
    gencost.at(cluster)->Write();
    bin1re.at(cluster)->Write(); 
    genpth.at(cluster)->Write();
    //
    fullmhh.at(cluster)->Write();
    fullpth1.at(cluster)->Write();
    fullpth2.at(cluster)->Write();
    fullcost.at(cluster)->Write();
    f1.Close();
    }
     */
    //
    
    //Njets_passing_kLooseID->Reset();
    //Nlep_passing_kLooseID->Reset();
    //genpth->Reset();
    //genmbb->Reset();
    //gencost->Reset();
    /////////    

    /////////
    // save to plot
    
    /////////////////////
    //TFile f2("Distros_V1_5p_20000ev_13sam_13TeV_all.root", "recreate");
    //f2.cd();
    //bin1->Write();    
    //f2.Close();
    //bin1->Reset();
    ////////////////////

    //  basicLeptons[0]->Reset();
    //  basicLeptons[1]->Reset();
    //  basicLeptons[2]->Reset();
    //  for(unsigned i=0;i<10;i++) basicHadtop[i]->Reset();
    //  for(unsigned i=0;i<13;i++) basicLeptop[i]->Reset();

    return 0;
}
/////////////////////////////////////////////////////////////////////////
int decla(int mass, int teste){
    
    /*
    //Distros_5p_20000ev_13sam_13TeV.root
    //string folder1_st = "0-851";// from makedistros
    TString filenameJHEP;    
    std::stringstream sstr;
    sstr << "Distros_5p_500000ev_12sam_13TeV_JHEP_500K.root";// "Distros_JHEP_5p_20000ev_12sam_13TeV.root"; // this come from makedistos
    filenameJHEP = sstr.str();
    TFile* fJHEP = new TFile(filenameJHEP);
    if (fJHEP == NULL) cout<<" something wrong with file "<<endl;
    cout<<" found the  file "<< filenameJHEP <<endl;
    //fJHEP->cd(folder1_st.c_str());    
    /////
    fJHEP->cd();
    for(unsigned int cluster=0; cluster<12; cluster++)  {
        
    string Result;          // string which will contain the result
    ostringstream convert;   // stream used for the conversion
    convert << cluster;      // insert the textual representation of 'Number' in the characters in the strea
    //Result = convert.str(); // set 'Result' to the contents of the stream
    /////
    TString fname = convert.str() +"_bin1;1";
    TString fnamemhh = convert.str() +"_mhh;1";
    TString fnamept = convert.str() +"_pt;1";
    TString fnamecost = convert.str() +"_hcths;1";
        
    //cout <<fname<<endl;
    //
    TH2D* histoJHEP_bin11 = (TH2D*)fJHEP->Get(fname);
    TH1D* histoJHEP_mhh1 = (TH1D*)fJHEP->Get(fnamemhh);
    TH1D* histoJHEP_pth1 = (TH1D*)fJHEP->Get(fnamept);
    TH1D* histoJHEP_cost1 = (TH1D*)fJHEP->Get(fnamecost);
        
    histoJHEP_bin1.push_back(histoJHEP_bin11);
    histoJHEP_mhh.push_back( histoJHEP_mhh1);
    histoJHEP_pth.push_back( histoJHEP_pth1);
    histoJHEP_cost.push_back( histoJHEP_cost1);

    histoJHEP_bin11->SetDirectory(0);
    histoJHEP_mhh1->SetDirectory(0);
    histoJHEP_pth1->SetDirectory(0);
    histoJHEP_cost1->SetDirectory(0);
        
    }//close cluster
    fJHEP->Close();
     */

    //////////////////////////////////////////////////////////////////////////
    // with the outlayers - 6 per cluster - 72 histograms
    //Distros_envelope_5p_20000ev_6sam_13TeV.root
    
    basicGen.resize(30);
    basicLeptons.resize(30);
    basicHadtop.resize(30);
    basicLeptop.resize(30);
    //////////
    Njets_passing_kLooseID.resize(30);
    Nlep_passing_kLooseID.resize(30);
    btagselected.resize(30);
    leptop.resize(30);
    hadtop.resize(30);
    fullmhh.resize(30);
    fullpth1.resize(30);
    fullpth2.resize(30);
    fullcost.resize(30);
    
    genmbb.resize(30); 
    genpth.resize(30);
    genpthsub.resize(30);
    gencost.resize(30);
    gencost2.resize(30);
    gendeta.resize(30);
    
    histoJHEP_mhh.resize(30);
    histoJHEP_pth.resize(30);
    histoJHEP_cost.resize(30);
    
    
    JHEPmhh.resize(30);
    JHEPcost.resize(30);
    JHEPpt.resize(30);
    
    REmhh.resize(30);
    REcost.resize(30);
    REpt.resize(30);
    toponept.resize(30);
    
    for(unsigned int cluster=0; cluster<19; cluster++) {
        //histoOutlayers->push_back();
        //stringstream clus;
        //clus << cluster+1;
        //string folder1_st = "clu"+clus.str()+";1";
        //fEnv->cd(folder1_st.c_str());
         
        

      
         
    
    const char* label="";
    
    TH1D *Njets_passing_kLooseID1 = new TH1D("njets_passing_kLooseID_ct4",  
                                      label, 
                                      17, 0.5, 18.5);
    Njets_passing_kLooseID1->GetYaxis()->SetTitle("");
    Njets_passing_kLooseID1->GetXaxis()->SetTitle("Njets after showering"); 
    Njets_passing_kLooseID[cluster]=(TH1D*) Njets_passing_kLooseID1->Clone();
    Njets_passing_kLooseID[cluster]->SetDirectory(0);

    TH1D *Nlep_passing_kLooseID1 = new TH1D("nlep_passing_kLooseID_ct4",  
                                      label, 
                                      13, -0.5, 12.5);
    Nlep_passing_kLooseID1->GetYaxis()->SetTitle("");
    Nlep_passing_kLooseID1->GetXaxis()->SetTitle("Nleptons after showering"); 
    Nlep_passing_kLooseID[cluster]=(TH1D*) Nlep_passing_kLooseID1->Clone();
    Nlep_passing_kLooseID[cluster]->SetDirectory(0);
        
    TH1D *btagselected1 = new TH1D("btagselected",  
                            label, 
                            13, -0.5, 12.5);
    btagselected1->GetYaxis()->SetTitle("");
    btagselected1->GetXaxis()->SetTitle("b-tagable b's on selected events");
    btagselected[cluster]=(TH1D*) btagselected1->Clone();
    btagselected[cluster]->SetDirectory(0);
    
    TH1D *genmbb1 = new TH1D("higgs_mhh",  
                      label, 
                      30,110.,1000.);
                            // 30,380.,3000.);
                            //  30,100.,3000.);
    genmbb1->GetYaxis()->SetTitle("");
    genmbb1->GetXaxis()->SetTitle("m_{hh}^{gen}");
            genmbb[cluster]=(TH1D*) genmbb1->Clone();
            genmbb[cluster]->SetDirectory(0);
            //genmbb1->SetDirectory(0);

        TH1D *toponept1 = new TH1D("toponept",  
                           label, 
                           30,0.,2000.);
        toponept1->GetYaxis()->SetTitle("");
        toponept1->GetXaxis()->SetTitle("p_{T}^{Q} (GeV)");
        toponept[cluster]=(TH1D*) toponept1->Clone();
        toponept[cluster]->SetDirectory(0);
        //genmbb1->SetDirectory(0);
        


    TH1D *genpth1 = new TH1D("higgs_pt",  
                      label, 
                      30,0.,2000.);
    genpth1->GetYaxis()->SetTitle("");
    genpth1->GetXaxis()->SetTitle("leading p_{T}^{h, gen}");     
        genpth[cluster]=(TH1D*) genpth1->Clone();
        genpth[cluster]->SetDirectory(0);

        TH1D *genpth1sub = new TH1D("higgs_ptsub",  
                                 label, 
                                 20,0.,1800.);
        genpth1sub->GetYaxis()->SetTitle("");
        genpth1sub->GetXaxis()->SetTitle("Sub-leading p_{T}^{h, gen}");     
        genpthsub[cluster]=(TH1D*) genpth1sub->Clone();
        genpthsub[cluster]->SetDirectory(0);
        
    TH1D *gencost1 = new TH1D("cost_h",  
                      label, 
                      //5,0.,1.);
                              20, 4, -4);
    gencost1->GetYaxis()->SetTitle("");
    //gencost1->GetXaxis()->SetTitle("higgs cos#theta*"); 
        gencost1->GetXaxis()->SetTitle("#eta_{h^{leading}}"); 
        gencost[cluster]=(TH1D*) gencost1->Clone();
        gencost[cluster]->SetDirectory(0);

        TH1D *gencost22 = new TH1D("cost2_h",  
                                  label, 
                                  15,0.,1.);
                                  //30, 4, -4);
        gencost22->GetYaxis()->SetTitle("");
        gencost22->GetXaxis()->SetTitle("higgs cos#theta^{*}"); 
        //gencost22->GetXaxis()->SetTitle("#eta_{h^{sub-leading}}"); 
        gencost2[cluster]=(TH1D*) gencost22->Clone();
        gencost2[cluster]->SetDirectory(0);

        TH1D *gendeta1 = new TH1D("deta_h",  
                                  label, 
                                  12,0.,6.);
        gendeta1->GetYaxis()->SetTitle("");
        gendeta1->GetXaxis()->SetTitle("higgs #Delta #eta_{hh}"); 
        gendeta[cluster]=(TH1D*) gendeta1->Clone();
        gendeta[cluster]->SetDirectory(0);
        
        /*
    TH2D *bin11 = new TH2D("all",  
                      label, 
                      90,0.,1800.,10,-1,1.);
    bin11->GetYaxis()->SetTitle("");
    bin11->GetXaxis()->SetTitle("mhh X cost*");
    bin1[cluster]=(TH2D*) bin11->Clone();
    bin1[cluster]->SetDirectory(0);
        
    TH2D *bin1re1 = new TH2D("allre",  
                    label, 
                    90,0.,1800.,10,-1,1.);
    bin1re1->GetYaxis()->SetTitle("");
    bin1re1->GetXaxis()->SetTitle("mhh X cost*"); 
    bin1re[cluster]=(TH2D*) bin1re1->Clone();
    bin1re[cluster]->SetDirectory(0);
        */
        
    TH1D *leptop1 = new TH1D("j1_pt",  
                      label, 
                      20, 4, -4);
    //leptop->SetLogY(1);    
    leptop1->GetYaxis()->SetTitle("");
    leptop1->GetXaxis()->SetTitle("Gen leading jet #eta"); 
        leptop[cluster]=(TH1D*) leptop1->Clone();
          leptop1->SetDirectory(0);
        
    TH1D *hadtop1 = new TH1D("j1_eta",  
                      label, 
                      20, 4, -4);
    hadtop1->GetYaxis()->SetTitle("");
    hadtop1->GetXaxis()->SetTitle("Gen sub-leading jet #eta"); 
        hadtop[cluster]=(TH1D*) hadtop1->Clone();
          hadtop1->SetDirectory(0);
        
    //////
    
    TH1D *fullmhh1 = new TH1D("higgs_mhh_full",  
                      label, 
                      60,0.,1800.);
    fullmhh1->GetYaxis()->SetTitle("");
    fullmhh1->GetXaxis()->SetTitle("di-higgs invariant mass (shower)"); 
    fullmhh[cluster]=(TH1D*) fullmhh1->Clone();
    fullmhh[cluster]->SetDirectory(0);

    TH1D *fullpth11 = new TH1D("higgs_pt1",  
                      label, 
                      50,0.,2000.);
    fullpth11->GetYaxis()->SetTitle("");
    fullpth11->GetXaxis()->SetTitle("Gen Leading jet p_{T}"); 
        fullpth1[cluster]=(TH1D*) fullpth11->Clone();
          fullpth11->SetDirectory(0);
        
    TH1D *fullpth21 = new TH1D("higgs_pt2",  
                        label, 
                        60,0.,2000.);
    fullpth21->GetYaxis()->SetTitle("");
    fullpth21->GetXaxis()->SetTitle("Gen sub-leading jet p_{T}"); 
        fullpth2[cluster]=(TH1D*) fullpth21->Clone();
          fullpth21->SetDirectory(0);
        
    TH1D *fullcost1 = new TH1D("cost_h_full",  
                       label, 
                       5,0.,1.);
    fullcost1->GetYaxis()->SetTitle("");
    fullcost1->GetXaxis()->SetTitle("higgs cos#theta^{*} (shower)"); 
        fullcost[cluster]=(TH1D*) fullcost1->Clone();    
          fullcost1->SetDirectory(0);
    ///////////////////////////////////////////////////////////////////////////////////
    }// close for outlayer

    /////////////////////////////////////////////////////////
    return 0;
}

void style (){
    TStyle *defaultStyle = new TStyle("defaultStyle","Default Style");
    defaultStyle->SetOptStat(0000);
    defaultStyle->SetOptFit(000); 
    defaultStyle->SetPalette(1);
    /////// pad ////////////
    defaultStyle->SetPadBorderMode(1);
    defaultStyle->SetPadBorderSize(1);
    defaultStyle->SetPadColor(0);
    defaultStyle->SetPadTopMargin(0.05);
    defaultStyle->SetPadBottomMargin(0.15);
    defaultStyle->SetPadLeftMargin(0.15);
    defaultStyle->SetPadRightMargin(0.16);
    /////// canvas /////////
    defaultStyle->SetCanvasBorderMode(0);
    defaultStyle->SetCanvasColor(0);
    defaultStyle->SetCanvasDefH(600);
    defaultStyle->SetCanvasDefW(600);
    /////// frame //////////
    defaultStyle->SetFrameBorderMode(0);
    defaultStyle->SetFrameBorderSize(1);
    defaultStyle->SetFrameFillColor(0); 
    defaultStyle->SetFrameLineColor(1);
    /////// label //////////
    defaultStyle->SetLabelOffset(0.005,"XY");
    defaultStyle->SetLabelSize(0.06,"XY");
    defaultStyle->SetLabelFont(46,"XY");
    /////// title //////////
    defaultStyle->SetTitleOffset(1.1,"X");
    defaultStyle->SetTitleSize(0.01,"X");
    defaultStyle->SetTitleOffset(1.25,"Y");
    defaultStyle->SetTitleSize(0.06,"Y");
    defaultStyle->SetTitleFont(44, "XYZ");
    /////// various ////////
    defaultStyle->SetNdivisions(505,"Y");
    defaultStyle->SetLegendBorderSize(0);  // For the axis titles:
    
    defaultStyle->SetTitleColor(1, "XYZ");
    defaultStyle->SetTitleFont(42, "XYZ");
    defaultStyle->SetTitleSize(0.055, "XYZ");
    
    // defaultStyle->SetTitleYSize(Float_t size = 0.02);
    defaultStyle->SetTitleXOffset(0.85);
    defaultStyle->SetTitleYOffset(1.05);
    // defaultStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset
    
    // For the axis labels:
    defaultStyle->SetLabelColor(1, "XYZ");
    defaultStyle->SetLabelFont(42, "XYZ");
    defaultStyle->SetLabelOffset(0.007, "XYZ");
    defaultStyle->SetLabelSize(0.042, "XYZ");
    
    // For the axis:
    defaultStyle->SetAxisColor(1, "XYZ");
    defaultStyle->SetStripDecimals(kTRUE);
    defaultStyle->SetTickLength(0.03, "XYZ");
    defaultStyle->SetNdivisions(7, "XYZ");
    defaultStyle->SetPadTickX(1);   // To get tick marks on the opposite side of the frame
    defaultStyle->SetPadTickY(1);
    defaultStyle->cd();
    return;
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
