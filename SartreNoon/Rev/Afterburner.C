
#ifdef __CLING__
#include <iostream>
#include <fstream>
#include <string>

#include "TROOT.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "NeutronGenerator.cxx+g"

#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <cmath>
#include "TMath.h"
#include "TDatabasePDG.h"

#include "NeutronGenerator.h"
#endif

void DrawXSection(TH1D* XS1, TH1D* XS2, TH1D* XS3, TH1D* XST);

void DrawEnergy(TH1D* E1, TH1D* E2, TH1D* E3);

void DrawPseudo(TH1D* Eta1,TH1D* Eta2,TH1D* Eta3);

void DrawNY(TH1D* NeutronY2, TH1D* NeutronY3);

void DrawNEta(TH1D* NeuttonEta2, TH1D* NeutronEta3);

void DrawNE(TH1D* NeutronE2, TH1D* NeutronE3);



void Afterburner()
{

    NeutronGenerator *gen = new NeutronGenerator();

    gen->SetStoreQA();
    gen->SetStoreGeneratorFunctions();
    gen->SetHadronicInteractionModel(NeutronGenerator::kGlauber);
    gen->Initialize();
    gen->SetRunMode(NeutronGenerator::kInterface);
    gen->ReadENDF(kTRUE);
    gen->LoadENDF("hENDF.root");
    gen->Setup();


    TFile *file = new TFile("example_Pb.root","READ");

    TTree *tree = dynamic_cast<TTree*>(file->Get("tree"));

    Long64_t nE = tree->GetEntries();

    TLorentzVector* lv = new TLorentzVector;

    TBranch *bc = tree->GetBranch("vm");

    bc->SetAddress(&lv);

    std::vector<Double_t> onon; //srotes rapiity
    std::vector<Double_t> onXn; // stores rapidity
    std::vector<Double_t> XnXn; // stores rapidity

    std::vector<Double_t> ononEta;//stores pseudo rapidity
    std::vector<Double_t> onXnEta;//stores pseudo rapidity
    std::vector<Double_t> XnXnEta;//stores pseudo rapidity

    std::vector<Int_t> on; //stores index
    std::vector<Int_t> ln; // stores index
    std::vector<Int_t> Xn; //stores index

    std::vector<Double_t> ononE;
    std::vector<Double_t> onXnE;
    std::vector<Double_t> XnXnE;

    std::vector<Double_t> onNeutronE;
    std::vector<Double_t> onNeutronEta;
    std::vector<Double_t> onNeutronY;
    
    std::vector<Double_t> lnNeutronE;
    std::vector<Double_t> lnNeutronEta;
    std::vector<Double_t> lnNeutronY;

    std::vector<Double_t> XnNeutronE;
    std::vector<Double_t> XnNeutronEta;
    std::vector<Double_t> XnNeutronY;



    Double_t Mv = 3.09;

    for(Long64_t index = 0; index < nE; index++)
    {
        bc->GetEntry(index);

        Double_t y = lv->Rapidity();

        Double_t eta= lv->Eta();

        Double_t k = 0.5*Mv*TMath::Exp(TMath::Abs(y));

        std::vector<Int_t> nNumbers = gen->runSartreNoon(k);

        if(nNumbers[0]==0 & nNumbers[1] == 0) 
        {
            onon.push_back(y);
            ononEta.push_back(eta);
            on.push_back(index);
            ononE.push_back(lv->Energy());
            gen->createSartreNeutrons(nNumbers[0],nNumbers[1],onNeutronE,onNeutronEta,onNeutronY);
        }
        else if(nNumbers[0]==0 || nNumbers[1]==0 )
        {
            onXn.push_back(y);
            onXnEta.push_back(eta);
            ln.push_back(index);
            onXnE.push_back(lv->Energy());
            gen->createSartreNeutrons(nNumbers[0],nNumbers[1],lnNeutronE,lnNeutronEta,lnNeutronY);
        }
        else if(nNumbers[0]!=0 & nNumbers[1]!=0) 
        {
            XnXn.push_back(y);
            XnXnEta.push_back(eta);
            Xn.push_back(index);
            XnXnE.push_back(lv->Energy());
            gen->createSartreNeutrons(nNumbers[0],nNumbers[1],XnNeutronE,XnNeutronEta,XnNeutronY);
        }
        else std::cout<<"Error in the Event"<<std::endl;
    }
    delete lv;
    delete tree;
    file->Close();
    delete file;

    std::cout<<"0n0n number of events: "<<onon.size()<<std::endl;
    std::cout<<"0nXn number of events: "<<onXn.size()<<std::endl;
    std::cout<<"XnXn number of events: "<<XnXn.size()<<std::endl;


    Double_t minRange = -6;
    Double_t maxRange = 6;
    // Create histograms for each rapidity list
    TH1D *hist1 = new TH1D("hist1", "Rapidity Distribution 1", 1000, minRange,maxRange);
    TH1D *hist2 = new TH1D("hist2", "Rapidity Distribution 2", 1000, minRange, maxRange);
    TH1D *hist3 = new TH1D("hist3", "Rapidity Distribution 3", 1000, minRange,maxRange);

    TH1D *Eta1 = new TH1D("Eta1", "PseudoRapidity Distribution 1", 1000,  -15,15);
    TH1D *Eta2 = new TH1D("Eta2", "PseudoRapidity Distribution 2", 1000, -15,15);
    TH1D *Eta3 = new TH1D("Eta3", "PseudoRapidity Distribution 3", 1000, -15,15);

    TH1D *E1 = new TH1D("E1", "Energy 1", 1000, 0 ,15);
    TH1D *E2 = new TH1D("E2", "Energy 2", 1000, 0, 15);
    TH1D *E3 = new TH1D("E3", "Energy 3", 1000, 0, 15);

    for (Int_t i; i<onon.size();i++) {
        
        // bc->GetEntry(on[i]);

        hist1->Fill(onon[i]);
        Eta1->Fill(ononEta[i]);
        E1->Fill(ononE[i]);

    }
    for (Int_t i; i<onXn.size();i++) {
        
        // bc->GetEntry(ln[i]);
        
        hist2->Fill(onXn[i]);
        Eta2->Fill(onXnEta[i]);
        E2->Fill(onXnE[i]);
    }
    
    for (Int_t i; i<XnXn.size();i++) {

        // bc->GetEntry(Xn[i]);

        hist3->Fill(XnXn[i]);
        Eta3->Fill(XnXnEta[i]);
        E3->Fill(XnXnE[i]);
    }

    // Create histograms for cross-section distribution
    TH1D *XS1 = new TH1D("XS1", "Xsection Distribution 1", 1000, minRange, maxRange);
    TH1D *XS2 = new TH1D("XS2", "Xsection Distribution 2", 1000, minRange,maxRange);
    TH1D *XS3 = new TH1D("XS3", "Xsection Distribution 3", 1000, minRange,maxRange);
    TH1D *XST = new TH1D("XST", "Total Cross section", 1000, minRange,maxRange);

    for(Int_t entry = 0; entry<1000; entry++)
    {
        XS1->SetBinContent(entry,(hist1->GetBinContent(entry)) * 524 * 1000 /100000);
        XS2->SetBinContent(entry,(hist2->GetBinContent(entry)) * 524 * 1000 /100000);
        XS3->SetBinContent(entry,(hist3->GetBinContent(entry)) * 524 * 1000 /100000);
        XST->SetBinContent(entry,XS1->GetBinContent(entry)+XS2->GetBinContent(entry)+XS3->GetBinContent(entry));
    }
    
    DrawXSection(XS1,XS2,XS3,XST);
    DrawEnergy(E1,E2,E3);
    DrawPseudo(Eta1,Eta2,Eta3);



    /**
     * 
     * Neutron Data in the following histograms
     * 
    */

//    TH1D *NeutronE1 = new TH1D("NeutronE1","Neutron energy",1000,0,5);
   TH1D *NeutronE2 = new TH1D("NeutronE2","Neutron energy",10000,0,5000);
   TH1D *NeutronE3 = new TH1D("NeutronE3","Neutron energy",10000,0,5000);
   
//    TH1D *NeutronEta1 = new TH1D("NeutronEta1", "Neutron Pseudo Rapidity",1000,-10,10);
   TH1D *NeutronEta2 = new TH1D("NeutronEta2", "Neutron Pseudo Rapidity",1000,-25,25);
   TH1D *NeutronEta3 = new TH1D("NeutronEta3", "Neutron Pseudo Rapidity",1000,-25,25);

//    TH1D *NeutronY1 = new TH1D("NeutronY1","Neutron Rapidity",1000,-6,6);
   TH1D *NeutronY2 = new TH1D("NeutronY2","Neutron Rapidity",1000,-25,25);
   TH1D *NeutronY3 = new TH1D("NeutronY3","Neutron Rapidity",1000,-25,25);


   for(Int_t i=0; i<lnNeutronE.size(); i++)
   {
    NeutronE2->Fill(lnNeutronE[i]);
    NeutronEta2->Fill(lnNeutronEta[i]);
    NeutronY2->Fill(lnNeutronEta[i]);
   }

   for(Int_t i=0; i<XnNeutronE.size(); i++)
   {
    NeutronE3->Fill(XnNeutronE[i]);
    NeutronEta3->Fill(XnNeutronEta[i]);
    NeutronY3->Fill(XnNeutronEta[i]);
   }

   std::cout<<"Number of Neutrons in onon: "<<onNeutronE.size()<<std::endl;
   std::cout<<"Number of Neutrons in onxn: "<<lnNeutronE.size()<<std::endl;
   std::cout<<"Number of Neutrons in xnxn: "<<XnNeutronE.size()<<std::endl;


   DrawNY(NeutronY2,NeutronY3);
   DrawNEta(NeutronEta2,NeutronEta3);
   DrawNE(NeutronE2,NeutronE3);
   

}

void DrawNY( TH1D* NeutronY2, TH1D* NeutronY3)
{
    TCanvas *c4 = new TCanvas("c4", "Cross-section Plot", 1200, 800);

    // NeutronY1->SetLineColor(kRed);
    NeutronY2->SetLineColor(kBlue);
    NeutronY3->SetLineColor(kGreen);

    // NeutronY1->GetXaxis()->SetTitle("Rapidity");
    // NeutronY1->GetYaxis()->SetTitle("Counts");
    // NeutronY1->Draw("L");
    NeutronY2->GetXaxis()->SetTitle("Rapidity");
    NeutronY2->GetYaxis()->SetTitle("Counts");
    NeutronY2->Draw("L");
    NeutronY3->GetXaxis()->SetTitle("Rapidity");
    NeutronY3->GetYaxis()->SetTitle("Counts");
    NeutronY3->Draw("SAME L");

    c4->SaveAs("Neutron_Rapidity_Boosted.png");
}

void DrawNEta( TH1D* NeutronEta2, TH1D* NeutronEta3)
{
    TCanvas *c5 = new TCanvas("c5", "Cross-section Plot", 1200, 800);

    // NeutronEta1->SetLineColor(kRed);
    NeutronEta2->SetLineColor(kBlue);
    NeutronEta3->SetLineColor(kGreen);

    // NeutronEta1->GetXaxis()->SetTitle("Pseudo Rapidity");
    // NeutronEta1->GetYaxis()->SetTitle("Counts");
    // NeutronEta1->Draw("L");
    NeutronEta2->GetXaxis()->SetTitle("Pseudo Rapidity");
    NeutronEta2->GetYaxis()->SetTitle("Counts");
    NeutronEta2->Draw("L");
    NeutronEta3->GetXaxis()->SetTitle("Pseudo Rapidity");
    NeutronEta3->GetYaxis()->SetTitle("Counts");
    NeutronEta3->Draw("SAME L");

    c5->SaveAs("Neutron_Eta_Boosted.png");
}


void DrawNE( TH1D* NeutronE2, TH1D* NeutronE3)
{
    TCanvas *c6 = new TCanvas("c6", "Cross-section Plot", 1200, 800);

    // NeutronE1->SetLineColor(kRed);
    NeutronE2->SetLineColor(kBlue);
    NeutronE3->SetLineColor(kGreen);

    // NeutronE1->GetXaxis()->SetTitle("Energy");
    // NeutronE1->GetYaxis()->SetTitle("Counts");
    // NeutronE1->Draw("L");
    // NeutronE2->GetXaxis()->SetRangeUser(0,15);
    NeutronE2->GetXaxis()->SetTitle("Energy");
    NeutronE2->GetYaxis()->SetTitle("Counts");
    NeutronE2->Draw("L");
    // NeutronE3->GetXaxis()->SetRangeUser(0,15);
    NeutronE3->GetXaxis()->SetTitle("Energy");
    NeutronE3->GetYaxis()->SetTitle("Counts");
    NeutronE3->Draw("SAME L");

    c6->SaveAs("Neutron_Energy_Boosted.png");
}





void DrawXSection(TH1D* XS1, TH1D* XS2, TH1D* XS3, TH1D* XST)
{// Create a canvas to draw histograms
    TCanvas *c1 = new TCanvas("c1", "Cross-section Plot", 1200, 800);

    // Set histogram styles
    XS1->SetLineColor(kRed);
    XS2->SetLineColor(kBlue);
    XS3->SetLineColor(kGreen);
    XST->SetLineColor(kBlack);

    // Draw histograms
    Int_t minRapidity = -6;
    Int_t maxRapidity = 6;
    Int_t maxY = 2000;
    XS1->GetXaxis()->SetTitle("Rapidity");
    XS1->GetYaxis()->SetTitle("Cross section in nb");
    XS1->GetXaxis()->SetRangeUser(minRapidity,maxRapidity);
    XS1->GetYaxis()->SetRangeUser(0,maxY);
    XS1->Draw("L");
    XS3->GetXaxis()->SetTitle("Rapidity");
    XS3->GetYaxis()->SetTitle("Cross section in nb");
    XS2->GetXaxis()->SetRangeUser(minRapidity,maxRapidity);
    XS2->GetYaxis()->SetRangeUser(0,maxY);
    XS2->Draw("SAME L");
    XS3->GetXaxis()->SetTitle("Rapidity");
    XS3->GetYaxis()->SetTitle("Cross section in nb");
    XS3->GetXaxis()->SetRangeUser(minRapidity,maxRapidity);
    XS3->GetYaxis()->SetRangeUser(0,maxY);
    XS3->Draw("SAME L");
    XST->GetXaxis()->SetTitle("Rapidity");
    XST->GetYaxis()->SetTitle("Cross section in nb");
    XST->GetXaxis()->SetRangeUser(minRapidity,maxRapidity);
    XST->GetYaxis()->SetRangeUser(0,maxY);
    XST->Draw("Same L");

    // Add legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(XS1, "0n0n", "l");
    legend->AddEntry(XS2, "0nXn", "l");
    legend->AddEntry(XS3, "XnXn", "l");
    legend->AddEntry(XST, "Total", "l");
    legend->Draw();

    // Draw canvas
    c1->Draw();
    c1->SaveAs("Classes_of_cross_section_NB.png");}

void DrawEnergy(TH1D* E1, TH1D* E2, TH1D* E3)
{
    TCanvas *c3 = new TCanvas("c3", "Energy", 1200, 800);

    E1->SetLineColor(kRed);
    E2->SetLineColor(kBlue);
    E3->SetLineColor(kGreen);


    E1->GetXaxis()->SetTitle("Energy");
    E1->GetYaxis()->SetTitle("Counts");
    E1->GetYaxis()->SetRangeUser(0, 500);
    E1->Draw("L");
    E2->GetXaxis()->SetTitle("Energy");
    E2->GetYaxis()->SetTitle("Counts");
    E2->GetYaxis()->SetRangeUser(0, 500);
    E2->Draw("SAME L");
    E3->GetXaxis()->SetTitle("Energy");
    E3->GetYaxis()->SetTitle("Counts");
    E3->GetYaxis()->SetRangeUser(0, 500);
    E3->Draw("SAME L");
    TLegend *legend3 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend3->AddEntry(E1, "0n0n", "l");
    legend3->AddEntry(E1, "0nXn", "l");
    legend3->AddEntry(E3, "XnXn", "l");
    // legend->AddEntry(EtaT, "Total", "l");
    legend3->Draw();

    c3->Draw();
    c3->SaveAs("Energy_NB.png");

}

void DrawPseudo(TH1D* Eta1,TH1D* Eta2,TH1D* Eta3)
{
    TCanvas *c2 = new TCanvas("c2", "Pseudo Rapidity", 1200, 800);

    
    Eta1->SetLineColor(kRed);
    Eta2->SetLineColor(kBlue);
    Eta3->SetLineColor(kGreen);

    // Draw histograms
    Eta1->GetXaxis()->SetTitle("Pseudo Rapidity");
    Eta1->GetYaxis()->SetTitle("Counts");
    Eta1->GetYaxis()->SetRangeUser(0, 500);
    Eta1->Draw("L");
    Eta2->GetXaxis()->SetTitle("Pseudo Rapidity");
    Eta2->GetYaxis()->SetTitle("Counts");
    Eta2->GetYaxis()->SetRangeUser(0, 500);
    Eta2->Draw("SAME L");
    Eta3->GetXaxis()->SetTitle("Pseudo Rapidity");
    Eta3->GetYaxis()->SetTitle("Counts");
    Eta3->GetYaxis()->SetRangeUser(0, 500);
    Eta3->Draw("SAME L");
    TLegend *legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend2->AddEntry(Eta1, "0n0n", "l");
    legend2->AddEntry(Eta1, "0nXn", "l");
    legend2->AddEntry(Eta3, "XnXn", "l");
    // legend->AddEntry(EtaT, "Total", "l");
    legend2->Draw();

    c2->Draw();
    c2->SaveAs("Pseudo_rapidity_NB.png");
}