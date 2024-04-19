
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

    TLorentzVector* lv = new TLorentzVector();

    TBranch *bc = tree->GetBranch("vm");

    bc->SetAddress(&lv);

    std::vector<Double_t> onon;
    std::vector<Double_t> onXn;
    std::vector<Double_t> XnXn;

    Double_t Mv = 3.09;

    for(Long64_t index = 0; index < nE; index++)
    {
        bc->GetEntry(index);

        Double_t y = lv->Rapidity();

        Double_t k = 0.5*Mv*TMath::Exp(TMath::Abs(y));

        std::vector<Int_t> nNumbers = gen->runSartreNoon(k);

        if(nNumbers[0]==0 & nNumbers[1] == 0) onon.push_back(y);
        else if(nNumbers[0]==0 || nNumbers[1]==0 )onXn.push_back(y);
        else if(nNumbers[0]!=0 & nNumbers[1]!=0) XnXn.push_back(y);
        else std::cout<<"Error in the Event"<<std::endl;
    }

    delete lv;
    delete tree;
    file->Close();
    delete file;

    std::cout<<"0n0n number of events: "<<onon.size()<<std::endl;
    std::cout<<"0nXn number of events: "<<onXn.size()<<std::endl;
    std::cout<<"XnXn number of events: "<<XnXn.size()<<std::endl;

    // Create histograms for each rapidity list
    TH1F *hist1 = new TH1F("hist1", "Rapidity Distribution 1", 100, -5.0, 5.0);
    TH1F *hist2 = new TH1F("hist2", "Rapidity Distribution 2", 100, -5.0, 5.0);
    TH1F *hist3 = new TH1F("hist3", "Rapidity Distribution 3", 100, -5.0, 5.0);

    // Fill histograms with rapidity values
    for (double value : onon) {
        hist1->Fill(value);
    }
    for (double value : onXn) {
        hist2->Fill(value);
    }
    for (double value : XnXn) {
        hist3->Fill(value);
    }

    // Create histograms for cross-section distribution
    TH1F *XS1 = new TH1F("XS1", "Xsection Distribution 1", 100, -3.5, 4.0);
    TH1F *XS2 = new TH1F("XS2", "Xsection Distribution 2", 100, -3.5, 4.0);
    TH1F *XS3 = new TH1F("XS3", "Xsection Distribution 3", 100, -3.5, 4.0);

    // Fill cross-section histograms
    for (double value : onon) {
        Int_t bin = hist1->FindBin(value);
        Double_t crossSection = hist1->GetBinContent(bin) * 524/100000;
        XS1->Fill(value, crossSection);
    }
    for (double value : onXn) {
        Int_t bin = hist2->FindBin(value);
        Double_t crossSection = hist2->GetBinContent(bin) * 524/100000;
        XS2->Fill(value, crossSection);
    }
    for (double value : XnXn) {
        Int_t bin = hist3->FindBin(value);
        Double_t crossSection = hist3->GetBinContent(bin) * 524/100000;
        XS3->Fill(value, crossSection);
    }

    // Create a canvas to draw histograms
    TCanvas *c1 = new TCanvas("c1", "Cross-section Plot", 1200, 800);

    // Set histogram styles
    XS1->SetLineColor(kRed);
    XS2->SetLineColor(kBlue);
    XS3->SetLineColor(kGreen);

    // Draw histograms
    XS1->Draw("");
    XS2->Draw("same");
    XS3->Draw("same");

    // Add legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(XS1, "0n0n", "l");
    legend->AddEntry(XS2, "0nXn", "l");
    legend->AddEntry(XS3, "XnXn", "l");
    legend->Draw();

    // Draw canvas
    c1->Draw();
    c1->SaveAs("Classes_of_cross_section.png");




}
