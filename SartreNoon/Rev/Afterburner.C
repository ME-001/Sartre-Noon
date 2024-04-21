
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

    std::vector<Double_t> ononEta;
    std::vector<Double_t> onXnEta;
    std::vector<Double_t> XnXnEta;

    std::vector<Int_t> on;
    std::vector<Int_t> ln;
    std::vector<Int_t> Xn;


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
        }
        else if(nNumbers[0]==0 || nNumbers[1]==0 )
        {
            onXn.push_back(y);
            onXnEta.push_back(eta);
            ln.push_back(index);
        }
        else if(nNumbers[0]!=0 & nNumbers[1]!=0) 
        {
            XnXn.push_back(y);
            XnXnEta.push_back(eta);
            Xn.push_back(index);
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

    // Create histograms for each rapidity list
    TH1F *hist1 = new TH1F("hist1", "Rapidity Distribution 1", 1000, -10.0, 10.0);
    TH1F *hist2 = new TH1F("hist2", "Rapidity Distribution 2", 1000, -10.0, 10.0);
    TH1F *hist3 = new TH1F("hist3", "Rapidity Distribution 3", 1000, -10.0, 10.0);

    TH1F *Eta1 = new TH1F("Eta1", "PseudoRapidity Distribution 1", 1000, -10.0, 10.0);
    TH1F *Eta2 = new TH1F("Eta2", "PseudoRapidity Distribution 2", 1000, -10.0, 10.0);
    TH1F *Eta3 = new TH1F("Eta3", "PseudoRapidity Distribution 3", 1000, -10.0, 10.0);

    TH1F *E1 = new TH1F("E1", "PseudoRapidity Distribution 1", 1000, -10.0, 10.0);
    TH1F *E2 = new TH1F("E2", "PseudoRapidity Distribution 2", 1000, -10.0, 10.0);
    TH1F *E3 = new TH1F("E3", "PseudoRapidity Distribution 3", 1000, -10.0, 10.0);

    // TLorentzVector *lv;

    // TTree *tre = dynamic_cast<TTree*>(file->Get("tree"));

    // TLorentzVector* l = new TLorentzVector();

    // TBranch *b = tre->GetBranch("vm");

    // bc->SetAddress(&lv);

    // Fill histograms with rapidity values
    for (Int_t i; i<onon.size();i++) {
        
        // bc->GetEntry(on[i]);

        hist1->Fill(onon[i]);
        Eta1->Fill(ononEta[i]);
        // E1->Fill(lv->E());

    }
    for (Int_t i; i<onXn.size();i++) {
        
        // bc->GetEntry(ln[i]);
        
        hist2->Fill(onXn[i]);
        Eta2->Fill(onXn[i]);
        // E2->Fill(lv->E());
    }
    // for (Int_t i; i<onXn.size();i++) {
        
    //     // bc->GetEntry(ln[i]);
        
    //     // hist2->Fill(onXn[i]);
    //     Eta2->SetBinContent(onXn[i],0.5*Eta2->GetBinContent(onXn[i]));
    //     // E2->Fill(lv->E());
    // }
    for (Int_t i; i<XnXn.size();i++) {

        // bc->GetEntry(Xn[i]);

        hist3->Fill(XnXn[i]);
        Eta3->Fill(XnXnEta[i]);
        // E3->Fill(lv->E());
    }

    // Create histograms for cross-section distribution
    TH1F *XS1 = new TH1F("XS1", "Xsection Distribution 1", 1000, -10.0, 10.0);
    TH1F *XS2 = new TH1F("XS2", "Xsection Distribution 2", 1000, -10.0, 10.0);
    TH1F *XS3 = new TH1F("XS3", "Xsection Distribution 3", 1000, -10.0, 10.0);
    TH1F *XST = new TH1F("XST", "Total Cross section", 1000, -10.0, 10.0);

    // // Fill cross-section histograms
    // for (double value : onon) {
    //     Int_t bin = hist1->FindBin(value);
    //     Double_t crossSection = hist1->GetBinContent(bin);
    //     XS1->Fill(value, crossSection);
    // }
    // for (double value : onXn) {
    //     Int_t bin = hist2->FindBin(value);
    //     Double_t crossSection = hist2->GetBinContent(bin) ;
    //     XS2->Fill(value, crossSection);
    // }
    // for (double value : XnXn) {
    //     Int_t bin = hist3->FindBin(value);
    //     Double_t crossSection = hist3->GetBinContent(bin) ;
    //     XS3->Fill(value, crossSection);
    // }
    // XS1 = hist1;
    // XS2 = hist2;
    // XS3 = hist3;

    for(Int_t entry = 0; entry<1000; entry++)
    {
        XS1->SetBinContent(entry,hist1->GetBinContent(entry) * 524 * 1000 /100000);
        XS2->SetBinContent(entry,hist2->GetBinContent(entry) * 524 * 1000 /100000);
        XS3->SetBinContent(entry,hist3->GetBinContent(entry) * 524 * 1000 /100000);
        XST->SetBinContent(entry,XS1->GetBinContent(entry)+XS2->GetBinContent(entry)+XS3->GetBinContent(entry));
    }

    // Create a canvas to draw histograms
    TCanvas *c1 = new TCanvas("c1", "Cross-section Plot", 1200, 800);

    // Set histogram styles
    XS1->SetLineColor(kRed);
    XS2->SetLineColor(kBlue);
    XS3->SetLineColor(kGreen);
    XST->SetLineColor(kBlack);

    // Draw histograms
    XS1->GetXaxis()->SetRangeUser(-4,4);
    XS1->GetYaxis()->SetRangeUser(0,5000);
    XS1->Draw("L");
    XS2->GetXaxis()->SetRangeUser(-4,4);
    XS2->GetYaxis()->SetRangeUser(0,5000);
    XS2->Draw("SAME L");
    XS3->GetXaxis()->SetRangeUser(-4,4);
    XS3->GetYaxis()->SetRangeUser(0,5000);
    XS3->Draw("SAME L");
    XST->GetXaxis()->SetRangeUser(-4,4);
    XST->GetYaxis()->SetRangeUser(0,5000);
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
    c1->SaveAs("Classes_of_cross_section.png");

    
    // delete legend;
    

    TCanvas *c2 = new TCanvas("c2", "Pseudo Rapidity", 1200, 800);

    
    Eta1->SetLineColor(kRed);
    Eta2->SetLineColor(kBlue);
    Eta3->SetLineColor(kGreen);

    // Draw histograms
    Eta1->GetYaxis()->SetRangeUser(0, 500);
    Eta1->Draw("L");
    Eta2->GetYaxis()->SetRangeUser(0, 500);
    Eta2->Draw("SAME L");
    Eta3->GetYaxis()->SetRangeUser(0, 500);
    Eta3->Draw("SAME L");

    // // Add legend
    // TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    // legend->AddEntry(Eta1, "0n0n", "l");
    // legend->AddEntry(Eta2, "0nXn", "l");
    // legend->AddEntry(Eta3, "XnXn", "l");
    // legend->Draw();

    c2->Draw();
    c2->SaveAs("Pseudo_rapidity.png");

    // std::cout<<"Pseudo Rapidity of 0n0n: "<<ononEta<<std::endl;
    // std::cout<<"Pseudo Rapidity of 0n0n: "<<onXnEta<<std::endl;
    // std::cout<<"Pseudo Rapidity of 0n0n: "<<XnXnEta<<std::endl;

    std::ofstream outFile("pseudo_rapidity.txt");

    if (!outFile.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return;
    }

    outFile << "Pseudo Rapidity of 0n0n:\tPseudo Rapidity of 0nXn:\tPseudo Rapidity of XnXn:" << std::endl;
    
    size_t maxSize = std::max({ononEta.size(), onXnEta.size(), XnXnEta.size()});

    for (size_t i = 0; i < maxSize; ++i) {
        // Output values from ononEta
        if (i < ononEta.size()) {
            outFile << ononEta[i] << "\t";
        } else {
            outFile << "\t"; // Empty column if vector is shorter
        }

        // Output values from onXnEta
        if (i < onXnEta.size()) {
            outFile << onXnEta[i] << "\t";
        } else {
            outFile << "\t"; // Empty column if vector is shorter
        }

        // Output values from XnXnEta
        if (i < XnXnEta.size()) {
            outFile << XnXnEta[i];
        }

        outFile << std::endl; // Move to the next line after each set of values
    }

    outFile.close();

    std::cout << "Data written to pseudo_rapidity.txt" << std::endl;

}
