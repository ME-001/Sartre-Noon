
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


void classes()
{
    #if defined(__CLING__)
        gROOT->LoadMacro("NeutronGenerator.cxx+g");
    #endif

    NeutronGenerator *gen = new NeutronGenerator();

    gen->SetStoreQA();
    gen->SetStoreGeneratorFunctions();
    gen->SetHadronicInteractionModel(NeutronGenerator::kGlauber);
    gen->Initialize();
    gen->SetRunMode(NeutronGenerator::kInterface);
    gen->ReadENDF(kTRUE);
    gen->LoadENDF("hENDF.root");
    gen->Setup();

    TFile *file = new TFile("example_Pb.root", "UPDATE");
    TTree *tree = dynamic_cast<TTree*>(file->Get("tree"));

    Long64_t nE = tree->GetEntries();

    TLorentzVector *lv;
    TBranch *bc = tree->GetBranch("vm");
    bc->SetAddress(&lv);

    // Create new TLorentzVectors for each branch
    TLorentzVector *onon = new TLorentzVector();
    TLorentzVector *onXn = new TLorentzVector();
    TLorentzVector *XnXn = new TLorentzVector();

    TTree *tre = new TTree("Clas","Classes");
    tre->Branch("onon", "TLorentzVector", &onon);
    tre->Branch("onXn", "TLorentzVector", &onXn);
    tre->Branch("XnXn", "TLorentzVector", &XnXn);

    Double_t Mv = 3.09;

    for(Long64_t index = 0; index < nE; index++)
    {
        bc->GetEntry(index);
        Double_t y = lv->Rapidity();
        Double_t eta = lv->Eta();
        Double_t k = 0.5 * Mv * TMath::Exp(TMath::Abs(y));

        std::vector<Int_t> nNumbers = gen->runSartreNoon(k);

        // Fill the TLorentzVectors based on nNumbers
        if(nNumbers[0] == 0 && nNumbers[1] == 0) 
        {
            *onon = *lv;  // Copy lv to onon
            *onXn = *lv;  // Copy lv to onXn
            *XnXn = *lv;  // Copy lv to XnXn
        }
        else if(nNumbers[0] == 0 || nNumbers[1] == 0)
        {
            *onon = *lv;  // Copy lv to onon
            *onXn = *lv;  // Copy lv to onXn
            *XnXn = *lv;  // Copy lv to XnXn
        }
        else if(nNumbers[0] != 0 && nNumbers[1] != 0) 
        {
            *onon = *lv;  // Copy lv to onon
            *onXn = *lv;  // Copy lv to onXn
            *XnXn = *lv;  // Copy lv to XnXn
        }
        else
        {
            std::cout << "Error in the Event" << std::endl;
        }

        // Fill the tree
        tre->Fill();
    }

    tre->Write();
    file->Close();

    delete gen;  // Clean up
    delete file;
    delete tree;
    delete tre;
    delete onon;
    delete onXn;
    delete XnXn;
}
