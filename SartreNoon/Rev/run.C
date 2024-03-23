/**
 * This features of this program will be implemented on sartre
 * 
 * AS of now it is used for the program to run separately
*/



#ifdef __CLING__

#include <iostream>
#include <stdio>

#include "TTree.h"
#include "<vector>"
#include "TFile.h"
#include "TROOT.h"
#include "NeutronGenerator.h"
#include "NeutronGenerator.cxx+g"
#include "TRandom.h"

#endif

//using namespace std;

void run()
{

    //Variables

    Double_t mv =3.09; //unit in Gev 

    Int_t nB = 1000;//number of bins of the histograms


    TFile *file = new TFile("example_Pb.root", "READ");

    TTree *tree = dynamic_cast<TTree*>(file->Get("tree"));

    TLorentzVector* lv = new TLorentzVector();

    TBranch *bc = tree->GetBranch("vm");


    bc->SetAddress(&lv);

    Long64_t nE = tree->GetEntries();

    std::vector<double> rp;
    std::vector<double> phk;

    for(Long64_t i = 0; i<nE; ++i)
    {
        bc->GetEntry(i);

        rp.push_back(lv->Rapidity());

        phk.push_back(0.5*mv*TMath::Exp(TMath::(Abs(lv->Rapidity()))));
    }

    delete lv;
    delete tree;
    file->Close();
    delete file;

    Double_t maxphk = *std::max_element(phk.begin(),phk.end());
    Double_t minphk = *std::max_element(phk.begin(),phk.end());

    Double_t maxrp = *std::max_element(rp.begin(),rp.end());
    Double_t minrp = *std::max_element(rp.begin(),rp.end());

    #if defined(__CLINT__)
        gROOT->LoadMacro("NeutronGenerator.cxx+g");
    #endif

    TH1D *phkHist = new TH1D("PhotonK","PhotonK",nB,minphk,maxphk);
    TH1D *rpHist = new TH1D("Rapidity", "Rapidity",nB,minrp,maxrp);

    for(Int_t i = 0; i<nE; ++i)
    {
        phkHist->Fill(phk[i]);

        rpHist->Fill(rp[i]);

    }

    TFile *file2 = new TFile("temp.root","RECREATE");

    phkHist->Write();
    rpHist->Write();

    file2->Close();
    delete file2;
//================================================================================================================
    /**
     * Till now I have created a histogram of the data 
     * Next I will select data random from the data.
    */
//=================================================================================

    NeutronGenerator *gen = new NeutronGenerator();

    gen->SetStoreQA();
    gen->SetStoreGeneratorFunctions();
    gen->SetHadronicInteractionModel(NeutronGenerator::kGlauber);
    gen->Initialize();
    gen->SetRunMode(NeutronGenerator::kInterface);
    gen->ReadENDF(kTRUE);
    gen->LoadENDF("hENDF.root");
    gen->Setup();

    TRandom rm;

    for(Int_t i = 0; i<nE; ++i)
    {   
        Int_t ei = rm.Integer(phk.size());
        Double_t phki = phk[ei];
        Double_t rpi = rp[ei];

        gen->GenerateEvent(phki);
        
    }

    






}