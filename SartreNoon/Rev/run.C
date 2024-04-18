/**
 * This features of this program will be implemented on sartre
 * 
 * AS of now it is used for the program to run separately
*/



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

//using namespace std;

void run()
{

    #if defined(__CLINT__)
        gROOT->LoadMacro("NeutronGenerator.cxx+g");
    #endif
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

    for(Long64_t i = 0; i<nE; ++i)
    {
        bc->GetEntry(i);

        rp.push_back(lv->Rapidity());
    }

    delete lv;
    delete tree;
    file->Close();
    delete file;

    Double_t maxrp = *std::max_element(rp.begin(),rp.end());
    Double_t minrp = *std::max_element(rp.begin(),rp.end());

    

    TH1D *rpHist = new TH1D("Rapidity", "Rapidity",nB,minrp,maxrp);

    for(Int_t i = 0; i<nE; ++i)
    {
        rpHist->Fill(rp[i]);
    }

    // TFile *file2 = new TFile("temp.root","RECREATE");
    // rpHist->Write();

//================================================================================================================
    /**
     * Till now I have created a histogram of the data 
     * I Probably don't need to create histogram
     * Will check that later
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

    std::vector<double> dp;
    std::vector<double> phk;

    for(Int_t i = 0; i<10; ++i)
    {  

        std::ofstream outfile("output.txt", std::ios_base::app);

        Double_t y = rpHist->GetRandom();
        
        Double_t k = 0.5*mv*TMath::Exp(TMath::Abs(y));

        gen->GenerateEvent(k); //number of neutrons and then their kinematics

        //desired: those kinematics idaelly in lab frame (boosted)

        //plot or print

        gen->runSartreNoon(k,i);

    }

    // gen->FinishProduction();

    // file2->Close();
    // delete file2;



}