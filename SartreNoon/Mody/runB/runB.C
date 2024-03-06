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
#endif

void computeModelBreakups();

void runB(){

  computeModelBreakups();

}

void computeModelBreakups(){

#if defined(__CINT__)
  gROOT->LoadMacro("NeutronGenerator.cxx+g");
#endif

  TClonesArray *fParticles = new TClonesArray("TParticle", 200);
  TTree *fEventTree = new TTree("fEventTree", "fEventTree");
  fEventTree ->Branch("fParticles", &fParticles);
  TClonesArray *fSLparticles = new TClonesArray("TParticle", 200);

  TMap *fSLparams = new TMap();

    gROOT->LoadMacro("NeutronGenerator.cxx+g");


    TClonesArray *fNGparticles = NULL;

    //NeutronGenerator *gen = new NeutronGenerator();

   NeutronGenerator *gen = new NeutronGenerator();
   gen->SetStoreQA();
   gen->SetStoreGeneratorFunctions();
  gen->SetHadronicInteractionModel(NeutronGenerator::kHardSphere); 
  gen->Initialize(); 
  gen->SetRunMode(NeutronGenerator::kInterface);  
  gen->ReadENDF(kFALSE);
  gen->Setup();
  //gen->Run(10000);
  
    TDatabasePDG pdgData;
    Double_t VMrapidity = 0;
    Double_t VMmass = 3.09;// jebspi mass
    Double_t photonK = 0;
    Double_t fRapMin = -4.0;
    Double_t fRapMax = 4.0;

    UInt_t nEvents = 10000;
    cout<<"Running production"<<endl; 

     TFile *inputFile2 = new TFile("Energy_Pb.root","READ");
     TH1D *Photonk = (TH1D*)inputFile2->Get("PhotonK");
     //Int_t nEvents = Photonk->GetNbinsX()+1;

    for(Int_t iEvent=0 ; iEvent<=nEvents;iEvent++)
    { 

      photonK = Photonk->GetBinContent(iEvent);
      //std::cout<<photonK<<std::endl;

      VMrapidity = TMath::Abs(TMath::Log(2*photonK/3.09)); //  change 3.09 to mass variable

      gen->GenerateEvent(photonK);

      gen->FinishEvent();


    }
    
    gen->FinishProduction();
    cout<<"Finish production"<<endl;
  
    // TFile *fOutputFile = new TFile("SLoutput.root","RECREATE");
    // fEventTree->Write();
    // fOutputFile->Close(); 
}
