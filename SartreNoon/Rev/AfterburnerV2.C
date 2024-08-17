/**
 * 
 * This is AfterburnerV2.c i.e. the 2nd itiration of the AFterburner.C
 * It is basically the same as previous version but with clearer notations and comments for user friendliness
 * 
*/


//Importing necessary packages

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
#include "TRandom.h"

#include "NeutronGenerator.h"
#endif

//Defining different fucntions to print the graphs

void DrawXSection(TH1D* XS1, TH1D* XS2, TH1D* XS3, TH1D* XST); // For drawing crss-section

void DrawEnergy(TH1D* E1, TH1D* E2, TH1D* E3); // for drawing the energies of vm 

void DrawPseudoRapidity(TH1D* Eta1,TH1D* Eta2,TH1D* Eta3); // for drawing pseudo rapidity of the vm

void DrawNeutronRapidity(TH1D* NeutronY2, TH1D* NeutronY3); // for drawing the rapidity of neutrons

void DrawNeutronPseudorapidity(TH1D* NeuttonEta2, TH1D* NeutronEta3); // for drawning the pseudo rapidities of the neutrons

void DrawNeutronEnergy(TH1D* NeutronE2, TH1D* NeutronE3); // for drawning the neutron energies

void DrawTotNeutronEnergy(TH1D *TotE); //for the neutron energy depositd graph

//Defining the main function

void AfterburnerV2()
{
    /**
     * creating a neutron generator object that will be used for neutron related calculations
     */

    NeutronGenerator *gen = new NeutronGenerator(); //it is a pointer varibale so every thing will be done in the NeutronGenerator programme

    //assigning different parameters of the neutron generator

    gen->SetStoreQA();
    gen->SetStoreGeneratorFunctions();
    gen->SetHadronicInteractionModel(NeutronGenerator::kGlauber);
    gen->Initialize();
    gen->SetRunMode(NeutronGenerator::kInterface);
    gen->ReadENDF(kTRUE);
    gen->LoadENDF("hENDF.root");
    gen->Setup();


    // Rminder: you have to maually set the beamgamma in the NeutronGenerator.CXX for each run with different parameter sets


    /**
     * Creating a TTree file and acessing the root objects fromt he example_Pb.root file 
     * Mind that both the Afterburner and the example_Pb.root files should be in the same location
    */
    TFile *file = new TFile("example_Pb.root","READ");

    TTree *tree = dynamic_cast<TTree*>(file->Get("tree"));

    Long64_t nE = tree->GetEntries(); // Get all the entries of the roor file

    TLorentzVector* lv = new TLorentzVector; // a LorentzVector object 

    TBranch *bc = tree->GetBranch("vm"); // Accessing the vector meson branch

    bc->SetAddress(&lv); // pointint the LorentzVector to the vector meson branch to access that data

    /**
     * Manually created various vectors(Lists/arrays) to store the data for different parameter values for vector meson and the neutrons
     */

    std::vector<Double_t> onon; //srotes rapiity for 0n0n neutron Class
    std::vector<Double_t> onXn; // stores rapidity for 0nXn neutron Class
    std::vector<Double_t> XnXn; // stores rapidity for XnXN neutron Class

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

    std::vector<Double_t> NeutronTotE;

    std::vector<Double_t> onEmptyE;

    std::vector<Double_t> onEmptyEta;
    
    std::vector<Double_t> onEmptyY;
    
    std::vector<Double_t> VmBoostE;

    //looping over all the entries of the root file

    for(Long64_t index = 0; index < nE; index++)
    {
        bc->GetEntry(index);  // get specific entry umber

        Double_t y = lv->Rapidity(); // store the rapidity corresponding to that entry

        Double_t eta= lv->Eta(); // store the pseudo rapidity corresponding to that entry

        Double_t k = 0.5*lv->M()*TMath::Exp(TMath::Abs(y)); //clculates the photon energy for vm production for that entry

        std::vector<Int_t> nNumbers = gen->runSartreNoon(k); // Gives number of neutrons in both beams


        if(nNumbers[0]==0 && nNumbers[1] == 0) 
        {
            onon.push_back(y);
            ononEta.push_back(eta);
            on.push_back(index);
            ononE.push_back(lv->Energy());
        }
        
        else if(nNumbers[0]==0 || nNumbers[1]==0 )
        {
            std::vector<Double_t> lnEmptyY;
            std::vector<Double_t> lnEmptyEta;
            std::vector<Double_t> lnEmptyE;

            onXn.push_back(y);
            onXnEta.push_back(eta);
            ln.push_back(index);
            onXnE.push_back(lv->Energy());
            
            gen->createSartreNeutrons(nNumbers[0],nNumbers[1],lnNeutronE,lnNeutronEta,lnNeutronY);
            
            if(gRandom->Rndm()<0.5)gen->createSartreNeutrons(nNumbers[0],nNumbers[1],lnEmptyE,lnEmptyEta,lnEmptyY);
            else gen->createSartreNeutrons(nNumbers[1],nNumbers[0],lnEmptyE,lnEmptyEta,lnEmptyY);

            Double_t TNE = 0;
            for(Int_t i=0;i<lnEmptyE.size();i++)
            {
                TNE +=lnEmptyE[i];
            }
            NeutronTotE.push_back(TNE);
        }


    }

}