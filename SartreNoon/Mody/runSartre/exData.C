/*
 * This program is for collecting and storing the data event by event
 * 
 * Te other progran "extract.C" collects the data and stores in histogram. 
 * 
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

void exData()
{
    TFile *file = new TFile("example_Pb.root","READ");

    TTree *tree = dynamic_cast<TTree*>(file->Get("tree"));

    TLorentzVector* lorentzVector = new TLorentzVector();
    TLorentzVector* lovec = new TLorentzVector();

    TBranch *branch = tree->GetBranch("vm");
    //TBranch *bch = tree->GetBranch("gamma");

    branch->SetAddress(&lorentzVector);
    //bch->SetAddress(&lovec);

    Long64_t nEntries = tree->GetEntries();

    
    std::vector<double> rapidityValues;
    std::vector<double> energyValues;

    for (Long64_t entry = 0; entry < nEntries; entry++) {
        branch->GetEntry(entry);
        //bch->GetEntry(entry);

        Double_t rapidity = lorentzVector->Rapidity();
        rapidityValues.push_back(rapidity);

        Double_t energy = lorentzVector->Energy();
        energyValues.push_back(energy);
        
        //PhotonK->Fill(energy);
        //Rapidity->Fill(rapidity);

    }
    delete lorentzVector;
    delete tree;
    file->Close();
    delete file;

    auto normalize = [](double& value) {value = std::round(value * 100.0) / 100.0;};

    // std::for_each(minRapidity, maxRapidity, normalize);
    // std::for_each(minEnergy, maxEnergy, normalize);

    for (auto &value : rapidityValues) normalize(value);
    for (auto &value : energyValues) normalize(value);

    Double_t maxRapidity = *std::max_element(rapidityValues.begin(), rapidityValues.end());
    Double_t minRapidity = *std::min_element(rapidityValues.begin(), rapidityValues.end());

    Double_t minEnergy = *std::min_element(energyValues.begin(), energyValues.end());
    Double_t maxEnergy = *std::max_element(energyValues.begin(), energyValues.end());

    #if defined(__CINT__)
        gROOT->LoadMacro("NeutronGenerator.cxx+g");
    #endif
    
    TClonesArray *fParticles = new TClonesArray("TParticle", energyValues.size());
    TTree *Tree = new TTree("Tree", "Tree");
    Tree ->Branch("fParticles", &fParticles);

    TClonesArray *fNGparticles = NULL;


    NeutronGenerator *gen = new NeutronGenerator();

    gen->SetStoreQA();
    gen->SetStoreGeneratorFunctions();
    gen->SetHadronicInteractionModel(NeutronGenerator::kGlauber);
    gen->Initialize();
    gen->SetRunMode(NeutronGenerator::kInterface);
    gen->ReadENDF(kTRUE);
    gen->LoadENDF("hENDF.root");
    gen->Setup();


    Double_t  VMR = 0;
    Double_t VMmass = 3.09; //  j/psi mass
    Double_t photonK = 0;

    cout<<"Running production"<<endl; 

    UInt_t nEvents = energyValues.size();

    for(Int_t iEvent = 0; iEvent < nEvents; ++iEvent)
    {   
        photonK = energyValues[iEvent];

        VMR = rapidityValues[iEvent];

        //NeutronGenerator::hRapidityVM->Fill(VMR);
        //NeutronGenerator::hMassVM->Fill(VMmass);

        gen->GenerateEvent(photonK);
        Int_t nTotalPart = 0; 
        fNGparticles = gen->ImportParticles();
        for(Int_t i = 0; i<fNGparticles->GetEntriesFast(); i++)
        {
            TParticle *part(dynamic_cast<TParticle*>(fNGparticles->At(i)));
            new((*fParticles)[nTotalPart++]) TParticle(*part);
        }
        gen->FinishEvent();
    
        Tree->Fill();
        fParticles->Clear("C");

        gen->FinishEvent();
    }

    gen->FinishProduction();
    cout<<"Finish production"<<endl;
    TFile *fOutputFile = new TFile("output.root","RECREATE");
    Tree->Write();


    //std::unordered_map<Double_t, Int_t> rapidityCounts;
    //std::unordered_map<Double_t, Int_t> energyCounts;

    //for(const auto &element : rapidityValues) rapidityCounts[element]+=1;
    //for(const auto &element :energyValues) energyCounts[element]+=1;

    TH1D *PhotonK = new TH1D("PhotonK","PhotonK",100,minEnergy,maxEnergy);
    TH1D *vmEnergy = new TH1D("vm Energy","Vm Energy",100,minEnergy,maxEnergy);
    TH1D *Rapidity = new TH1D("Rapidity","Rapidity",100,minRapidity,maxRapidity);

    

    for (Long64_t entry = 0; entry < nEntries; entry++) {
        
        //Double_t energy = energyValues[entry];
        Double_t rapidity = rapidityValues[entry];
        Double_t energy = 0.5*3.09*TMath::Exp(TMath::Abs(rapidity));
        Double_t vmE = energyValues[entry];


        PhotonK->Fill(energy);
        Rapidity->Fill(rapidity);
        vmEnergy->Fill(vmE);

    }

    TFile *file2 = new TFile("exData.root","RECREATE");

    PhotonK->Write();
    Rapidity->Write();
    vmEnergy->Write();


    TH1D *xSection = new TH1D("xSection","xSection",100,minRapidity,maxRapidity);

    const Int_t cs = 524;//cross section 524 nb
    const Int_t ne = 10000;//Number of events 1M

    Double_t xSectionValue = 0;

    for(Long64_t index = 1; index <=100; index++)
    {
        Double_t binValue = Rapidity->GetBinContent(index);
        //Double_t binWidth = Rapidity->GetBinWidth(index);

        xSectionValue = binValue * cs / ne ;//binWidth ;

        xSection->SetBinContent(index, xSectionValue);
    }

    xSection->Write();

    TCanvas* c1 = new TCanvas("c1","c1");
        Rapidity->Draw();
    TCanvas* c2 = new TCanvas("c2","c2");
        xSection->Draw();
    TCanvas* c3 = new TCanvas("c3","c3");
        PhotonK->Draw();
    TCanvas* c4 = new TCanvas("c4","c4");
        vmEnergy->Draw();

    c1->SaveAs("Rapidity.png");
    c2->SaveAs("xSection.png");
    c3->SaveAs("PhotonK.png");
    c4->SaveAs("vmEnergy.png");

    file2->Close();
    delete file2;
    fOutputFile->Close();

    


    //delete bch;
    //delete lovec;
    delete c1;
    delete c2;
    delete c3;
    delete c4;


    


}
