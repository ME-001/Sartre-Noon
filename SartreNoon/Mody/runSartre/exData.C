/*
 * This program is for collecting and storing the data event by event
 * 
 * Te other progran "extract.C" collects the data and stores in histogram. 
 * 
*/



#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include <cmath>
#include "TCanvas.h"
#include "TGraph.h"

void exData()
{
    TFile *file = new TFile("example_Pb.root","READ");

    TTree *tree = dynamic_cast<TTree*>(file->Get("tree"));

    TLorentzVector* lorentzVector = new TLorentzVector();

    TBranch *branch = tree->GetBranch("vm");

    branch->SetAddress(&lorentzVector);

    Long64_t nEntries = tree->GetEntries();

    

    
    std::vector<double> rapidityValues;
    std::vector<double> energyValues;

    for (Long64_t entry = 0; entry < nEntries; entry++) {
        branch->GetEntry(entry);

        Double_t rapidity = lorentzVector->Rapidity();
        rapidityValues.push_back(rapidity);

        Double_t energy = lorentzVector->Gamma();
        energyValues.push_back(energy);
        
        //PhotonK->Fill(energy);
        //Rapidity->Fill(rapidity);

    }

    auto normalize = [](double& value) {value = std::round(value * 1000.0) / 1000.0;};

    // std::for_each(minRapidity, maxRapidity, normalize);
    // std::for_each(minEnergy, maxEnergy, normalize);

    for (auto &value : rapidityValues) normalize(value);
    for (auto &value : energyValues) normalize(value);

    Double_t maxRapidity = *std::max_element(rapidityValues.begin(), rapidityValues.end());
    Double_t minRapidity = *std::min_element(rapidityValues.begin(), rapidityValues.end());

    Double_t minEnergy = *std::min_element(energyValues.begin(), energyValues.end());
    Double_t maxEnergy = *std::max_element(energyValues.begin(), energyValues.end());


    //std::unordered_map<Double_t, Int_t> rapidityCounts;
    //std::unordered_map<Double_t, Int_t> energyCounts;

    //for(const auto &element : rapidityValues) rapidityCounts[element]+=1;
    //for(const auto &element :energyValues) energyCounts[element]+=1;

    TH1D *PhotonK = new TH1D("PhotonK","PhotonK",1000,minEnergy,maxEnergy);
    TH1D *Rapidity = new TH1D("Rapidity","Rapidity",1000,minRapidity,maxRapidity);


    for (Long64_t entry = 0; entry < nEntries; entry++) {
        
        Double_t energy = energyValues[entry];
        Double_t rapidity = rapidityValues[entry];

        PhotonK->Fill(energy);
        Rapidity->Fill(rapidity);

    }

    TFile *file2 = new TFile("exData.root","RECREATE");

    PhotonK->Write();
    Rapidity->Write();


    TH1D *xSection = new TH1D("xSection","xSection",1000,minRapidity,maxRapidity);

    const Int_t cs = 524;//cross section 524 nb
    const Int_t ne = 1e6;//Number of events 1M

    Double_t xSectionValue = 0;

    for(Long64_t index = 1; index <=1000; index++)
    {
        Double_t binValue = Rapidity->GetBinContent(index);
        //Double_t binWidth = Rapidity->GetBinWidth(index);

        xSectionValue = binValue * cs / ne ;//binWidth ;

        xSection->SetBinContent(index, xSectionValue);
    }

    xSection->Write();

    file->Close();
    file2->Close();

    delete file;
    delete file2;
    delete lorentzVector;
    delete branch;
    


}
