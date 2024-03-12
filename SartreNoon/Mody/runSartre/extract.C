#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include <cmath>

TH1D* exRapidity(TTree* tree, TBranch* branch);
TH1D* exEnergy(TTree* tree, TBranch* branch);
TH1D* exSection(TH1D* histogram);

void extract()
{
    TFile *inputFile = new TFile("example_Pb.root", "READ");

    if (!inputFile->IsOpen() || inputFile->IsZombie()) {
        std::cerr << "Error: Unable to open ROOT file '" << inputFile << "'.\n";
        return;
    }

    TTree *tree = dynamic_cast<TTree*>(inputFile->Get("tree"));

    if (!tree) {
        std::cerr << "Error: TTree 'tree' not found in the ROOT file.\n";
        inputFile->Close();
        return;
    }

    TFile *outputRootFile = new TFile("extract.root", "RECREATE");

    exRapidity(tree, tree->GetBranch("vm"))->Write();
    exEnergy(tree, tree->GetBranch("Gamma"))->Write();
    exSection(exRapidity(tree, tree->GetBranch("vm")))->Write();

    outputRootFile->Close();
    inputFile->Close();
}

TH1D* exRapidity(TTree* tree, TBranch* branch)
{
    TLorentzVector* lorentzVector = new TLorentzVector();
    branch->SetAddress(&lorentzVector);

    std::vector<double> rapidityValues;

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        branch->GetEntry(entry);
        Double_t rapidity = lorentzVector->Rapidity();
        rapidityValues.push_back(rapidity);
    }

    Double_t minRapidity = *std::min_element(rapidityValues.begin(), rapidityValues.end());
    Double_t maxRapidity = *std::max_element(rapidityValues.begin(), rapidityValues.end());

    TH1D* histogram = new TH1D("Rapidity", "Rapidity", 100, minRapidity, maxRapidity);

    for (const auto& rapidity : rapidityValues) {
        Double_t roundedRapidity = std::round(rapidity * 100.0) / 100.0;
        histogram->Fill(roundedRapidity);
    }

    delete lorentzVector;
    return histogram;
}

TH1D* exEnergy(TTree* tree, TBranch* branch)
{
    TLorentzVector* lorentzVector = new TLorentzVector();
    branch->SetAddress(&lorentzVector);

    std::vector<double> EnergyValues;

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        branch->GetEntry(entry);
        Double_t Energy = lorentzVector->Energy();
        EnergyValues.push_back(Energy);
    }

    Double_t minEnergy = *std::min_element(EnergyValues.begin(), EnergyValues.end());
    Double_t maxEnergy = *std::max_element(EnergyValues.begin(), EnergyValues.end());

    TH1D* histogram = new TH1D("PhotonK", "PhotonK", 100, minEnergy, maxEnergy);

    for (const auto& Energy : EnergyValues) {
        Double_t roundedEnergy = std::round(Energy * 100.0) / 100.0;
        histogram->Fill(roundedEnergy);
    }

    delete lorentzVector;
    return histogram;
}

TH1D* exSection(TH1D* histogram)
{
    TH1D* xSectionHistogram = new TH1D("XSection", "X Section", histogram->GetNbinsX(), histogram->GetXaxis()->GetXmin(), histogram->GetXaxis()->GetXmax());

    double xSectionValue = 0;

    for (int binIndex = 1; binIndex <= histogram->GetNbinsX(); ++binIndex) {
        double binContent = histogram->GetBinContent(binIndex);
        double binWidth = histogram->GetBinWidth(binIndex);

        xSectionValue = (binContent * 524 / 100) / binWidth;

        xSectionHistogram->SetBinContent(binIndex, xSectionValue);
    }

    xSectionHistogram->Write();
    return xSectionHistogram;
}

