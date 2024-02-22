#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include <cmath>

void createRoundedEnergyHistogram(const char* inputFile, const char* outputFile) {
    // Open the ROOT file
    TFile *file = new TFile(inputFile, "READ");

    if (!file->IsOpen() || file->IsZombie()) {
        std::cerr << "Error: Unable to open ROOT file '" << inputFile << "'.\n";
        return;
    }

    // Get the TTree named "tree" from the file
    TTree *tree = dynamic_cast<TTree*>(file->Get("tree"));

    if (!tree) {
        std::cerr << "Error: TTree 'tree' not found in the ROOT file.\n";
        file->Close();
        return;
    }

    // Create a TLorentzVector to store the data
    TLorentzVector* lorentzVector = new TLorentzVector();

    // Set up the TBranch named "vm" from the TTree
    TBranch *branch = tree->GetBranch("vm");

    if (!branch) {
        std::cerr << "Error: TBranch 'vm' not found in the TTree.\n";
        file->Close();
        return;
    }

    // Set the branch address to the TLorentzVector
    branch->SetAddress(&lorentzVector);

    // Create a vector to store recorded rapidity values
    std::vector<double> EnergyValues;

    // Loop over all entries in the tree and record rapidity values
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        branch->GetEntry(entry);  // Get the TLorentzVector data for this entry
        Double_t Energy = lorentzVector->Et();  // Access rapidity data
        EnergyValues.push_back(Energy);  // Record the rapidity value
    }

    // Create a 1D histogram for storing rounded rapidity values
    TH1D *histogram = new TH1D("PhotonK", "PhotonK", 10000, 0, 5);

    // Round off the recorded rapidity values and fill the histogram
    for (const auto& Energy : EnergyValues) {
        Double_t roundedEnergy = std::round(Energy * 10000.0) / 10000.0;
        histogram->Fill(roundedEnergy);
    }

    // Open a new ROOT file for writing
    TFile *outputRootFile = new TFile(outputFile, "RECREATE");

    // Write the histogram to the output ROOT file
    histogram->Write(); // why??

    // Close the files
    outputRootFile->Close();
    file->Close();
}

void test8() {
    createRoundedEnergyHistogram("example_Pb.root", "Energy_Pb.root");
    // Modify the 1st parameter to your exampl.root file 
}
