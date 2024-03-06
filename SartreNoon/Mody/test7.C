#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include <cmath>

void createRoundedRapidityHistogram(const char* inputFile, const char* outputFile) {
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
    std::vector<double> rapidityValues;

    // Loop over all entries in the tree and record rapidity values
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        branch->GetEntry(entry);  // Get the TLorentzVector data for this entry
        Double_t rapidity = lorentzVector->Rapidity();  // Access rapidity data
        rapidityValues.push_back(rapidity);  // Record the rapidity value
    }
    Double_t minRapidity = *std::min_element(rapidityValues.begin(), rapidityValues.end());
    Double_t maxRapidity = *std::max_element(rapidityValues.begin(), rapidityValues.end());

    // Create a 1D histogram for storing rounded rapidity values
    TH1D *histogram = new TH1D("Rapidity", "Rapidity", 100, minRapidity, maxRapidity);

    // Round off the recorded rapidity values and fill the histogram
    for (const auto& rapidity : rapidityValues) {
        Double_t roundedRapidity = std::round(rapidity * 100.0) / 100.0;
        histogram->Fill(roundedRapidity);
    }

    // Open a new ROOT file for writing
    TFile *outputRootFile = new TFile(outputFile, "RECREATE");

    // Write the histogram to the output ROOT file
    histogram->Write(); // why??

    // Close the files
    outputRootFile->Close();
    file->Close();
}

void test7() {
    createRoundedRapidityHistogram("example_Pb.root", "Rapidity_Pb.root");
    // Modify the 1st parameter to your exampl.root file 
}
