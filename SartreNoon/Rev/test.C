#include "TTree.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1F.h" // Include TH1F for histograms
#include <random>

void test()
{
    // Open a ROOT file for writing
    TFile *file = new TFile("test.root", "RECREATE");
    // Create a TTree
    TTree *tree = new TTree("tree", "My tree");

    // Create a TClonesArray to hold arrays of TLorentzVector objects
    TClonesArray *arrayOfArrays = new TClonesArray("TClonesArray", 5); // Adjust the size as needed
    // Add a branch to the TTree to hold the TClonesArray
    tree->Branch("ArrayOfArraysBranch", "TClonesArray", &arrayOfArrays, 32000, 0);

    // Create histograms for px, py, pz, and e
    TH1F *h_px = new TH1F("h_px", "Distribution of px", 100, -1.0, 1.0);
    TH1F *h_py = new TH1F("h_py", "Distribution of py", 100, -1.0, 1.0);
    TH1F *h_pz = new TH1F("h_pz", "Distribution of pz", 100, -1.0, 1.0);
    TH1F *h_energy = new TH1F("h_energy", "Distribution of energy", 100, 0.0, 2.0);

    // Use std::random_device to seed the random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(-1.0, 1.0); // Adjust the range as needed

    // Loop through each array in the array of arrays
    for (Int_t index = 0; index < arrayOfArrays->GetSize(); index++)
    {
        TClonesArray *particleArray = new TClonesArray("TLorentzVector", 8); // Adjust the size as needed
        for (Int_t particleindex = 0; particleindex < particleArray->GetSize(); particleindex++)
        {
            // Generate random values for px, py, pz, and e within the specified range
            float px = dis(gen);
            float py = dis(gen);
            float pz = dis(gen);
            float energy = dis(gen);

            // Create a TLorentzVector object and add it to the particleArray
            TLorentzVector *myParticle = new ((*particleArray)[particleindex]) TLorentzVector(px, py, pz, energy);

            // Fill histograms
            h_px->Fill(px);
            h_py->Fill(py);
            h_pz->Fill(pz);
            h_energy->Fill(energy);
        }
        // Add the particleArray to the arrayOfArrays
        new ((*arrayOfArrays)[index]) TClonesArray(*particleArray);
    }

    // Fill the TTree with data
    tree->Fill();

    // Write histograms and TTree to the ROOT file
    h_px->Write();
    h_py->Write();
    h_pz->Write();
    h_energy->Write();
    tree->Write();

    // Close the file
    file->Close();
}
