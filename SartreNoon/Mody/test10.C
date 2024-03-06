#include "TFile.h"
#include "TH1D.h"

void calculateAndSaveXSectionForAllBins(const char* inputFile, const char* outputFile) {
    // Open the input ROOT file
    TFile* inputFilePtr = TFile::Open(inputFile, "READ");

    if (!inputFilePtr || inputFilePtr->IsZombie()) {
        std::cerr << "Error: Could not open input file " << inputFile << std::endl;
        return;
    }

    // Get the histogram from the input file (assuming the histogram name is "histRapidity")
    TH1D* histogram = dynamic_cast<TH1D*>(inputFilePtr->Get("Rapidity"));

    if (!histogram) {
        std::cerr << "Error: Histogram not found in the input file." << std::endl;
        inputFilePtr->Close();
        return;
    }

    // Create a new ROOT file for the output
    TFile* outputFilePtr = TFile::Open(outputFile, "RECREATE");

    // Create a new histogram to store the cross-section values
    TH1D* xSectionHistogram = new TH1D("XSection", "X Section", histogram->GetNbinsX()+1, histogram->GetXaxis()->GetXmin(), histogram->GetXaxis()->GetXmax());

    //histogram->SetBinContent(1,0);
    double xSectionValue = 0 ;
    xSectionHistogram->SetBinContent(0,0);
    // Loop over all bins
    for (int binIndex = 1; binIndex <= histogram->GetNbinsX(); ++binIndex) {
        // Get the content and width of the bin
        //double prev = histogram->GetBinContent(binIndex-1);
        double binContent = histogram->GetBinContent(binIndex);
        double binWidth = 0.08;

        // Calculate the cross-section value for the bin
        
        xSectionValue = (binContent * 524/0.08);

        // Set the value in the new histogram
        xSectionHistogram->SetBinContent(binIndex, xSectionValue);
    }

    // Write the new histogram to the output file
    xSectionHistogram->Write();

    // Close the files
    inputFilePtr->Close();
    outputFilePtr->Close();

    std::cout << "XSection values calculated and saved to " << outputFile << std::endl;
}

void test10() {
    const char* inputFile = "Rapidity_Pb.root";
    const char* outputFile = "XSection.root";

    calculateAndSaveXSectionForAllBins(inputFile, outputFile);
}
