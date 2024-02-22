#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH1F.h"
#include <cmath>

void replaceXsection(const char* fil1, const char* fil2) {

    TFile *file2 = new TFile(fil2, "READ");

    TFile *file1 = new TFile(fil1, "UPDATE");

    TH1D *rehist = dynamic_cast<TH1D*>(file2->Get("Rapidity"));

    TH1D *hist = dynamic_cast<TH1D*>(file1->Get("xSectionHist"));

    //if (rehist && hist) {
        // Copy the content of rehist to hist
        hist->Reset();
        hist->Add(rehist);

        // Update the histogram in the file
        file1->cd();
        hist->Write("", TObject::kOverwrite);

        std::cout << "Histogram replaced successfully." << std::endl;
    // } else {
    //     std::cerr << "Error: Histograms not found in one or both files." << std::endl;
    // }

    file1->Close();
    file2->Close();
}

void test9() {
    replaceXsection("ExampleTheory.root", "Rapidity_Pb.root");
}
