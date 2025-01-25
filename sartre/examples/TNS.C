/**
 * General convention is E= energy; Y,y = Rapidity; Eta= pseudo rapidity 
*/

//Importing necessary packages

#ifdef __CLING__

// Copied from runBreakup.C as is, from noon
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

#include "TROOT.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"


// Included by us (DSR, BPN) for our use-case
#include <vector>  
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <cmath>

#include "TMath.h"
#include "TDatabasePDG.h"
#include "TRandom.h"


#endif

void DrawTotNeutronEnergy(TH1D *TotE);

void TNS()
{
	TFile *file = new TFile("example_Pb.root","READ");              // The root file name should be "example_Pb.root" if not please change it here

	TTree *tree = dynamic_cast<TTree*>(file->Get("tree"));

	Long64_t nE = tree->GetEntries();  
	
	printf("Number of entries %lld", nE);
	
	//TLorentzVector* lv = new TLorentzVector; 
	//TBranch *bc = tree->GetBranch("neutrons1");                            // Accessing the vector meson branch

    	//bc->SetAddress(&lv);
    	
    	std::vector<Double_t> nESum;
    	//std::vector<Double_t> nEta;
    	
    	 TClonesArray *neutrons = nullptr;
    	tree->SetBranchAddress("neutrons2", &neutrons);

    	for (Long64_t i = 0; i < nE; i++) {
        	tree->GetEntry(i);
        	Double_t E = 0;

        	if (neutrons) {
            		Int_t size = neutrons->GetEntriesFast();
            		printf(" Event %lld: number of neutrons %d\n", i, size);
			
			
			
            		for (Int_t j = 0; j < size; j++) {
                	TLorentzVector *neutron = dynamic_cast<TLorentzVector*>(neutrons->At(j));
                		if (neutron) {
                    		//printf("  Neutron %d: Px = %.2f, Py = %.2f, Pz = %.2f, E = %.2f\n", j, neutron->Px(), neutron->Py(), neutron->Pz(), neutron->E());
                    		E += neutron->E();
                    		//nEta.push_back(neutron->Eta());
                		}
            		}

            		
        	} else {
            		printf(" Event %lld: neutrons array is null\n", i);
        		}
        printf(" Total neutron Energy: %f", E);
    	nESum.push_back(E*0.001);
    	}
    	file->Close();
    	
    	TH1D *TotE = new TH1D("TotE","Total Neutron Energy in each Event",1000,0,30000);

    	for(Int_t i=0; i<nESum.size();i++)
    	{
        TotE->Fill(nESum[i]);
    	}

    	DrawTotNeutronEnergy(TotE);
    }
    
    
    void DrawTotNeutronEnergy(TH1D *TotE)
{
    TCanvas *c7 = new TCanvas("c7","TOtal Neutron Energy per Event side 1",1200,800);

    TotE->GetXaxis()->SetTitle("Energy(GeV)");
    TotE->SetAxisRange(0,10000,"X");
    TotE->GetYaxis()->SetTitle("Counts");
    TotE->Draw("L");
    c7->SetLogy();
    c7->SaveAs("Total_Neutron_E_5.02_side1_ss.png");
}

