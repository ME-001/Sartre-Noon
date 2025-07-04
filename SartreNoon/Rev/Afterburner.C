/**
 * General convention is E= energy; Y,y = Rapidity; Eta= pseudo rapidity 
*/

//Importing necessary packages

#ifdef __CLING__

// Copied from runBreakup.C as is, from noon
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

// Included by us (DSR, BPN) for our use-case
#include <vector>  
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

void DrawXSection(TH1D* XS1, TH1D* XS2, TH1D* XS3, TH1D* XST);          // For drawing cross-section

void DrawEnergy(TH1D* E1, TH1D* E2, TH1D* E3);                          // for drawing the energies of vm 

void DrawPseudoRapidity(TH1D* Eta1,TH1D* Eta2,TH1D* Eta3);              // for drawing pseudo rapidity of the vm

void DrawNeutronRapidity(TH1D* NeutronY2, TH1D* NeutronY3);             // for drawing the rapidity of neutrons

void DrawNeutronPseudorapidity(TH1D* NeuttonEta2, TH1D* NeutronEta3);   // for drawning the pseudo rapidities of the neutrons

void DrawNeutronEnergy(TH1D* NeutronE2, TH1D* NeutronE3);               // for drawning the neutron energies

void DrawTotNeutronEnergy(TH1D *TotE);                                  //for the neutron energy depositd graph

//Defining the main function

void Afterburner()
{
    /**
     * creating a neutron generator object that will be used for neutron related calculations
     */

    NeutronGenerator *gen = new NeutronGenerator();                     //it is a pointer varibale so every thing will be done in the NeutronGenerator programme

    //assigning different parameters of the neutron generator
    // below were taken as is from runBreakup.C in noon files
    gen->SetStoreQA();                                                  // sets kStoreQA = True; refer neutronGenerator.h and neutronGenerator.cxx files
    gen->SetStoreGeneratorFunctions();                                  //similar to above
    gen->SetHadronicInteractionModel(NeutronGenerator::kGlauber);       //selects model for ion interaction
    gen->Initialize();                                                  //loads default data used in noon (photon flux, cross-section etc.)
    gen->SetRunMode(NeutronGenerator::kInterface);                      // runMode = kInterface is plain neutron generation
    gen->ReadENDF(kTRUE);                                               //enables reading of ENDF data tables
    gen->LoadENDF("hENDF.root");                                        // data table for energy distribution amongst created neutrons
    gen->Setup();


    /** !!!
     * reminder: you have to maually set the beamgamma in the NeutronGenerator.CXX for each run with different parameter sets
     * You can check for the beamgamma in the root file
        !!!
    */ 


    /**
     * Creating a TTree file and acessing the root objects fromt he example_Pb.root file 
     * Mind that both the Afterburner and the example_Pb.root files should be in the same location
    */
    TFile *file = new TFile("example_Pb.root","READ");              // The root file name should be "example_Pb.root" if not please change it here

    TTree *tree = dynamic_cast<TTree*>(file->Get("tree"));

    Long64_t nE = tree->GetEntries();                               // Get all the entries of the roor file

    TLorentzVector* lv = new TLorentzVector;                        // a LorentzVector object 

    TBranch *bc = tree->GetBranch("vm");                            // Accessing the vector meson branch

    bc->SetAddress(&lv);                                            // pointint the LorentzVector to the vector meson branch to access that data

    /**
     * Manually created various vectors(Lists/arrays) to store the data for different parameter values for vector meson and the neutrons
     * 0n0n -> No breakup neutrons from both beams
     * 0nXn -> Breakup neutrons from either one of two beams 
     * XnXn -> Breakup neutrons from both beams
     */


    std::vector<Double_t> onon; //srotes rapiity for 0n0n neutron Class
    std::vector<Double_t> onXn; // stores rapidity for 0nXn neutron Class
    std::vector<Double_t> XnXn; // stores rapidity for XnXN neutron Class

    // Pseudo-rapidities
    std::vector<Double_t> ononEta;      //stores pseudo rapidity 
    std::vector<Double_t> onXnEta;      //stores pseudo rapidity
    std::vector<Double_t> XnXnEta;      //stores pseudo rapidity

    // Indices
    std::vector<Int_t> on; 
    std::vector<Int_t> ln; 
    std::vector<Int_t> Xn; 

    // Energies VM
    std::vector<Double_t> ononE;
    std::vector<Double_t> onXnE;
    std::vector<Double_t> XnXnE;

    // 0n0n class Neutron kinematics
    std::vector<Double_t> onNeutronE;   // energy
    std::vector<Double_t> onNeutronEta; // pseudorapidity
    std::vector<Double_t> onNeutronY;   // rapidity
    
    // 0nXn class neutron kinematics
    std::vector<Double_t> lnNeutronE;
    std::vector<Double_t> lnNeutronEta;
    std::vector<Double_t> lnNeutronY;

    // XnXn class neutron kinematics
    std::vector<Double_t> XnNeutronE;
    std::vector<Double_t> XnNeutronEta;
    std::vector<Double_t> XnNeutronY;

    // ZDC Calorimeter; total energy of neutrons
    std::vector<Double_t> NeutronTotE;

    std::vector<Double_t> onEmptyE;

    std::vector<Double_t> onEmptyEta;
    
    std::vector<Double_t> onEmptyY;
    
    // Boosted vm energies
    std::vector<Double_t> VmBoostE;

    //looping over all the entries of the root file

    for(Long64_t index = 0; index < nE; index++)
    {
        bc->GetEntry(index);                                    // get specific entry umber

        Double_t y = lv->Rapidity();                            // store the rapidity corresponding to that entry

        Double_t eta= lv->Eta();                                // store the pseudo rapidity corresponding to that entry

        Double_t k = 0.5*lv->M()*TMath::Exp(TMath::Abs(y));     //clculates the photon energy for vm production for that entry

        std::vector<Int_t> nNumbers = gen->runSartreNoon(k);    // Gives number of neutrons in both beams


        if(nNumbers[0]==0 && nNumbers[1] == 0)                  // loop to collect data for the 0n0n class of events
        {
            onon.push_back(y);
            ononEta.push_back(eta);
            on.push_back(index);
            ononE.push_back(lv->Energy());
        }

        else if(nNumbers[0]==0 || nNumbers[1]==0 )              // loop to collect data for 0nXn events
        {
            // these are to store neutrons for individual events only and resets at each event
            std::vector<Double_t> lnEmptyY; 
            std::vector<Double_t> lnEmptyEta;
            std::vector<Double_t> lnEmptyE;

            onXn.push_back(y);
            onXnEta.push_back(eta);
            ln.push_back(index);
            onXnE.push_back(lv->Energy());
            
            gen->createSartreNeutrons(nNumbers[0],nNumbers[1],lnNeutronE,lnNeutronEta,lnNeutronY);

            // the below code is for collecting neutron eneries for a single beam side.

            /**
             * While plotting for ZDC calorimeter neutron energies, manually put the beamside to either 0 or 1 in the NeutronGenerator.cxx 
             * You can do this by changing the condition at line 310 of NeutronGenerator.cxx, i.e. the for loop.
             * comment the below lines of code if not calculating ZDC calorimeter data
             */
            gen->createSartreNeutrons(nNumbers[0],nNumbers[1],lnEmptyE,lnEmptyEta,lnEmptyY);

            Double_t TNE = 0; //total neutron energy (for ZDC Calorimeter)
            for(Int_t i=0;i<lnEmptyE.size();i++)
            {
                TNE +=lnEmptyE[i];
            }
            NeutronTotE.push_back(TNE);
        }

        else if(nNumbers[0]!=0 && nNumbers[1]!=0) 
        {
            // these are to store neutrons for individual events only and resets at each event
            std::vector<Double_t> XnEmptyY;
            std::vector<Double_t> XnEmptyE;
            std::vector<Double_t> XnEmptyEta;

            XnXn.push_back(y);
            XnXnEta.push_back(eta);
            Xn.push_back(index);
            XnXnE.push_back(lv->Energy());
            gen->createSartreNeutrons(nNumbers[0],nNumbers[1],XnNeutronE,XnNeutronEta,XnNeutronY);

            // the below code is for collecting neutron eneries for a single beam side.

            /**
             * While plotting for ZDC calorimeter neutron energies, manually put the beamside to either 0 or 1 in the NeutronGenerator.cxx 
             * You can do this by changing the condition at line 310 of NeutronGenerator.cxx, i.e. the for loop.
             * comment the below lines of code if not calculating ZDC calorimeter data
             */
            gen->createSartreNeutrons(nNumbers[0],nNumbers[1],XnEmptyE,XnEmptyEta,XnEmptyY);

            Double_t TNE = 0;
            for(Int_t i=0;i<XnEmptyE.size();i++)
            {
                TNE +=XnEmptyE[i];
            }
            NeutronTotE.push_back(TNE);
        }
        else std::cout<<"Error in the Event"<<std::endl;
    }

    /**close the used objects 
     * It's a good practice and prevent errors or data leaks on the code.
     */
    delete lv;
    delete tree;
    file->Close();
    delete file;

    // Print the number of events in each class of events.
    std::cout<<"0n0n number of events: "<<onon.size()<<std::endl;
    std::cout<<"0nXn number of events: "<<onXn.size()<<std::endl;
    std::cout<<"XnXn number of events: "<<XnXn.size()<<std::endl;

    // Histograms for different kinematics
    TH1D *NeutronEln = new TH1D("NeutronE2","Neutron energy",10000,0,5000); // neutron energy values for 0nXn class of events
    TH1D *NeutronEXn = new TH1D("NeutronE3","Neutron energy",10000,0,5000); // neutron energy values for XnXn class of events
   
    //    TH1D *NeutronEta0n = new TH1D("NeutronEta1", "Neutron Pseudo Rapidity",1000,-10,10);
    TH1D *NeutronEtaln = new TH1D("NeutronEta2", "Neutron Pseudo Rapidity",1000,-25,25);
    TH1D *NeutronEtaXn = new TH1D("NeutronEta3", "Neutron Pseudo Rapidity",1000,-25,25);

    //    TH1D *NeutronYXn = new TH1D("NeutronY1","Neutron Rapidity",1000,-6,6);
    TH1D *NeutronYln = new TH1D("NeutronY2","Neutron Rapidity",1000,-25,25);
    TH1D *NeutronYXn = new TH1D("NeutronY3","Neutron Rapidity",1000,-25,25);


    // Filling the histograms
    for(Int_t i=0; i<lnNeutronE.size(); i++)
    {
        NeutronEln->Fill(lnNeutronE[i]);
        NeutronEtaln->Fill(lnNeutronEta[i]);
        NeutronYln->Fill(lnNeutronY[i]);
    }

    for(Int_t i=0; i<XnNeutronE.size(); i++)
    {
        NeutronEXn->Fill(XnNeutronE[i]);
        NeutronEtaXn->Fill(XnNeutronEta[i]);
        NeutronYXn->Fill(XnNeutronY[i]);
    }

    // Plotting tthe Neutron energies
    DrawNeutronEnergy(NeutronEln,NeutronEXn);

    // printing number of neutrons in each class of events 
    std::cout<<"Number of Neutrons in onon: "<<onNeutronE.size()<<std::endl;
    std::cout<<"Number of Neutrons in onxn: "<<lnNeutronE.size()<<std::endl;
    std::cout<<"Number of Neutrons in xnxn: "<<XnNeutronE.size()<<std::endl;

    /**
     * This part is for plotting the ZDC calorimeter datarecreation. 
     * Please commnet it out when not is use, in use comment other plots and functions: to save time and unnecessary calculations
     * While using please plot for ither one side and not both in NeutronGenerator.cxx
     */
    TH1D *TotE = new TH1D("TotE","Total Neutron Energy in each Event",1000,0,30000);

    for(Int_t i=0; i<NeutronTotE.size();i++)
    {
        TotE->Fill(NeutronTotE[i]);
    }

    DrawTotNeutronEnergy(TotE);


    // vector meson kinematic histograms

    /**
     * The convention is  1 for onon class of events; 2 for onXn class of events; 3 for XnXn class of events; T for total
     */

    Double_t minRange = -6; // for rapidity values
    Double_t maxRange = 6;
    // Create histograms for each rapidity list
    TH1D *hist1 = new TH1D("hist1", "vm Rapidity Distribution ", 1000, minRange,maxRange);
    TH1D *hist2 = new TH1D("hist2", "vm Rapidity Distribution ", 1000, minRange, maxRange);
    TH1D *hist3 = new TH1D("hist3", "vm Rapidity Distribution ", 1000, minRange,maxRange);
    
    // Histogram for pseudo rapidity 
    TH1D *Eta1 = new TH1D("Eta1", "vm PseudoRapidity Distribution ", 1000,  -15,15);
    TH1D *Eta2 = new TH1D("Eta2", "vm PseudoRapidity Distribution ", 1000, -15,15);
    TH1D *Eta3 = new TH1D("Eta3", "vm PseudoRapidity Distribution ", 1000, -15,15);

    // Histogram for energy
    TH1D *E1 = new TH1D("E1", "vm Energy ", 1000, 0 ,300);
    TH1D *E2 = new TH1D("E2", "vm Energy ", 1000, 0, 300);
    TH1D *E3 = new TH1D("E3", "vm Energy ", 1000, 0, 300);


    // filling in the histograms with kinematic data
    
    // onon class of events
    for (Int_t i; i<onon.size();i++) {

        hist1->Fill(onon[i]);
        Eta1->Fill(ononEta[i]);
        E1->Fill(ononE[i]);

    }
    
    // onXn class of events
    for (Int_t i; i<onXn.size();i++) {
                
        hist2->Fill(onXn[i]);
        Eta2->Fill(onXnEta[i]);
        E2->Fill(onXnE[i]);
    }
    
    // XnXn class of events
    for (Int_t i; i<XnXn.size();i++) {

        hist3->Fill(XnXn[i]);
        Eta3->Fill(XnXnEta[i]);
        E3->Fill(XnXnE[i]);
    }

    // Create histograms for cross-section distribution
    TH1D *XS1 = new TH1D("XS1", "Cros-Section", 1000, minRange, maxRange);
    TH1D *XS2 = new TH1D("XS2", "Cros-Section", 1000, minRange,maxRange);
    TH1D *XS3 = new TH1D("XS3", "Cros-Section", 1000, minRange,maxRange);
    TH1D *XST = new TH1D("XST", "Cros-Section", 1000, minRange,maxRange);

    // Filling in the cross section data for different class of events
    for(Int_t entry = 0; entry<1000; entry++)
    {
        XS1->SetBinContent(entry,(hist1->GetBinContent(entry)) * 523 * 1000 /100000);
        XS2->SetBinContent(entry,(hist2->GetBinContent(entry)) * 523 * 1000 /100000);
        XS3->SetBinContent(entry,(hist3->GetBinContent(entry)) * 523 * 1000 /100000);
        XST->SetBinContent(entry,XS1->GetBinContent(entry)+XS2->GetBinContent(entry)+XS3->GetBinContent(entry));
    }
    
    /**
     * Plotting and Saving as png
    */
    DrawXSection(XS1,XS2,XS3,XST);
    DrawEnergy(E1,E2,E3);
    DrawPseudoRapidity(Eta1,Eta2,Eta3);

    /**
     * Plotting and Saving as png
    */
}



// Definitions for different fucntions

/**
 * The genral structure of plot is 
 * create a canvas
 * set X and Y axis ranges 
 * then draw
 * Additionally we created .png files of the canvas
 */

void DrawTotNeutronEnergy(TH1D *TotE)
{
    TCanvas *c7 = new TCanvas("c7","TOtal Neutron Energy per Event side 1",1200,800);

    TotE->GetXaxis()->SetTitle("Energy(GeV)");
    TotE->SetAxisRange(0,5000,"X");
    TotE->GetYaxis()->SetTitle("Counts");
    TotE->Draw("L");
    c7->SetLogy();
    c7->SaveAs("Total_Neutron_E_2.76_side2_ss.png");
}

void DrawNeutronRapidity( TH1D* NeutronY2, TH1D* NeutronY3)
{
    TCanvas *c4 = new TCanvas("c4", "Cross-section Plot", 1200, 800);

    // NeutronY1->SetLineColor(kRed);
    NeutronY2->SetLineColor(kBlue);
    NeutronY3->SetLineColor(kGreen);

    // NeutronY1->GetXaxis()->SetTitle("Rapidity");
    // NeutronY1->GetYaxis()->SetTitle("Counts");
    // NeutronY1->Draw("L");
    
    NeutronY3->GetXaxis()->SetTitle("Rapidity");
    NeutronY3->GetYaxis()->SetTitle("Counts");
    NeutronY3->Draw("L");
    NeutronY2->GetXaxis()->SetTitle("Rapidity");
    NeutronY2->GetYaxis()->SetTitle("Counts");
    NeutronY2->Draw("SAME L");
    TLegend *legend4 = new TLegend(0.7, 0.7, 0.9, 0.9);
    // legend36->AddEntry(E1, "0n0n", "l");
    legend4->AddEntry(NeutronY2, "0nXn", "l");
    legend4->AddEntry(NeutronY3, "XnXn", "l");
    // legend->AddEntry(EtaT, "Total", "l");
    legend4->Draw();

    c4->SaveAs("Neutron_Rapidity_Boosted_2.76.png");
}

void DrawNeutronPseudorapidity( TH1D* NeutronEta2, TH1D* NeutronEta3)
{
    TCanvas *c5 = new TCanvas("c5", "Cross-section Plot", 1200, 800);

    // NeutronEta1->SetLineColor(kRed);
    NeutronEta2->SetLineColor(kBlue);
    NeutronEta3->SetLineColor(kGreen);

    // NeutronEta1->GetXaxis()->SetTitle("Pseudo Rapidity");
    // NeutronEta1->GetYaxis()->SetTitle("Counts");
    // NeutronEta1->Draw("L");
    
    NeutronEta3->GetXaxis()->SetTitle("Pseudo Rapidity");
    NeutronEta3->GetYaxis()->SetTitle("Counts");
    NeutronEta3->Draw("L");
    NeutronEta2->GetXaxis()->SetTitle("Pseudo Rapidity");
    NeutronEta2->GetYaxis()->SetTitle("Counts");
    NeutronEta2->Draw("SAME L");
    TLegend *legend5 = new TLegend(0.7, 0.7, 0.9, 0.9);
    // legend36->AddEntry(E1, "0n0n", "l");
    legend5->AddEntry(NeutronEta2, "0nXn", "l");
    legend5->AddEntry(NeutronEta3, "XnXn", "l");
    // legend->AddEntry(EtaT, "Total", "l");
    legend5->Draw();

    c5->SaveAs("Neutron_Eta_Boosted_2.76.png");
}


void DrawNeutronEnergy( TH1D* NeutronE2, TH1D* NeutronE3)
{
    TCanvas *c6 = new TCanvas("c6", "Cross-section Plot", 1200, 800);

    // NeutronE1->SetLineColor(kRed);
    NeutronE2->SetLineColor(kBlue);
    NeutronE3->SetLineColor(kGreen);

    // NeutronE1->GetXaxis()->SetTitle("Energy");
    // NeutronE1->GetYaxis()->SetTitle("Counts");
    // NeutronE1->Draw("L");
    // NeutronE2->GetXaxis()->SetRangeUser(0,15);
    
    // NeutronE3->GetXaxis()->SetRangeUser(0,15);
    NeutronE3->GetXaxis()->SetTitle("Energy(GeV)");
    NeutronE3->GetYaxis()->SetTitle("Counts");
    NeutronE3->Draw("L");
    NeutronE2->GetXaxis()->SetTitle("Energy(GeV)");
    NeutronE2->GetYaxis()->SetTitle("Counts");
    NeutronE2->Draw("SAME L");
    TLegend *legend6 = new TLegend(0.7, 0.7, 0.9, 0.9);
    // legend36->AddEntry(E1, "0n0n", "l");
    legend6->AddEntry(NeutronE2, "0nXn", "l");
    legend6->AddEntry(NeutronE3, "XnXn", "l");
    // legend->AddEntry(EtaT, "Total", "l");
    legend6->Draw();

    c6->SaveAs("Neutron_Energy_Boosted_5.02.png");
}



void DrawXSection(TH1D* XS1, TH1D* XS2, TH1D* XS3, TH1D* XST)
{// Create a canvas to draw histograms
    TCanvas *c1 = new TCanvas("c1", "Cross-section Plot", 1200, 800);

    // Set histogram styles
    XS1->SetLineColor(kRed);
    XS2->SetLineColor(kBlue);
    XS3->SetLineColor(kGreen);
    XST->SetLineColor(kBlack);

    // Draw histograms
    Int_t minRapidity = -6;
    Int_t maxRapidity = 6;
    Int_t maxY = 2000;
    XS1->GetXaxis()->SetTitle("Rapidity");
    XS1->GetYaxis()->SetTitle("Cross section in nb");
    XS1->GetXaxis()->SetRangeUser(minRapidity,maxRapidity);
    XS1->GetYaxis()->SetRangeUser(0,maxY);
    XS1->Draw("L");
    XS3->GetXaxis()->SetTitle("Rapidity");
    XS3->GetYaxis()->SetTitle("Cross section in nb");
    XS2->GetXaxis()->SetRangeUser(minRapidity,maxRapidity);
    XS2->GetYaxis()->SetRangeUser(0,maxY);
    XS2->Draw("SAME L");
    XS3->GetXaxis()->SetTitle("Rapidity");
    XS3->GetYaxis()->SetTitle("Cross section in nb");
    XS3->GetXaxis()->SetRangeUser(minRapidity,maxRapidity);
    XS3->GetYaxis()->SetRangeUser(0,maxY);
    XS3->Draw("SAME L");
    XST->GetXaxis()->SetTitle("Rapidity");
    XST->GetYaxis()->SetTitle("Cross section in nb");
    XST->GetXaxis()->SetRangeUser(minRapidity,maxRapidity);
    XST->GetYaxis()->SetRangeUser(0,maxY);
    XST->Draw("Same L");

    // Add legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(XS1, "0n0n", "l");
    legend->AddEntry(XS2, "0nXn", "l");
    legend->AddEntry(XS3, "XnXn", "l");
    legend->AddEntry(XST, "Total", "l");
    legend->Draw();

    // Draw canvas
    c1->Draw();
    c1->SaveAs("Classes_of_cross_section_Boosted_2.76.png");
    
    
}



void DrawEnergy(TH1D* E1, TH1D* E2, TH1D* E3)
{
    TCanvas *c3 = new TCanvas("c3", "Energy", 1200, 800);

    E1->SetLineColor(kRed);
    E2->SetLineColor(kBlue);
    E3->SetLineColor(kGreen);


    E1->GetXaxis()->SetTitle("Energy(GeV)");
    E1->GetYaxis()->SetTitle("Counts");
    E1->GetYaxis()->SetRangeUser(0,2500);
    E1->Draw("L");
    E2->GetXaxis()->SetTitle("Energy(GeV)");
    E2->GetYaxis()->SetTitle("Counts");
    E2->GetYaxis()->SetRangeUser(0, 2500);
    E2->Draw("SAME L");
    E3->GetXaxis()->SetTitle("Energy(GeV)");
    E3->GetYaxis()->SetTitle("Counts");
    E3->GetYaxis()->SetRangeUser(0,2500);
    E3->Draw("SAME L");
    TLegend *legend3 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend3->AddEntry(E1, "0n0n", "l");
    legend3->AddEntry(E2, "0nXn", "l");
    legend3->AddEntry(E3, "XnXn", "l");
    // legend->AddEntry(EtaT, "Total", "l");
    legend3->Draw();

    c3->Draw();
    c3->SaveAs("Energy_Boostedd_2.76.png");

}

void DrawPseudoRapidity(TH1D* Eta1,TH1D* Eta2,TH1D* Eta3)
{
    TCanvas *c2 = new TCanvas("c2", "Pseudo Rapidity", 1200, 800);

    
    Eta1->SetLineColor(kRed);
    Eta2->SetLineColor(kBlue);
    Eta3->SetLineColor(kGreen);

    // Draw histograms
    Eta1->GetXaxis()->SetTitle("Pseudo Rapidity");
    Eta1->GetYaxis()->SetTitle("Counts");
    Eta1->GetYaxis()->SetRangeUser(0, 250);
    Eta1->Draw("L");
    Eta2->GetXaxis()->SetTitle("Pseudo Rapidity");
    Eta2->GetYaxis()->SetTitle("Counts");
    Eta2->GetYaxis()->SetRangeUser(0, 250);
    Eta2->Draw("SAME L");
    Eta3->GetXaxis()->SetTitle("Pseudo Rapidity");
    Eta3->GetYaxis()->SetTitle("Counts");
    Eta3->GetYaxis()->SetRangeUser(0, 250);
    Eta3->Draw("SAME L");
    TLegend *legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend2->AddEntry(Eta1, "0n0n", "l");
    legend2->AddEntry(Eta2, "0nXn", "l");
    legend2->AddEntry(Eta3, "XnXn", "l");
    // legend->AddEntry(EtaT, "Total", "l");
    legend2->Draw();

    c2->Draw();
    c2->SaveAs("Pseudo_rapidity_Boosted_2.76.png");
}
