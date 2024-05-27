
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
#include "TRandom.h"

#include "NeutronGenerator.h"
#endif

void DrawXSection(TH1D* XS1, TH1D* XS2, TH1D* XS3, TH1D* XST);

void DrawEnergy(TH1D* E1, TH1D* E2, TH1D* E3);

void DrawPseudo(TH1D* Eta1,TH1D* Eta2,TH1D* Eta3);

void DrawNY(TH1D* NeutronY2, TH1D* NeutronY3);

void DrawNEta(TH1D* NeuttonEta2, TH1D* NeutronEta3);

void DrawNE(TH1D* NeutronE2, TH1D* NeutronE3);

void DrawTotNE(TH1D *TotE);


void Afterburner()
{

    NeutronGenerator *gen = new NeutronGenerator();

    gen->SetStoreQA();
    gen->SetStoreGeneratorFunctions();
    gen->SetHadronicInteractionModel(NeutronGenerator::kGlauber);
    gen->Initialize();
    gen->SetRunMode(NeutronGenerator::kInterface);
    gen->ReadENDF(kTRUE);
    gen->LoadENDF("hENDF.root");
    gen->Setup();


    TFile *file = new TFile("example_Pb.root","READ");

    TTree *tree = dynamic_cast<TTree*>(file->Get("tree"));

    Long64_t nE = tree->GetEntries();

    TLorentzVector* lv = new TLorentzVector;

    TBranch *bc = tree->GetBranch("vm");

    bc->SetAddress(&lv);

    std::vector<Double_t> onon; //srotes rapiity
    std::vector<Double_t> onXn; // stores rapidity
    std::vector<Double_t> XnXn; // stores rapidity

    std::vector<Double_t> ononEta;//stores pseudo rapidity
    std::vector<Double_t> onXnEta;//stores pseudo rapidity
    std::vector<Double_t> XnXnEta;//stores pseudo rapidity

    std::vector<Int_t> on; //stores index
    std::vector<Int_t> ln; // stores index
    std::vector<Int_t> Xn; //stores index

    std::vector<Double_t> ononE;
    std::vector<Double_t> onXnE;
    std::vector<Double_t> XnXnE;

    std::vector<Double_t> onNeutronE;
    std::vector<Double_t> onNeutronEta;
    std::vector<Double_t> onNeutronY;
    
    std::vector<Double_t> lnNeutronE;
    std::vector<Double_t> lnNeutronEta;
    std::vector<Double_t> lnNeutronY;

    std::vector<Double_t> XnNeutronE;
    std::vector<Double_t> XnNeutronEta;
    std::vector<Double_t> XnNeutronY;

    std::vector<Double_t> NeutronTotE;

    std::vector<Double_t> onEmptyE;
    
    

    std::vector<Double_t> onEmptyEta;
    
    

    std::vector<Double_t> onEmptyY;
    
    

    std::vector<Double_t> VmBoostE;


    // Double_t Mv = 3.09;
    Double_t beamGamma = 1471;

    for(Long64_t index = 0; index < nE; index++)
    {
        bc->GetEntry(index);

        Double_t y = lv->Rapidity();

        Double_t eta= lv->Eta();

        Double_t k = 0.5*lv->M()*TMath::Exp(TMath::Abs(y));

        std::vector<Int_t> nNumbers = gen->runSartreNoon(k); // Gives number of neutrons in both beams

        if(nNumbers[0]==0 && nNumbers[1] == 0) 
        {
            onon.push_back(y);
            ononEta.push_back(eta);
            on.push_back(index);
            ononE.push_back(lv->Energy());
            // gen->createSartreNeutrons(nNumbers[0],nNumbers[1],onNeutronE,onNeutronEta,onNeutronY);

            // NeutronTotE.push_back(0.0);
            
        }
        else if(nNumbers[0]==0 || nNumbers[1]==0 )
        {
            std::vector<Double_t> lnEmptyY;
            std::vector<Double_t> lnEmptyEta;
            std::vector<Double_t> lnEmptyE;

            onXn.push_back(y);
            onXnEta.push_back(eta);
            ln.push_back(index);
            onXnE.push_back(lv->Energy());
            
            gen->createSartreNeutrons(nNumbers[0],nNumbers[1],lnNeutronE,lnNeutronEta,lnNeutronY);
            
            if(gRandom->Rndm()<0.5)gen->createSartreNeutrons(nNumbers[0],nNumbers[1],lnEmptyE,lnEmptyEta,lnEmptyY);
            else gen->createSartreNeutrons(nNumbers[1],nNumbers[0],lnEmptyE,lnEmptyEta,lnEmptyY);

            Double_t TNE = 0;
            for(Int_t i=0;i<lnEmptyE.size();i++)
            {
                TNE +=lnEmptyE[i];
            }
            NeutronTotE.push_back(TNE);
        }
        else if(nNumbers[0]!=0 && nNumbers[1]!=0) 
        {
            std::vector<Double_t> XnEmptyY;
            std::vector<Double_t> XnEmptyE;
            std::vector<Double_t> XnEmptyEta;

            XnXn.push_back(y);
            XnXnEta.push_back(eta);
            Xn.push_back(index);
            XnXnE.push_back(lv->Energy());
            gen->createSartreNeutrons(nNumbers[0],nNumbers[1],XnNeutronE,XnNeutronEta,XnNeutronY);
            if(gRandom->Rndm()<0.5)gen->createSartreNeutrons(nNumbers[0],nNumbers[1],XnEmptyE,XnEmptyEta,XnEmptyY);
            else gen->createSartreNeutrons(nNumbers[1],nNumbers[0],XnEmptyE,XnEmptyEta,XnEmptyY);

            Double_t TNE = 0;
            for(Int_t i=0;i<XnEmptyE.size();i++)
            {
                TNE +=XnEmptyE[i];
            }
            NeutronTotE.push_back(TNE);
        }
        else std::cout<<"Error in the Event"<<std::endl;

        // if(gRandom->Rndm()<0.5)lv->Boost(0,0,1*TMath::Sqrt(1.0-1.0/beamGamma/beamGamma));
        // // else lv->Boost(0,0,1*TMath::Sqrt(1.0-1.0/beamGamma/beamGamma));

        // VmBoostE.push_back(lv->Energy());

    }
    delete lv;
    delete tree;
    file->Close();
    delete file;

    std::cout<<"0n0n number of events: "<<onon.size()<<std::endl;
    std::cout<<"0nXn number of events: "<<onXn.size()<<std::endl;
    std::cout<<"XnXn number of events: "<<XnXn.size()<<std::endl;

//--------------------------------------------------------------------------------------------------------
    // Double_t minRange = -6;
    // Double_t maxRange = 6;
    // // Create histograms for each rapidity list
    // TH1D *hist1 = new TH1D("hist1", "vm Rapidity Distribution ", 1000, minRange,maxRange);
    // TH1D *hist2 = new TH1D("hist2", "vm Rapidity Distribution ", 1000, minRange, maxRange);
    // TH1D *hist3 = new TH1D("hist3", "vm Rapidity Distribution ", 1000, minRange,maxRange);

    // TH1D *Eta1 = new TH1D("Eta1", "vm PseudoRapidity Distribution ", 1000,  -15,15);
    // TH1D *Eta2 = new TH1D("Eta2", "vm PseudoRapidity Distribution ", 1000, -15,15);
    // TH1D *Eta3 = new TH1D("Eta3", "vm PseudoRapidity Distribution ", 1000, -15,15);

    // TH1D *E1 = new TH1D("E1", "vm Energy ", 1000, 0 ,300);
    // TH1D *E2 = new TH1D("E2", "vm Energy ", 1000, 0, 300);
    // TH1D *E3 = new TH1D("E3", "vm Energy ", 1000, 0, 300);

    // for (Int_t i; i<onon.size();i++) {
        
    //     // bc->GetEntry(on[i]);

    //     hist1->Fill(onon[i]);
    //     Eta1->Fill(ononEta[i]);
    //     E1->Fill(ononE[i]);

    // }
    // for (Int_t i; i<onXn.size();i++) {
        
    //     // bc->GetEntry(ln[i]);
        
    //     hist2->Fill(onXn[i]);
    //     Eta2->Fill(onXnEta[i]);
    //     E2->Fill(onXnE[i]);
    // }
    
    // for (Int_t i; i<XnXn.size();i++) {

    //     // bc->GetEntry(Xn[i]);

    //     hist3->Fill(XnXn[i]);
    //     Eta3->Fill(XnXnEta[i]);
    //     E3->Fill(XnXnE[i]);
    // }

//     // Create histograms for cross-section distribution
//     TH1D *XS1 = new TH1D("XS1", "Cros-Section", 1000, minRange, maxRange);
//     TH1D *XS2 = new TH1D("XS2", "Cros-Section", 1000, minRange,maxRange);
//     TH1D *XS3 = new TH1D("XS3", "Cros-Section", 1000, minRange,maxRange);
//     TH1D *XST = new TH1D("XST", "Cros-Section", 1000, minRange,maxRange);

//     for(Int_t entry = 0; entry<1000; entry++)
//     {
//         XS1->SetBinContent(entry,(hist1->GetBinContent(entry)) * 523 * 1000 /100000);
//         XS2->SetBinContent(entry,(hist2->GetBinContent(entry)) * 523 * 1000 /100000);
//         XS3->SetBinContent(entry,(hist3->GetBinContent(entry)) * 523 * 1000 /100000);
//         XST->SetBinContent(entry,XS1->GetBinContent(entry)+XS2->GetBinContent(entry)+XS3->GetBinContent(entry));
//     }
    
//     /**
//      * Plotting and Saving as png
//     */
//     // DrawXSection(XS1,XS2,XS3,XST);
//     // DrawEnergy(E1,E2,E3);
    // DrawPseudo(Eta1,Eta2,Eta3);

//     /**
//      * Plotting and Saving as png
//     */

//----------------------------------------------------------------------------------------------------------

//     /**
//      * 
//      * Neutron Data in the following histograms
//      * 
//     */

    //    TH1D *NeutronE1 = new TH1D("NeutronE1","Neutron energy",1000,0,5);
   TH1D *NeutronE2 = new TH1D("NeutronE2","Neutron energy",10000,0,5000);
   TH1D *NeutronE3 = new TH1D("NeutronE3","Neutron energy",10000,0,5000);
   
    //    TH1D *NeutronEta1 = new TH1D("NeutronEta1", "Neutron Pseudo Rapidity",1000,-10,10);
   TH1D *NeutronEta2 = new TH1D("NeutronEta2", "Neutron Pseudo Rapidity",1000,-25,25);
   TH1D *NeutronEta3 = new TH1D("NeutronEta3", "Neutron Pseudo Rapidity",1000,-25,25);

    //    TH1D *NeutronY1 = new TH1D("NeutronY1","Neutron Rapidity",1000,-6,6);
   TH1D *NeutronY2 = new TH1D("NeutronY2","Neutron Rapidity",1000,-25,25);
   TH1D *NeutronY3 = new TH1D("NeutronY3","Neutron Rapidity",1000,-25,25);


   for(Int_t i=0; i<lnNeutronE.size(); i++)
   {
    NeutronE2->Fill(lnNeutronE[i]);
    NeutronEta2->Fill(lnNeutronEta[i]);
    NeutronY2->Fill(lnNeutronY[i]);
   }

   for(Int_t i=0; i<XnNeutronE.size(); i++)
   {
    NeutronE3->Fill(XnNeutronE[i]);
    NeutronEta3->Fill(XnNeutronEta[i]);
    NeutronY3->Fill(XnNeutronY[i]);
   }

    /**
     * Plotting the graphs and saving as png file
    */
//    DrawNY(NeutronY2,NeutronY3);
//    DrawNEta(NeutronEta2,NeutronEta3);
   DrawNE(NeutronE2,NeutronE3);
   /**
    *  End of Plotting
    */
   
//-------------------------------------------------------------------------------------------------------------------------------

   std::cout<<"Number of Neutrons in onon: "<<onNeutronE.size()<<std::endl;
   std::cout<<"Number of Neutrons in onxn: "<<lnNeutronE.size()<<std::endl;
   std::cout<<"Number of Neutrons in xnxn: "<<XnNeutronE.size()<<std::endl;



  TH1D *TotE = new TH1D("TotE","Total Neutron Energy in each Event",1000,0,30000);
//   TH1D *TotEL = new TH1D("ThotEL","Log scale of Neutron Energy per Event", 1000,0,10000);

  for(Int_t i=0; i<NeutronTotE.size();i++)
  {
    TotE->Fill(NeutronTotE[i]);
  }
//   for(Int_t i=0; i<100;i++)
//   {
//     TotEL->SetBinContent(i,TMath::Log10(TotE->GetBinContent(i)));
//   }

  DrawTotNE(TotE);
// // //-------------------------------------------------------------------------------------------------------
    // TH1D *VmB = new TH1D("VmB","Boosted Vectormeson energy",100,0,5000);
    // TH1D *VmBL = new TH1D("VmBLog","Log scale",100,0,5000);
    // for(Int_t i=0; i<VmBoostE.size();i++)
    // {
    //     VmB->Fill(VmBoostE[i]);
    // }
    // for(Int_t i=0; i<1000;i++)
    // {
    //     VmBL->SetBinContent(i,TMath::Log10(VmB->GetBinContent(i)));   
    // } 

    // DrawTotNE(VmBL);

}


/**
 * 
 * 
 * The Following functions are for ploting the data
 * 
 * 
*/

void DrawTotNE(TH1D *TotE)
{
    TCanvas *c7 = new TCanvas("c7","TOtal Neutron Energy per Event side 1",1200,800);

    TotE->GetXaxis()->SetTitle("Energy(GeV)");
    TotE->SetAxisRange(0,5000,"X");
    TotE->GetYaxis()->SetTitle("Counts");
    TotE->Draw("L");
    c7->SetLogy();
    c7->SaveAs("Total_Neutron_E_2.76_side2_ss.png");
}

void DrawNY( TH1D* NeutronY2, TH1D* NeutronY3)
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

void DrawNEta( TH1D* NeutronEta2, TH1D* NeutronEta3)
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


void DrawNE( TH1D* NeutronE2, TH1D* NeutronE3)
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

void DrawPseudo(TH1D* Eta1,TH1D* Eta2,TH1D* Eta3)
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