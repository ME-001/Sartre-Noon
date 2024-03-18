/*
 * This program is for collecting and storing the data event by event
 * 
 * The other program "extract.C" collects the data and stores in histogram. 
 * 
*/

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

#include "NeutronGenerator.h"
#endif

void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07,int columns=2)
{
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  currentLegend->SetNColumns(columns);
  return;
}

void runSartre()
{   
    Double_t  VMR = 0;
    Double_t VMmass = 3.09; //  j/psi mass
    Double_t photonK = 0;
    Int_t numberOfbins = 100;
    // Double_t cs;
    // Double_t ne;

    const Int_t cs = 524;//cross section 524 nb
    //const Int_t ne = 100000;//Number of events 1M

    // std::cout<<"Enter vector meson mass "<<endl;
    // std::cin>>VMmass;

    // std::cout<<"Enter number of events "<<endl;
    // std::cin>>ne;
    
    // std::cout<<"Enter cross section in nb "<<endl;
    // std::cin>>cs;



    TFile *file = new TFile("example_Pb.root","READ");

    TTree *tree = dynamic_cast<TTree*>(file->Get("tree"));

    TLorentzVector* lorentzVector = new TLorentzVector();
    TLorentzVector* lovec = new TLorentzVector();

    TBranch *branch = tree->GetBranch("vm");
    //TBranch *bch = tree->GetBranch("gamma");

    branch->SetAddress(&lorentzVector);
    //bch->SetAddress(&lovec);

    Long64_t nEntries = tree->GetEntries();
    Long64_t ne = nEntries;
    
    std::vector<double> rapidityValues;
    std::vector<double> energyValues;
    std::vector<double> phk;

    for (Long64_t entry = 0; entry < nEntries; entry++) {
        branch->GetEntry(entry);
        //bch->GetEntry(entry);

        Double_t rapidity = lorentzVector->Rapidity();
        rapidityValues.push_back(rapidity);

        Double_t energy = lorentzVector->Energy();
        energyValues.push_back(energy);

        Double_t k = 0.5*VMmass*TMath::Exp(TMath::Abs(rapidity));

        phk.push_back(k);
        
        //PhotonK->Fill(energy);
        //Rapidity->Fill(rapidity);

    }
    delete lorentzVector;
    delete tree;
    file->Close();
    delete file;

    Double_t maxRapidity = *std::max_element(rapidityValues.begin(), rapidityValues.end());
    Double_t minRapidity = *std::min_element(rapidityValues.begin(), rapidityValues.end());

    Double_t minEnergy = *std::min_element(energyValues.begin(), energyValues.end());
    Double_t maxEnergy = *std::max_element(energyValues.begin(), energyValues.end());

    #if defined(__CINT__)
        gROOT->LoadMacro("NeutronGenerator.cxx+g");
    #endif
    
    TClonesArray *fParticles = new TClonesArray("TParticle", energyValues.size());
    TTree *Tree = new TTree("Tree", "Tree");
    Tree ->Branch("fParticles", &fParticles);

    TClonesArray *fNGparticles = NULL;


    NeutronGenerator *gen = new NeutronGenerator();

    gen->SetStoreQA();
    gen->SetStoreGeneratorFunctions();
    gen->SetHadronicInteractionModel(NeutronGenerator::kGlauber);
    gen->Initialize();
    gen->SetRunMode(NeutronGenerator::kInterface);
    gen->ReadENDF(kTRUE);
    gen->LoadENDF("hENDF.root");
    gen->Setup();


    

    cout<<"Running production"<<endl; 

    UInt_t nEvents = energyValues.size();

    

    for(Int_t iEvent = 0; iEvent < nEvents; ++iEvent)
    {   
        photonK = phk[iEvent];

        VMR = rapidityValues[iEvent];

        //NeutronGenerator::hRapidityVM->Fill(VMR);
        //NeutronGenerator::hMassVM->Fill(VMmass);

        gen->GenerateEvent(photonK);
        Int_t nTotalPart = 0; 
        fNGparticles = gen->ImportParticles();
        for(Int_t i = 0; i<fNGparticles->GetEntriesFast(); i++)
        {
            TParticle *part(dynamic_cast<TParticle*>(fNGparticles->At(i)));
            new((*fParticles)[nTotalPart++]) TParticle(*part);
        }
        gen->FinishEvent();
    
        Tree->Fill();
        fParticles->Clear("C");

        gen->FinishEvent();

        
    }

    gen->FinishProduction();
    cout<<"Finish production"<<endl;
    TFile *fOutputFile = new TFile("output1.root","RECREATE");
    Tree->Write();
    fOutputFile->Close();


    //std::unordered_map<Double_t, Int_t> rapidityCounts;
    //std::unordered_map<Double_t, Int_t> energyCounts;

    //for(const auto &element : rapidityValues) rapidityCounts[element]+=1;
    //for(const auto &element :energyValues) energyCounts[element]+=1;

    TH1D *PhotonK = new TH1D("PhotonK","PhotonK",numberOfbins,minEnergy,maxEnergy);
    TH1D *vmEnergy = new TH1D("vm Energy","Vm Energy",numberOfbins,minEnergy,maxEnergy);
    TH1D *Rapidity = new TH1D("Rapidity","Rapidity",numberOfbins,minRapidity,maxRapidity);

    TH1D *massHist = new TH1D("massHist","massHist",1000,3,3.2);
    Double_t mean = VMmass;
    Double_t sd = 0.01;

    for (Long64_t entry = 0; entry < nEntries; entry++) 
    {
        
        //Double_t energy = energyValues[entry];
        Double_t rapidity = rapidityValues[entry];
        // Double_t energy = 0.5*3.09*TMath::Exp(TMath::Abs(rapidity));
        Double_t k = phk[entry];
        Double_t vmE = energyValues[entry];

        Double_t massValue = gRandom->Gaus(mean,sd);

        PhotonK->Fill(k);
        Rapidity->Fill(rapidity);
        vmEnergy->Fill(vmE);

        massHist->Fill(massValue);

    }

    TFile *file2 = new TFile("exData.root","RECREATE");

    PhotonK->Write();
    Rapidity->Write();
    vmEnergy->Write();
    massHist->Write();

    TH1D *xSection = new TH1D("xSection","xSection",numberOfbins,minRapidity,maxRapidity);

    Double_t xSectionValue = 0;

    //TH1D *hInputMass = massHist;
    //TH1D *hInputRapidity = xSection;

    TString breakups[] = {"All","0n0n","Xn0n","XnXn"};
    TH1D *hRapidityBreakup[4];

    Double_t prob[4];

    for(Long64_t index = 1; index <=numberOfbins; index++)
    {
        Double_t binValue = Rapidity->GetBinContent(index);
        Double_t rpdt = Rapidity->GetBinCenter(index);

        //Double_t binWidth = Rapidity->GetBinWidth(index);

        xSectionValue = binValue * cs / ne ;//binWidth ;

        xSection->SetBinContent(index, xSectionValue);

        for(Int_t k=0; k<4 ;++k)
        {
            prob[k] = 0;

            for(Int_t iEvent = 0; iEvent<100  ;++ne)
            {
                Double_t vmm = massHist->GetRandom();
                Double_t pk = 0.5*vmm*TMath::Exp(TMath::Abs(rpdt));
                if(k == 0)
                {
  	                prob[k] += 1.0;
                    
                }
                if(k == 1)
                {   
                    
  	                prob[k] += gen->GetBreakupProbability(pk, 0, 0); 
                    
  	            }
                if(k == 2)
                {
  	                prob[k] += gen->GetBreakupProbability(pk, -1, 0); 
                    
  	            }
                if(k == 3)
                {
  	                prob[k] += gen->GetBreakupProbability(pk, -1, -1); 
                    
  	            }
                

            }
            prob[k] = prob[k]/numberOfbins;
            hRapidityBreakup[k]->SetBinContent(index, xSectionValue*prob[k]);

        }

        
    }

    xSection->Write();
    file2->Close();
    delete file2;
    

    TCanvas* d1 = new TCanvas("d1","d1");
        Rapidity->Draw();
    TCanvas* d2 = new TCanvas("d2","d2");
        xSection->Draw();
    TCanvas* d3 = new TCanvas("d3","d3");
        PhotonK->Draw();
    TCanvas* d4 = new TCanvas("d4","d4");
        vmEnergy->Draw();

    d1->SaveAs("Rapidity.png");
    d2->SaveAs("xSection.png");
    d3->SaveAs("PhotonK.png");
    d4->SaveAs("vmEnergy.png");

    

    delete d1;
    delete d2;
    delete d3;
    delete d4;

    

    TCanvas *c1 = new TCanvas("c1","c1",0,0,800,800);
    TCanvas *c2 = new TCanvas("c2","c2",0,0,800,800);
    TCanvas *c3 = new TCanvas("c3","c3",0,0,800,800);
 
    c1->cd();  
    hRapidityBreakup[0]->GetYaxis()->SetRangeUser(0,40);
    hRapidityBreakup[0]->DrawCopy();

    for(Int_t k=1; k<4; k++) hRapidityBreakup[k]->DrawCopy("same");

    gPad->SetGridy();gPad->SetGridx();
    c2->cd();

    for(Int_t k=1; k<4; k++)
    {
        hRapidityBreakup[k]->Divide(hRapidityBreakup[0]);
        hRapidityBreakup[k]->GetYaxis()->SetRangeUser(0.001,0.9);
        if(k == 1)hRapidityBreakup[k]->DrawCopy();
        hRapidityBreakup[k]->DrawCopy("same");
    }
    gPad->SetGridy();gPad->SetGridx();
    c3->cd();
    gPad->SetLogy();

    for(Int_t k=0; k<4; k++)
    { 
        hRapidityBreakup[k]->GetYaxis()->SetRangeUser(0.01,1);
        if(k == 0)hRapidityBreakup[k]->DrawCopy();
        hRapidityBreakup[k]->DrawCopy("same");
    }
  
    TLegend *myLegend1 = new TLegend(0.42,0.29,0.69,0.48);
    myLegendSetUp(myLegend1,0.04,1);
    myLegend1->AddEntry(hRapidityBreakup[0],"All","l");
    myLegend1->AddEntry(hRapidityBreakup[1],"0n0n","l");
    myLegend1->AddEntry(hRapidityBreakup[2],"0nXn","l");
    myLegend1->AddEntry(hRapidityBreakup[3],"XnXn","l");
    c1->cd(); myLegend1->Draw();
  
    TLatex *noontitle = new TLatex(0.45,0.83,"Sartre + #bf{n_{O}^{O}n}");//"Hot-spot model + #bf{n_{O}^{O}n}");
    noontitle->SetNDC();
    noontitle->SetTextFont(42);
    noontitle->SetTextSize(0.04);
    noontitle->Draw();
  
    TLegend *myLegend2 = new TLegend(0.42,0.29,0.69,0.48);
    myLegendSetUp(myLegend2,0.04,1);
    myLegend2->AddEntry(hRapidityBreakup[1],"0n0n/All","l");
    myLegend2->AddEntry(hRapidityBreakup[2],"0nXn/All","l");
    myLegend2->AddEntry(hRapidityBreakup[3],"XnXn/All","l");
    c2->cd(); myLegend2->Draw();
    c3->cd(); myLegend2->Draw();



    
    



}
