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
#endif


TH1D* exRapidity(TTree* tree, TBranch* branch);
TH1D* exEnergy(TTree* tree, TBranch* branch);
TH1D* exSection(TH1D* histogram);

void extract();
void runSartre();
void computeModelBreakups();

void run(){

    extract();
    runSartre();
    computeModelBreakups();

}

void runSartre(){

#if defined(__CINT__)
  gROOT->LoadMacro("NeutronGenerator.cxx+g");
#endif

  TClonesArray *fParticles = new TClonesArray("TParticle", 200);
  TTree *fEventTree = new TTree("fEventTree", "fEventTree");
  fEventTree ->Branch("fParticles", &fParticles);
  TClonesArray *fSLparticles = new TClonesArray("TParticle", 200);

  TMap *fSLparams = new TMap();

    gROOT->LoadMacro("NeutronGenerator.cxx+g");


    TClonesArray *fNGparticles = NULL;

    //NeutronGenerator *gen = new NeutronGenerator();

  NeutronGenerator *gen = new NeutronGenerator();
  gen->SetStoreQA();
  gen->SetStoreGeneratorFunctions();
  gen->SetHadronicInteractionModel(NeutronGenerator::kHardSphere); 
  gen->Initialize(); 
  gen->SetRunMode(NeutronGenerator::kInterface);  
  gen->ReadENDF(kFALSE);
  gen->Setup();
  //gen->Run(10000);
  
    // TDatabasePDG pdgData;
    Double_t VMrapidity = 0;
    Double_t VMmass = 3.09;// jebspi mass
    Double_t photonK = 0;
    // Double_t fRapMin = -4.0;
    // Double_t fRapMax = 4.0;

    UInt_t nEvents = 10000;
    cout<<"Running production"<<endl; 

    // TFile *inputFile = new TFile("Rapidity_Pb.root","READ");
    // TFile *inputFile2 = new TFile("Energy_Pb.root","READ");
    
    TFile *inputfile = new TFile("extract.root","READ");

    TH1D *Photonk = (TH1D*)inputfile->Get("PhotonK");
    //TH1D *rapidity = (TH1D*)inputfile->Get("Rapidity");

     //Int_t nEvents = Photonk->GetNbinsX()+1;

    for(Int_t iEvent = 1 ; iEvent<=nEvents;iEvent++)
    { 

      photonK = Photonk->GetBinContent(iEvent);
      //std::cout<<photonK<<std::endl;

      VMrapidity = rapidity->GetBinContent(iEvent); //TMath::Abs(TMath::Log(2*photonK/3.09)); //  change 3.09 to mass variable

      gen->GenerateEvent(photonK);

      gen->FinishEvent();


    }
    
    gen->FinishProduction();
    cout<<"Finish production"<<endl;
 
}

void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07,int columns=2){
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

void computeModelBreakups(){

  #if defined(__CINT__)
  gROOT->LoadMacro("NeutronGenerator.cxx+g");
#endif

  NeutronGenerator *gen = new NeutronGenerator();
  
  gen->SetStoreQA();
  gen->SetStoreGeneratorFunctions();
  gen->SetHadronicInteractionModel(NeutronGenerator::kHardSphere); 
  gen->Initialize(); 
  gen->SetRunMode(NeutronGenerator::kInterface);  
  gen->ReadENDF(kFALSE);
  gen->Setup();
  

//   TFile *inputFile = new TFile("XSection.root","READ");
//   TFile *inputFile2 = new TFile("Energy_Pb.root","READ");
  // TH1D *hInputMass = (TH1D*)inputFile->Get("massHist");
  TFile *inputfile = new TFile("extract.root","READ");
  TH1D *hInputRapidity = (TH1D*)inputfile->Get("XSection");
  TH1D *Photonk = (TH1D*)inputfile->Get("PhotonK");

  TString breakups[] = {"All","0n0n","Xn0n","XnXn"};
  TH1D *hRapidityBreakup[4];
  
  Double_t VMrapidity = 0;
  Double_t VMmass = 3.09;
  //Double_t photonK_Low = 0;
  //Double_t photonK_High = 0;
  Double_t probLow[4];
  Double_t probHigh[4];
  Double_t photonK = 0;
  Double_t prob[4];
  
  //for(Int_t j=1; j<=41; j++)hInputRapidity->SetBinContent(j,hInputRapidity->GetBinContent(j)*gen->GetTotalFlux(0.5*3.09*TMath::Exp(hInputRapidity->GetBinCenter(j))));
  
  hInputRapidity->SetLineWidth(2);
  hInputRapidity->SetLineColor(kBlack);
  hInputRapidity->SetStats(kFALSE);
  hInputRapidity->SetTitle("");
  hInputRapidity->GetXaxis()->SetTitle("y");
  hInputRapidity->GetYaxis()->SetTitle("d#sigma/dy[nb]");
  hInputRapidity->GetYaxis()->SetTitleOffset(1.5);
  
  for(Int_t k=0; k<4; k++){
    TString breakupName = breakups[k].Data();
    hRapidityBreakup[k] = (TH1D*)hInputRapidity->Clone(breakupName.Data());
    hRapidityBreakup[k]->SetLineWidth(2);
    hRapidityBreakup[k]->SetLineColor(1+k);
    hRapidityBreakup[k]->SetStats(kFALSE);
    hRapidityBreakup[k]->GetXaxis()->SetTitle("y");
    hRapidityBreakup[k]->GetYaxis()->SetTitle("#sigma [nb]");
    hRapidityBreakup[k]->GetYaxis()->SetTitleOffset(1.5);
    }
  Int_t nBinsInput = hInputRapidity->GetNbinsX()+1;
  for(Int_t j=1; j<=nBinsInput/2; j++){
    VMrapidity = hInputRapidity->GetBinCenter(j);
    photonK = Photonk->GetBinContent(j);
    for(Int_t k=0; k<4; k++){
      //probLow[k] = 0;
      //probHigh[k] = 0;
      prob[k] = 0;
      for(Int_t iEvent = 0; iEvent<10000; iEvent++){
        //VMmass = hInputMass->GetRandom();
        //photonK_Low = 0.5*VMmass*TMath::Exp(VMrapidity);
        //photonK_High = 0.5*VMmass*TMath::Exp(-1*VMrapidity);
        if(k == 0){
  	      //probLow[k] += 1.0;
          //probHigh[k] += 1.0;
          prob[k] += 1.0;

        }
        if(k == 1){
  	      //probLow[k] += gen->GetBreakupProbability(photonK_Low, 0, 0); 
          //probHigh[k] += gen->GetBreakupProbability(photonK_High, 0, 0);
          prob[k] += gen->GetBreakupProbability(photonK,0,0);

  	    }
        if(k == 2){
  	      //probLow[k] += gen->GetBreakupProbability(photonK_Low, -1, 0); 
          //probHigh[k] += gen->GetBreakupProbability(photonK_High, -1, 0);
          prob[k] += gen->GetBreakupProbability(photonK,-1,0);
  	    }
        if(k == 3){
  	      //probLow[k] += gen->GetBreakupProbability(photonK_Low, -1, -1); 
          //probHigh[k] += gen->GetBreakupProbability(photonK_High, -1, -1);
          prob[k] += gen->GetBreakupProbability(photonK,-1,-1);
  	    }
      }
      //probLow[k] /= 1000;
      //probHigh[k] /= 1000;
      prob[k] /= 10000; 
      //hRapidityBreakup[k]->SetBinContent(j,hInputRapidity->GetBinContent(j)*probLow[k]+hInputRapidity->GetBinContent(nBinsInput-j)*probHigh[k]);
      //hRapidityBreakup[k]->SetBinContent(nBinsInput-j,hRapidityBreakup[k]->GetBinContent(j));
      hRapidityBreakup[k]->SetBinContent(j,hInputRapidity->GetBinContent(j)*prob[k]);
      hRapidityBreakup[k]->SetBinContent(nBinsInput-j,hRapidityBreakup[k]->GetBinContent(j));
    }
  }

  TCanvas *c1 = new TCanvas("c1","c1",0,0,800,800);
  TCanvas *c2 = new TCanvas("c2","c2",0,0,800,800);
  TCanvas *c3 = new TCanvas("c3","c3",0,0,800,800);
 
  c1->cd();  
  hRapidityBreakup[0]->GetYaxis()->SetRangeUser(0,15000);
  hRapidityBreakup[0]->DrawCopy();
  for(Int_t k=1; k<4; k++) hRapidityBreakup[k]->DrawCopy("same");
  gPad->SetGridy();gPad->SetGridx();
  c2->cd();
  for(Int_t k=1; k<4; k++){ 
    hRapidityBreakup[k]->Divide(hRapidityBreakup[0]);
    hRapidityBreakup[k]->GetYaxis()->SetRangeUser(0.001,0.9);
    if(k == 1)hRapidityBreakup[k]->DrawCopy();
    hRapidityBreakup[k]->DrawCopy("same");
    }
   gPad->SetGridy();gPad->SetGridx();
   c3->cd();
   gPad->SetLogy();
   for(Int_t k=0; k<4; k++){ 
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
  
  TLatex *noontitle = new TLatex(0.45,0.83,"Hot-spot model + #bf{n_{O}^{O}n}");
  noontitle->SetNDC();
  noontitle->SetTextFont(42);
  noontitle->SetTextSize(0.04);
  //noontitle->Draw();
  
  TLegend *myLegend2 = new TLegend(0.42,0.29,0.69,0.48);
  myLegendSetUp(myLegend2,0.04,1);
  myLegend2->AddEntry(hRapidityBreakup[1],"0n0n/All","l");
  myLegend2->AddEntry(hRapidityBreakup[2],"0nXn/All","l");
  myLegend2->AddEntry(hRapidityBreakup[3],"XnXn/All","l");
  c2->cd(); myLegend2->Draw();
  c3->cd(); myLegend2->Draw();

}


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
    exEnergy(tree, tree->GetBranch("vm"))->Write();
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