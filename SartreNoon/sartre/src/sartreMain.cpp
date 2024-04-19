//==============================================================================
//  sartreMain.cpp
//
//  Copyright (C) 2010-2019 Tobias Toll and Thomas Ullrich
//
//  This file is part of Sartre.
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation.
//  This program is distributed in the hope that it will be useful,
//  but without any warranty; without even the implied warranty of
//  merchantability or fitness for a particular purpose. See the
//  GNU General Public License for more details.
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//  Author: Thomas Ullrich
//  Last update:
//  $Date: 2019-03-08 14:13:19 -0500 (Fri, 08 Mar 2019) $
//  $Author: ullrich $
//==============================================================================
//
//  Example main program. Use to get started.
//
//==============================================================================
#include <iostream>
#include <cmath>
#include "TTree.h"
#include "TFile.h"
#include "Sartre.h"
#include "DipoleModelParameters.h"
#include "TGenPhaseSpace.h"
#include "Settings.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

#define PR(x) cout << #x << " = " << (x) << endl;
//_________________________________________________
//_______________________________________________________

// #include <iostream>
// #include <fstream>
// #include <string>

// #include "TROOT.h"
// #include "TH1.h"
// #include "TH1D.h"
// #include "TH2D.h"
// #include "TString.h"
// #include "TGraph.h"
// #include "TCanvas.h"
// #include "TLegend.h"
// #ifdef __CLING__
// #include "NeutronGenerator.cxx+g"
// #endif


// #include <vector>
// #include "TFile.h"
// #include "TTree.h"
// #include "TLorentzVector.h"
// #include <cmath>

// #include "TFile.h"
// #include "TTree.h"
// #include "TLorentzVector.h"
// #include <cmath>
// #include "TMath.h"
// #include "TDatabasePDG.h"

// #include "NeutronGenerator.h"
// // #include "NeutronGenerator.cxx"





//
//_________________________________________________________________
using namespace std;

struct rootSartreEvent {
    double t;
    double Q2;
    double x;
    double s;
    double y;
    double W;
    double xpom;
    int    iEvent;
    int    pol;      // 0=transverse or 1=longitudinal
    int    dmode;    // 0=coherent, 1=Incoherent
};

rootSartreEvent myRootSartreEvent;
void randomlyReverseBeams(Event* );  // for UPC mode

int main(int argc, char *argv[])
{

    // #if defined(__CLINT__)
    //     gROOT->LoadMacro("NeutronGenerator.cxx+g");
    // #endif
    // NeutronGenerator *gen = new NeutronGenerator();

    // gen->SetStoreQA();
    // gen->SetStoreGeneratorFunctions();
    // gen->SetHadronicInteractionModel(NeutronGenerator::kGlauber);
    // gen->Initialize();
    // gen->SetRunMode(NeutronGenerator::kInterface);
    // gen->ReadENDF(kTRUE);
    // gen->LoadENDF("hENDF.root");
    // gen->Setup();


    //
    //  Check command line arguments
    //
    if (! (argc == 2 || argc == 3) ) {
        cout << "Usage:  " << argv[0] << " runcard [rootfile]" << endl;
        return 2;
    }
    string runcard = argv[1];
    
    //
    //  Create the generator and initialize it.
    //  Once initialized you cannot (should not) change
    //  the settings w/o re-initialing sartre.
    //
    Sartre sartre;
    bool ok = sartre.init(runcard);
    if (!ok) {
        cerr << "Initialization of sartre failed." << endl;
        return 1;
    }

    EventGeneratorSettings* settings = sartre.runSettings();

    //
    //  ROOT file
    //  Use the one from command line unless not given
    //  in which case the one in the runcard is used.
    //
    string rootfile;
    if (argc == 3) {
        rootfile = argv[2];
        settings->setRootfile(argv[2]);
    }
    else
        rootfile = settings->rootfile();
 
    //
    //  Print out all Sartre settings
    //
    settings->list();
    
    TFile *hfile = 0;
    if (rootfile.size()) {
        hfile  = new TFile(rootfile.c_str(),"RECREATE");
        cout << "ROOT file is '" <<  rootfile.c_str() << "'." << endl;
    }
    
    //
    //  Setup ROOT tree
    //
    TLorentzVector *eIn = new TLorentzVector;
    TLorentzVector *pIn = new TLorentzVector;
    TLorentzVector *vm = new TLorentzVector;
    TLorentzVector *eOut = new TLorentzVector;
    TLorentzVector *pOut = new TLorentzVector;
    TLorentzVector *gamma = new TLorentzVector;
    TLorentzVector *vmDaughter1 = new TLorentzVector;
    TLorentzVector *vmDaughter2 = new TLorentzVector;
    TLorentzVector *neutrons1 = new TLorentzVector;
    TLorentzVector *neutrons2 = new TLorentzVector;
    
    TClonesArray protons("TLorentzVector");
    TClonesArray neutrons("TLorentzVector");
    TClonesArray remnants("TLorentzVector");
    // TClonesArray beam1("TClonesArray");
    // TClonesArray beam2("TClonesArray");
    
    
    TTree tree("tree","sartre");
    tree.Branch("event", &myRootSartreEvent.t,
                "t/D:Q2/D:x/D:s/D:y/D:W/D:xpom/D:iEvent/I:pol/I:dmode/I");
    tree.Branch("eIn",  "TLorentzVector", &eIn, 32000, 0);
    tree.Branch("pIn",  "TLorentzVector", &pIn, 32000, 0);
    tree.Branch("vm",   "TLorentzVector", &vm, 32000, 0);
    tree.Branch("eOut", "TLorentzVector", &eOut, 32000, 0);
    tree.Branch("pOut", "TLorentzVector", &pOut, 32000, 0);
    tree.Branch("gamma","TLorentzVector", &gamma, 32000, 0);
    tree.Branch("vmDaughter1", "TLorentzVector", &vmDaughter1, 32000, 0);
    tree.Branch("vmDaughter2", "TLorentzVector", &vmDaughter2, 32000, 0);
    //_________________________________________________________________________________________
    // Modification section  Start
    //__________________________________________________________________________________
    /**
     * Creating two more branches
     * They will store array of particles data
    */
    
    tree.Branch("n1","TLorentzVector",&neutrons1,32000,0);
    tree.Branch("n2","TLorentzVector",&neutrons2,32000,0); 
    // tree.Branch("n1") // n1 must have 
    // tree.Branch("n2")
    // tree.Branch("neutronIndices") (1)

    //________________________________________________________________________________________________________________________
    //  Modification section END
    //__________________________________________________________________________________________________________________
    if(settings->enableNuclearBreakup()){
        tree.Branch("protons", &protons);
        tree.Branch("neutrons", &neutrons);
        tree.Branch("nuclearRemnants", &remnants);
    }
    
    //

    //  Prepare event generation
    //
    TGenPhaseSpace decay;  // for VM decays
    int daughterID = settings->userInt();
    double daughterMasses[2] = {0, 0};
    bool doPerformDecay = false;
    if (daughterID && settings->vectorMesonId() != 22) {
        doPerformDecay = true;
        daughterMasses[0] = settings->lookupPDG(daughterID)->Mass();
        daughterMasses[1] = settings->lookupPDG(-daughterID)->Mass();
        cout << "Will decay vector meson: ";
        cout << settings->lookupPDG(settings->vectorMesonId())->GetName();
        cout << " -> ";
        cout << settings->lookupPDG(daughterID)->GetName();
        cout << " ";
        cout << settings->lookupPDG(-daughterID)->GetName();
        cout << endl;
    }
    
    int nPrint;
    if (settings->timesToShow())
        nPrint = settings->numberOfEvents()/settings->timesToShow();
    else
        nPrint = 0;
    
    unsigned long maxEvents = settings->numberOfEvents();
    
    cout << "Generating " << maxEvents << " events." << endl << endl;
    
    //
    //  Event generation
    //
    for (unsigned long iEvent = 0; iEvent < maxEvents; iEvent++) {
        //
        //  Generate one event
        //
        Event *event = sartre.generateEvent();
        if (nPrint && (iEvent+1)%nPrint == 0 && iEvent != 0) {
            cout << "processed " << iEvent+1 << " events" << endl;
        }
        
        //
        //  If Sartre is run in UPC mode, half of the events needs to be
        //  rotated around and axis perpendicular to z:
        //
        if(settings->UPC() and settings->A()==settings->UPCA()){
            randomlyReverseBeams(event);
        }
        
        //
        //  Print out (here only for the first few events)
        //
        if (iEvent < 10) event->list();
        
        //
        //  Fill ROOT tree
        //
        myRootSartreEvent.iEvent = event->eventNumber;
        myRootSartreEvent.t = event->t;
        myRootSartreEvent.Q2 = event->Q2;
        myRootSartreEvent.x = event->x;
        myRootSartreEvent.y = event->y;
        myRootSartreEvent.s = event->s;
        myRootSartreEvent.W = event->W;
        myRootSartreEvent.xpom = event->xpom;
        myRootSartreEvent.pol = event->polarization == transverse ? 0 : 1;
        myRootSartreEvent.dmode = event->diffractiveMode == coherent ? 0 : 1;
        eIn     = &event->particles[0].p;
        pIn     = &event->particles[1].p;
        eOut    = &event->particles[2].p;
        pOut    = &event->particles[6].p;
        vm      = &event->particles[4].p;
        gamma   = &event->particles[3].p;

        // Event *e = sartre.GenerateNeutronEvent(event,vm);
        
        // for(int i = 0; i<2;i++)
        // {
        //     for(int j = 0; j<e->idx[i];j++)
        //     {   
        //         if(i==0)neutrons1 = &e->particles[7+j].p;
        //         else if(i==1)neutrons2 = &e->particles[7+e->idx[0]+j].p;         
        //     }
        // }

    
    //_______________________________________________________________________________________________________________________________

        // Double_t y = vm->Rapidity();
        // Double_t Mv = 3.09;
        // Double_t k = 0.5*Mv*TMath::Exp(TMath::Abs(y));
        // cout<<"Check Before neutron Generation"<<endl;
        
        // neutronArray EventNeutrons = gen->runSartreNoon(k);

        // for(Int_t side = 0; side <=1;side++)
        // {   
        //     vector<TLorentzVector> Ni;
        //     if(side == 0) Ni = EventNeutrons.n1Array;
        //     else Ni = EventNeutrons.n2Array;
        //     for(Long_t n = 0; n < sizeof(Ni); n++)
        //     {   
        //         cout<<"vec Parameters"<<endl;
        //         Ni[n].Print();
        //     }
        // }

//___________________________________________________________________________________________________________________________________
        
        //If the event is incoherent, and nuclear breakup is enabled, fill the remnants to the tree
        if(settings->enableNuclearBreakup() and event->diffractiveMode == incoherent){
            protons.Clear();
            neutrons.Clear();
            remnants.Clear();
            for(unsigned int iParticle=7; iParticle < event->particles.size(); iParticle++){
                if(event->particles[iParticle].status == 1) {  // Final-state particle
                    const Particle& particle = event->particles[iParticle];
                    switch (abs(particle.pdgId)) {
                        case 2212:  // (Anti-)proton
                            new(protons[protons.GetEntries()]) TLorentzVector(particle.p);
                            break;
                        case 2112:  // (Anti-)neutron
                            new(neutrons[neutrons.GetEntries()]) TLorentzVector(particle.p);
                            break;
                        default:  // Any other remnant
                            new(remnants[remnants.GetEntries()]) TLorentzVector(particle.p);
                            break;
                    }  // switch
                }  // if
            }  // for
        }  // if
        
        //
        // Decay the vector meson and fill the decay roducts in the tree:
        //
        if(doPerformDecay) {
            if( decay.SetDecay(*vm, 2, daughterMasses) ){
                double weight = decay.Generate(); // weight is always 1 here
                if ( (weight-1) > FLT_EPSILON) {
                    cout << "sartreMain: Warning weight != 1, weight = " << weight << endl;
                }
                vmDaughter1 = decay.GetDecay(0);
                vmDaughter2 = decay.GetDecay(1);
            }
            else {
                cout << "sartreMain: Warning: Kinematics of Vector Meson does not allow decay!" << endl;
            }
        }
        
        tree.Fill();
        // delete EventNeutrons;
    }
    cout << "All events processed\n" << endl;
    cout<<"Chekc noon"<<endl;
    //
    //  That's it, finish up
    //
    double totalCS=sartre.totalCrossSection();
    TH1D* histoForCSandNumberOfEvents = new TH1D("histoForCSandNumberOfEvents", "Cross-section and Number of Events", 2, 0., 1.);
    histoForCSandNumberOfEvents->SetBinContent(1, totalCS);
    histoForCSandNumberOfEvents->SetBinContent(2, maxEvents);
    
    double runTime = sartre.runTime();
    hfile->Write();
    cout << rootfile.c_str() << " written." << endl;
    cout << "Total cross-section: " << totalCS << " nb" << endl;
    sartre.listStatus();   
    cout << "CPU Time/event: " << 1000*runTime/maxEvents << " msec/evt" << endl;   
    
    return 0;   
}   

// UPC only
void randomlyReverseBeams(Event* myEvent) {
    
    TRandom3 *random = EventGeneratorSettings::randomGenerator();
    
    if(random->Uniform(1) > 0.5){
        for(unsigned int i=0; i<myEvent->particles.size(); i++)
            myEvent->particles.at(i).p.RotateX(M_PI);
    }
}


