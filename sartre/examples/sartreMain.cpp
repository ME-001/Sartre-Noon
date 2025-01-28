//==============================================================================
//  sartreMain.cpp
//
//  Copyright (C) 2010-2021 Tobias Toll and Thomas Ullrich
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
//  $Date: 2021-10-07 19:18:01 -0400 (Thu, 07 Oct 2021) $
//  $Author: ullrich $
//==============================================================================
//
//  Example main program. Useful to get started and easy to modify to fit
//  the user need. Data are stored in a simple but complete Root TTree.
//
//==============================================================================
#include <iostream>
#include <cmath>
#include "TTree.h"
#include "TFile.h"
#include "TGenPhaseSpace.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "Sartre.h"
#include "Settings.h"
#include "DipoleModelParameters.h"
#include "TwoBodyVectorMesonDecay.h"
#include "EicSmearFormatWriter.h"

// Neutron Generator related code
#include <fstream>
#include <string>
#include <vector>  

#include "TH1.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TDatabasePDG.h"
#include "TRandom.h"
#include "NeutronGenerator.h"

#define PR(x) cout << #x << " = " << (x) << endl;

//#define EIC_SMEAR_OUTPUT 1  // uncomment for eic-smear output

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
    int    nBeam1;		// number of neutrons in e beamside
    int    nBeam2;		// number of neutrons in p beamside
};

rootSartreEvent myRootSartreEvent;
void randomlyReverseBeams(Event* );  // for UPC mode

int main(int argc, char *argv[])
{
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
        cerr << "Error: Initialization of sartre failed." << endl;
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
    //  Derive filename for eic smear compatible output
    //  from root file
    //
    #if defined(EIC_SMEAR_OUTPUT)
    string eicSmearFilename;
    eicSmearFilename = rootfile.substr(0, rootfile.find(".root"));
    eicSmearFilename += ".out";
    EicSmearFormatWriter eicSmearWriter;
    eicSmearWriter.open(eicSmearFilename, settings->enableNuclearBreakup());
    cout << "eic-smear compatible output file is '" <<  eicSmearFilename.c_str() << "'." << endl;
    #endif

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
    
    TClonesArray protons("TLorentzVector");
    TClonesArray neutrons("TLorentzVector");
    TClonesArray remnants("TLorentzVector");
    TClonesArray neutrons1("TLorentzVector");
    TClonesArray neutrons2("TLorentzVector");
    
    TTree tree("tree","sartre");
    tree.Branch("event", &myRootSartreEvent.t,
                "t/D:Q2/D:x/D:s/D:y/D:W/D:xpom/D:iEvent/I:pol/I:dmode/I:nBeam1/I:nBeam2/I");
    tree.Branch("eIn",  "TLorentzVector", &eIn, 32000, 0);
    tree.Branch("pIn",  "TLorentzVector", &pIn, 32000, 0);
    tree.Branch("vm",   "TLorentzVector", &vm, 32000, 0);
    tree.Branch("eOut", "TLorentzVector", &eOut, 32000, 0);
    tree.Branch("pOut", "TLorentzVector", &pOut, 32000, 0);
    tree.Branch("gamma","TLorentzVector", &gamma, 32000, 0);
    tree.Branch("vmDaughter1", "TLorentzVector", &vmDaughter1, 32000, 0);
    tree.Branch("vmDaughter2", "TLorentzVector", &vmDaughter2, 32000, 0);
    
    // conditioned such that neutron branches are created for UPC only
    if(settings->UPC() and settings->A()==settings->UPCA()){   
    	tree.Branch("neutrons1", &neutrons1);
    	tree.Branch("neutrons2",&neutrons2);
    	//tree.Branch("Neutron_number",&myRootSartreEvent.nBeam1, "nBeam1/I:nBeam2/I"); // used for some testing 
    	
    }

    if(settings->enableNuclearBreakup()) {
        tree.Branch("protons", &protons);
        tree.Branch("neutrons", &neutrons);
        tree.Branch("nuclearRemnants", &remnants);
    }
    
    //
    //  Prepare event generation
    //
    TwoBodyVectorMesonDecay decayEngine;
    pair<TLorentzVector, TLorentzVector> daughters;
    int  daughterID = settings->userInt();
    bool doPerformDecay = false;
    
    if (daughterID && settings->vectorMesonId() != 22) {
        if (abs(daughterID) == 11 || abs(daughterID) == 13 ||
            abs(daughterID) == 211 || abs(daughterID) == 321) {
            doPerformDecay = true;
            cout << "Will decay vector meson: ";
            cout << settings->lookupPDG(settings->vectorMesonId())->GetName();
            cout << " -> ";
            cout << settings->lookupPDG(daughterID)->GetName();
            cout << " ";
            cout << settings->lookupPDG(-daughterID)->GetName();
            cout << endl;
        }
        else {
            cerr << "Error: Cannot decay vector meson to daughters with ID=" << daughterID << "." << endl;
            return 1;
        }
    }
    
    //
    //  Events and how often to show status
    //
    int nPrint;
    if (settings->timesToShow()){
        nPrint = settings->numberOfEvents()/settings->timesToShow();
    }
    else{
        nPrint = 0;
    }
    unsigned long maxEvents = settings->numberOfEvents();
    
    cout << "Generating " << maxEvents << " events." << endl << endl;
    
    //
    //creating a vector to store the generator object that will be created for only once in the loop
    // The idea is that it will create the required tables to use in the subsequent iterations.
    //
    
    vector<NeutronGenerator*> nGenVec; // vector to store the neutron generator object i.e. created in the first loop only
    
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
        if(settings->UPC() && settings->A() == settings->UPCA()){
            randomlyReverseBeams(event);
        }
        
        //
        //  Fill ROOT tree
        //
        myRootSartreEvent.iEvent = event->eventNumber;
        myRootSartreEvent.t  = event->t;
        myRootSartreEvent.Q2 = event->Q2;
        myRootSartreEvent.x  = event->x;
        myRootSartreEvent.y  = event->y;
        myRootSartreEvent.s  = event->s;
        myRootSartreEvent.W  = event->W;
        myRootSartreEvent.xpom  = event->xpom;
        myRootSartreEvent.pol   = event->polarization == transverse ? 0 : 1;
        myRootSartreEvent.dmode = event->diffractiveMode == coherent ? 0 : 1;
        myRootSartreEvent.nBeam1 = 0;
        myRootSartreEvent.nBeam2 = 0;
        *eIn     = event->particles[0].p;
        *pIn     = event->particles[1].p;
        *eOut    = event->particles[2].p;
        *pOut    = event->particles[6].p;
        *vm      = event->particles[4].p;
        *gamma   = event->particles[3].p;
        
        
        //
        //   Decay the vector meson and fill the decay products in the tree.
        //   For the decay we use the version that is sensitive to the
        //   polarization of the virtual photon (on a statistical basis).
        //   Also update "event" with the decay product.
        //

        if (doPerformDecay) {
            //
            //  Check if decay is possible. If not redo event.
            //
            double m = settings->lookupPDG(daughterID)->Mass();
            if (vm->M() < 2*m) {
                iEvent--;
                sartre.postAbortEvent();  // for Sartre internal accounting
                continue;
            }
            
            daughters = decayEngine.decayVectorMeson(*vm, *event, daughterID);
            *vmDaughter1 = daughters.first;
            *vmDaughter2 = daughters.second;
            
            Particle vmD1, vmD2;
            vmD1.index = 7;
            vmD1.pdgId = daughterID;
            vmD1.status = 1;
            vmD1.p = daughters.first;
            vmD1.parents.push_back(4);
            
            vmD2.index = 8;
            vmD2.pdgId = -daughterID;
            vmD2.status = 1;
            vmD2.p = daughters.second;
            vmD2.parents.push_back(4);

            event->particles[4].status = 2;  // mark J/psi as decayed
            event->particles[4].daughters.push_back(7);
            event->particles[4].daughters.push_back(8);
            
            //
            //  If there are no breakup fragments we simply
            //  add the decay daughters to the particle list.
            //  Otherwise we insert them before the fragments.
            //
            if (settings->enableNuclearBreakup() && event->diffractiveMode == incoherent) {
                vector<Particle>::iterator it = event->particles.begin();
                event->particles.insert(it+7, vmD1);
                it = event->particles.begin(); // it might have become invalid if reallocated
                event->particles.insert(it+8, vmD2);
            }
            else {
                event->particles.push_back(vmD1);
                event->particles.push_back(vmD2);
            }
        }
        
        
        if(settings->UPC() and settings->A()==settings->UPCA()){ // this is to ensure that this part of the code run only for UPC
            //loading the Neutron generator only once 
            if(iEvent < 1){
        	
        	    NeutronGenerator * nGen = new NeutronGenerator(eIn->Gamma());           // nGen = neutron Generator
        	    //gen->SetStoreQA();                                                  // Set kStoreQA = True; 
    		    //gen->SetStoreGeneratorFunctions();                                  // Store generator functions
    		    nGen->SetHadronicInteractionModel(NeutronGenerator::kGlauber);       // Selects model for ion interaction
    		    nGen->Initialize();                                                  // Loads default data (photon flux, cross-section, etc.)
    		    nGen->SetRunMode(NeutronGenerator::kInterface);                      // Run mode = kInterface for plain neutron generation
    		    //gen->ReadENDF(kTRUE); // By Default True in new file.               // Enables reading of ENDF data tables
    		    nGen->LoadENDF("hENDF.root");   // Need modification of this file.   // Data table for energy distribution among created neutrons
    		    //gen->Setup(); 
    		
    		    nGenVec.push_back(nGen);	
            }
        
            neutrons1.Clear();
            neutrons2.Clear();
        
            double y = vm->Rapidity();                                   // Store the rapidity corresponding to that entry
            //Double_t eta = vm->Eta();                                      // Store the pseudo rapidity corresponding to that entry

            double k = 0.5 * vm->M() * TMath::Exp(TMath::Abs(y));        // Calculate photon energy for vm production

            vector<int> nNumbers = nGenVec[0]->runSartreNoon(k);              // Gives the number of neutrons in both beams
        	
            event->particles.resize(7+nNumbers[0]+nNumbers[1]);            // resizing particle number
        
            //cout<<nNumbers[0]<<" "<<nNumbers[1]<<endl;
        
            //
            // Newly added code for storing neutron nombers in the event tree
            //
       
            myRootSartreEvent.nBeam1  = nNumbers[0];
            myRootSartreEvent.nBeam2  = nNumbers[1];
        
        
        
            if(nNumbers[0]!=0) {	
                vector<TLorentzVector> Vec; 
        	    nGenVec[0]->neutronRecord(nNumbers[0], nNumbers[1], 0, Vec);    // get the neutron record from the neutron generator for beam side 0
        	    for(int ineutron = 0; ineutron < nNumbers[0]; ineutron++){
	        	    Particle& particle = event->particles[7+ineutron];
        		    TLorentzVector* neutron1 = &Vec[ineutron];
        		    particle.p = *neutron1;
        		    particle.pdgId = 2112;
        		    particle.parents = {0};
        		    particle.status = 1;
			        new (neutrons1[neutrons1.GetEntriesFast()]) TLorentzVector(particle.p);
		        }
	        }
            if(nNumbers[1]!=0){	
        	    vector<TLorentzVector> Vec; 
        	    nGenVec[0]->neutronRecord(nNumbers[0], nNumbers[1], 1, Vec);  // get the neutron record from the neutron generator for beam side 1
        	    for(int ineutron = 0;ineutron < nNumbers[1]; ineutron++){
        	        Particle& particle = event->particles[6+nNumbers[0]+ineutron];
        	        TLorentzVector* neutron2= &Vec[ineutron];
        	        particle.p = *neutron2;
        	        particle.pdgId = 2112;
        	        particle.parents = {1};
        	        particle.status = 1;
		            new (neutrons2[neutrons2.GetEntriesFast()]) TLorentzVector(particle.p);
                }
	        }
	       
	    if (iEvent == maxEvents - 1) {
            	nGenVec.clear();  // Clear the vector of pointers
            	hfile->cd();      
        	}
	    }
	

        //
        //  If the event is incoherent, and nuclear breakup is enabled, fill the remnants to the tree
        //
        if (settings->enableNuclearBreakup() && event->diffractiveMode == incoherent){
            protons.Clear();
            neutrons.Clear();
            remnants.Clear();
            unsigned int startIndex = doPerformDecay ? 9 : 7;
            for(unsigned int iParticle=startIndex; iParticle < event->particles.size(); iParticle++){
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
        //  Print out (here only for the first few events)
        //
        if (iEvent < 10) event->list();
        
        //
        //  Fill and write event to file
        //
        tree.Fill();
        
        #if defined(EIC_SMEAR_OUTPUT)
        eicSmearWriter.writeEvent(event);
        #endif
    }
    // tree.Write();
    cout << "All events processed\n" << endl;
    
    //
    //  That's it, finish up
    //
    
    double totalCS=sartre.totalCrossSection();
    TH1D* histoForCSandNumberOfEvents = new TH1D("histoForCSandNumberOfEvents", "Cross-section and Number of Events", 2, 0., 1.);
    histoForCSandNumberOfEvents->SetBinContent(1, totalCS);
    histoForCSandNumberOfEvents->SetBinContent(2, maxEvents);
    
    double runTime = sartre.runTime();
    hfile->Write();
    cout << "File '" << rootfile << "' written." << endl;
    #if defined(EIC_SMEAR_OUTPUT)
    eicSmearWriter.close();
    cout << "File '" << eicSmearWriter.filename() << "' written." << endl;
    #endif
    cout << "Total cross-section: " << totalCS << " nb" << endl;
    sartre.listStatus();   
    cout << "CPU Time/event: " << 1000*runTime/maxEvents << " msec/evt" << endl;   
    
    delete histoForCSandNumberOfEvents;
    hfile->Close();
    return 0;   
}   

// UPC only
void randomlyReverseBeams(Event* myEvent)
{
    TRandom3 *random = EventGeneratorSettings::randomGenerator();
    
    if(random->Uniform(1) > 0.5){
        for(unsigned int i=0; i<myEvent->particles.size(); i++)
            myEvent->particles.at(i).p.RotateX(M_PI);
    }
}


