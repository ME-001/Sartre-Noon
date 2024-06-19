//==============================================================================
//  createBSatBDependenceTable.cpp
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
//  $Date: 2019-05-31 12:40:34 -0400 (Fri, 31 May 2019) $
//  $Author: ullrich $
//==============================================================================
//   
//  Main program to create impact parameter lookup table.   
//  [Developer only]  
//==============================================================================   
#include <iostream>  
#include <cmath>  
#include <sstream>  
#include <cstdlib>  
#include <vector>  
#include "TableGeneratorNucleus.h"  
#include "TableGeneratorSettings.h"  
#include "TFile.h"  
#include "TH2F.h"  
#include "TVector3.h"  
#include "Constants.h"  
#include "DipoleModelParameters.h"

#define PR(x) cout << #x << " = " << (x) << endl;

using namespace std;  

double overlapFunctionT(double);  
void createLookUpTableOnFile(unsigned int, unsigned int);  

DipoleModelParameters *gParameters;

int main(int argc, const char* argv[])  
{  
    if (argc != 3) {
        cout << "Usage:  " << argv[0] << " A  numberOfConfigurations" << endl;
        return 2;
    }
    
    gParameters = new DipoleModelParameters(bSat, KMW);
    
    //To explicitly delete all new histograms myself:
    TH1::AddDirectory(false);
    unsigned int A=atoi(argv[1]);
    unsigned int numConf=atoi(argv[2]);
    cout<<"Create table for A="<<A<<" and for "<<numConf<<" configurations."<<endl;
    createLookUpTableOnFile(A, numConf);
    return 0;
}  

double overlapFunctionT(double b)  
{  
    // Gaussian shape for proton
    // b in GeV
    double BG = gParameters->BG(); // GeV^-2
    double arg = (b*b/hbarc2) / (2*BG);
    return 1/(2*M_PI*BG) * exp(-arg);
}  


void createLookUpTableOnFile(unsigned int A, unsigned int numberOfConfigurations)  
{  
    //
    // A function to create a lookup table of nucleon configurations and save them on a file
    // This function is supposed to be independent from the rest of the code
    // and be called in the Main program with the nucleus number as an argument.
    //
    TableGeneratorNucleus myNucleus(A);
    //    myNucleus.init(A);
    TableGeneratorSettings::instance()->setSeed(42); // sets also seed for ROOT
    
    vector<TH2F*> hConfigurations(numberOfConfigurations);
    int numBbins=1e3, numAngleBins=1e3;
    double bRange=upperIntegrationLimit*myNucleus.radius()*1.1, angleRange=2*M_PI;
    ostringstream histoName;
    
    for(unsigned int iNuclei=0; iNuclei<numberOfConfigurations; iNuclei++) {
        histoName.str( "" );
        histoName << "Configuration_" << iNuclei;
        hConfigurations[iNuclei] = new TH2F(histoName.str().c_str(), "Sum_i^A T_p(b-b_i)",
                                            numBbins, 0., bRange, numAngleBins, 0., angleRange);
        //Generate a configuration:
        while(!myNucleus.generate()){}

        for(int ib=1; ib <= numBbins; ib++){
            double b=hConfigurations[iNuclei]->GetXaxis()->GetBinCenter(ib);
            for(int iphi=1; iphi <= numAngleBins; iphi++) {
                double phi=hConfigurations[iNuclei]->GetYaxis()->GetBinCenter(iphi);
                TVector3 bvector=TVector3(b*cos(phi), b*sin(phi), 0.);
                double sumOfT=0.;
                for(unsigned int iA=0; iA < A; iA++){
                    sumOfT+=overlapFunctionT((bvector-myNucleus.configuration.at(iA).position()).Perp());
                } //for iA
                hConfigurations[iNuclei]->SetBinContent(ib, iphi, sumOfT);
            }//for iphi
        } //for ib
        if ((numberOfConfigurations/10) && iNuclei % (numberOfConfigurations/10) == 0) {
            cout<<double(iNuclei)/numberOfConfigurations*100<<"% done."<<endl;
        }
    } //iNuclei
    cout<<"100% done."<<endl<<endl;
    cout<<"Writing configurations to file...";
    //Open the file and write to it:
    TFile *lufile = 0;
    ostringstream filename;
    filename.str("");
    filename << "bSat_bDependence_A" << A <<".root";
    lufile = new TFile(filename.str().c_str(), "RECREATE");
    for(unsigned int i=0; i<numberOfConfigurations; i++){
        hConfigurations[i]->Write();
        delete hConfigurations[i];
        hConfigurations[i] = 0;
    }
    lufile->Close();
    cout<<" done, bye!"<<endl;
}  

