//==============================================================================
//  PhotonFlux.h
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
//  $Date: 2019-11-23 17:21:53 -0500 (Sat, 23 Nov 2019) $
//  $Author: ullrich $
//==============================================================================
//         
//  Functor class.         
//  Photon flux is given in: d2sig/(dQ2 dW2)         
//==============================================================================       
#ifndef PhotonFlux_h         
#define PhotonFlux_h         
#include "Enumerations.h"
#include "EventGeneratorSettings.h"
#include "TH1D.h"
         
class Nucleus;

class PhotonFlux {         
public:         
    PhotonFlux();         
    PhotonFlux(double s);         
             
    void setS(double);         
             
    double operator()(double Q2, double W2, GammaPolarization p) const;         
    
    // UPC only:
    double operator()(double Egamma) const;         
    double nuclearPhotonFlux(double Egamma) const;
    
protected:
    double fluxTransverse(double Q2, double W2) const;
    double fluxLongitudinal(double Q2, double W2) const;
    double TAA(double);
    double TAAForIntegration(const double*) const;
    void   setupTAAlookupTable();
    
    // UPC only:
    double unintegratedNuclearPhotonFlux(double*, double*) const;
    double sigma_nn(double) const;

private:
    double mS;
    bool   mSIsSet;
    bool   mIsUPC;
    double mEBeamEnergy;
    double mB;
    TH1D*  mTAA_of_b;

    Nucleus* mNucleus;
    Nucleus* mNucleusUPC;
    
    EventGeneratorSettings *mSettings;
};         
#endif         
