//==============================================================================
//  Nucleus.cpp
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
//  $Date: 2021-07-02 14:27:46 -0400 (Fri, 02 Jul 2021) $
//  $Author: ullrich $
//==============================================================================
#include "Nucleus.h"    
#include "Constants.h"    
#include "TF1.h"    
#include "TH1.h"    
#include "TH1D.h"
#include <cmath>
#include <cstdlib>    
  
#define PR(x) cout << #x << " = " << (x) << endl;    

Nucleus::Nucleus()
{    
    mA = mZ = 0;    
    mRadius = 0;    
    mSpin = 0;    
    mSurfaceThickness = 0;    
    mRho0 = 0;    
    mOmega = 0;    
    mMass = 0;
    mHulthenA = 0;
    mHulthenB = 0;
}
    
Nucleus::Nucleus(unsigned int A)     
{
    mA = mZ = 0;
    mRadius = 0;
    mSpin = 0;
    mSurfaceThickness = 0;
    mRho0 = 0;
    mOmega = 0;
    mMass = 0;
    mHulthenA = 0;
    mHulthenB = 0;
    
    Nucleus::init(A);
}    
    
Nucleus& Nucleus::operator=(const Nucleus& n)  
{  
    if (this != &n) {
        mA = n.mA;
        mZ = n.mZ;    
        mRadius = n.mRadius;    
        mSpin = n.mSpin;    
        mSurfaceThickness = n.mSurfaceThickness;    
        mRho0 = n.mRho0;    
        mOmega = n.mOmega;    
        mMass = n.mMass;      
        mName = n.mName;      
        mHulthenA = n.mHulthenA;
        mHulthenB = n.mHulthenB;
        mLookupTable = unique_ptr<TH1D>(new TH1D(*(n.mLookupTable)));
        mLookupTable->SetDirectory(0);  
    }  
    return *this;  
}  
  
Nucleus::Nucleus(const Nucleus& n)  
{  
    mA = n.mA;
    mZ = n.mZ;    
    mRadius = n.mRadius;    
    mSpin = n.mSpin;    
    mSurfaceThickness = n.mSurfaceThickness;    
    mRho0 = n.mRho0;    
    mOmega = n.mOmega;    
    mMass = n.mMass;      
    mName = n.mName;      
    mHulthenA = n.mHulthenA;
    mHulthenB = n.mHulthenB;
    mLookupTable = unique_ptr<TH1D>(new TH1D(*(n.mLookupTable)));
    mLookupTable->SetDirectory(0);
}  
  
Nucleus::~Nucleus()
{
    //
    //  This is a temorary fix of a problem in the code that is not understood.
    //  In the table generation code after all is calculated the ~TH1 of the
    //  look up table histogram causes a crash. This appears (?) to be related
    //  to issues in GSL and ROOT. Releasing the pointer causes a memory leak
    //  but this is not a big problem since it happens at the end of the program.
    //
    mLookupTable.release();
}
    
void Nucleus::init(unsigned int A)     
{    
    mA = A;    
        
    //    
    // Set the parameters of the nucleus    
    // Woods-Saxon parameter are from Ramona Vogt's paper on nuclear geometry (table 1).
    // The parametrization is a 3 parameter fermi model (or 2 parameter with omega=0).
    // One exception is the deuteron where a Hulthen function is used (nucl-ex/0701025v1).
    // All version are defined such that Integral(d3r rho(r)) = A
    // Masses from https://wwwndc.jaea.go.jp/NuC/
    //
    const double unifiedAtomicMass = 0.931494013;    
    switch (mA) {    
        case 208:   // Pb    
            mZ = 82;    
            mRadius = 6.624;     
            mSurfaceThickness = 0.549;     
            mOmega = 0;    
            mRho0 = 0.160;     
            mName = "Pb";    
            mMass = 207.976651189*unifiedAtomicMass-mZ*electronMass;
            mSpin = 0;
            break;    
        case 197:   // Au    
            mZ = 79;    
            mRadius = 6.38;     
            mSurfaceThickness = 0.535;     
            mOmega = 0;    
            mRho0 = 0.1693;     
            mName = "Au";    
            mMass = 196.966568812*unifiedAtomicMass-mZ*electronMass;
            mSpin = 3./2.;    
            break;    
        case 110: //  Cd    
            mZ = 48;    
            mRadius = 5.33;     
            mSurfaceThickness = 0.535;      
            mOmega = 0;    
            mRho0 = 0.1577;     
            mName = "Cd";    
            mMass = 109.903004927*unifiedAtomicMass-mZ*electronMass;
            mSpin = 0;
            break;
        case 90:   // Zr
            mZ = 40;
            mRadius = 4.9736093; // from R. Vogt 1.19*A^(1/3)-1.61*A^(-1/3)
            mSurfaceThickness = 0.54; // from R. Vogt
            mOmega = 0; // from R. Vogt
            mRho0 = 0.15643742;
            mName = "Zr";
            mMass = 89.904696939*unifiedAtomicMass-mZ*electronMass;
            mSpin = 0;
            break;
        case 63: //  Cu    
            mZ = 29;    
            mRadius = 4.214;     
            mSurfaceThickness = 0.586;      
            mOmega = 0;    
            mRho0 = 0.1701;     
            mName = "Cu";    
            mMass = 62.929597770*unifiedAtomicMass-mZ*electronMass;
            mSpin = 3./2.;    
            break;    
        case 40: // Ca    
            mZ = 20;    
            mRadius = 3.766;     
            mSurfaceThickness = 0.586;      
            mOmega = -0.161;    
            mRho0 = 0.1699;     
            mName = "Ca";    
            mMass= 39.962590863*unifiedAtomicMass-mZ*electronMass;
            mSpin = 0;
            break;    
        case 27: //  Al    
            mZ = 13;    
            mRadius = 3.07;     
            mSurfaceThickness = 0.519;      
            mOmega = 0;    
            mRho0 = 0.1739;     
            mName = "Al";    
            mMass = 26.981538578*unifiedAtomicMass-mZ*electronMass;
            mSpin = 5./2.;    
            break;    
        case 16: //  O    
            mZ = 8;    
            mRadius = 2.608;     
            mSurfaceThickness = 0.513;      
            mOmega = -0.051;    
            mRho0 = 0.1654;     
            mName = "O";    
            mMass = 15.99491461957*unifiedAtomicMass-mZ*electronMass;
            mSpin = 0;
            break;    
        case 2: // D
            mZ = 1;
            mHulthenA = 0.456;
            mHulthenB = 2.36;
            mRho0 = 0.39946312;
            mName = "D";
            mMass = 2.01410177812*unifiedAtomicMass-mZ*electronMass;
            mSpin = 1;
            mRadius = 2.116;  // not used for density function - taken from Jager et al.
            break;
        case 1:
            // When generating this nucleus, only one nucleon is put at (0, 0, 0)    
            mZ = 1;    
            mRadius = 1.;    
            mMass = protonMass; // 1.00794*unifiedAtomicMass    
            mName = "p";    
            //The following 3 are dummies, needed for radial distribution etc.:    
            mSurfaceThickness = 0.513;      
            mOmega = -0.051;    
            mRho0 = 0.1654;     
            mSpin = 1./2.;    
            break;    
        default:    
            cout << "Nucleus::init(): Error, cannot handle A=" << mA << ", no parameters defined for this mass number." << endl;    
            exit(1);    
            break;    
    }       
        
    //    
    //  Generate lookup table for overlap integral T    
    //  (b and z in fm)    
    //    
    double range = 3*mRadius;    
    TF1* funcToIntegrate = new TF1("funcToIntegrate",this, &Nucleus::rhoForIntegration, -range, range, 1);
    int nbins = 5000; // size of lookup table
    mLookupTable = unique_ptr<TH1D>(new TH1D("mLookupTable","T lookup table", nbins, 0, range));
    mLookupTable->SetDirectory(0);
    for (int i=1; i<=nbins; i++) {
        double b = mLookupTable->GetBinCenter(i);    
        funcToIntegrate->SetNpx(1000);    
        funcToIntegrate->SetParameter(0, b);    
        double res = funcToIntegrate->Integral(-range, range, 1e-8);    
        mLookupTable->SetBinContent(i, res);        
    }    
    delete funcToIntegrate;
}
    
unsigned int Nucleus::A() const { return mA; }    
    
unsigned int Nucleus::Z() const { return mZ; }    
    
float Nucleus::spin() const { return mSpin; }    
    
double Nucleus::nuclearMass() const { return mMass; }    
    
string Nucleus::name() const { return mName; }    
    
double Nucleus::radius() const { return mRadius; }    
    
double Nucleus::rho(double b, double z) // returns value in units of fm^-3    
{    
    // b and z in fm    
    double r = sqrt(b*b+z*z);
    if (mA == 2) {
        double term = (exp(-mHulthenA*r)/r) - (exp(-mHulthenB*r)/r);  // Hulthen
        return mRho0*term*term;
    }
    else {
        return mRho0*(1+mOmega*(r/mRadius)*(r/mRadius))/(1+exp((r-mRadius)/mSurfaceThickness)); //  3pF
    }
}    
    
double Nucleus::rhoForIntegration(double *x, double* par)    
{    
    // b and z in fm    
    double b = par[0];    
    double z = *x;    
    return rho(b, z);    
}    
    
double Nucleus::T(double b) const    
{    
    //     
    //  Returns overlap integral     
    //  b in fm, overlap function T in GeV^2    
    //  Returns 0 if b exceeds table range    
    //    
    if (mA == 0) {    
        cout << "Nucleus::T(): Error, instance is not initialized - cannot calculate overlap integral." << endl;    
        return 0;    
    }    
    int bin = mLookupTable->FindBin(b);
    double res = mLookupTable->GetBinContent(bin);
    return res*hbarc2;
}    
    
double Nucleus::TForIntegration(double *x, double*) const  
{    
    double b = x[0];    
    return 2*b*M_PI*T(b)/hbarc2;    
}    
    
void Nucleus::normalizationOfT(double eps)
{    
    //    
    //  This is for internal checks only.    
    //  Normalization of rho/T is such that fully integrated    
    //  function should yield A.    
    //    
    double range = 3*mRadius;    
    TF1* func = new TF1("func",this, &Nucleus::TForIntegration, 0, range, 0);    
    double res = func->Integral(0, range, eps);
    cout << "Normalization of T: = " << res << endl;
    delete func;
}    
    
double Nucleus::rho0() const { return mRho0; }

double Nucleus::TofProton(double b)    
{    
    //    
    //  Gaussian shape for proton    
    //  Units:    
    //  b in fm    
    //    
    const double BG = 4; // GeV^-2    
    double arg = b*b/(2*BG);    
    arg /= hbarc2;    
    return 1/(2*M_PI*BG) * exp(-arg);        
}    
    
int Nucleus::pdgID(int Z, int A) const    
{    
    // PDG for nuclei 10LZZZAAAI    
    // L = number of strange quark (for hypernuclei)    
    // I = isomer level, with I = 0 corresponding to the     
    //     ground state and I > 0 to excitations,    
    if (Z==1 && A==1)     
        return 2212;    
    else if (Z==0 && A==1)    
        return 2112;    
    else if (A > 1)    
        return 1000000000 + 10000*Z + 10*A;      
    else     
        return 0;    
}    
    
int Nucleus::pdgID() const    
{    
    return pdgID(mZ, mA);    
}    
  
