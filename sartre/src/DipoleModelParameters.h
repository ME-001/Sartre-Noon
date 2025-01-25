//==============================================================================
//  DipoleModelParameters.h
//
//  Copyright (C) 2016-2019 Tobias Toll and Thomas Ullrich
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
//  $Date: 2019-09-06 11:36:10 -0400 (Fri, 06 Sep 2019) $
//  $Author: ullrich $
//==============================================================================
#ifndef DipoleModelParameters_h
#define DipoleModelParameters_h
#include "Enumerations.h"
#include <iostream>
#include <vector>

using namespace std;

class Settings;

class DipoleModelParameters {
    
public:
    DipoleModelParameters(Settings*);
    DipoleModelParameters(DipoleModelType, DipoleModelParameterSet);
    
    void setDipoleModelType(DipoleModelType);
    void setDipoleModelParameterSet(DipoleModelParameterSet);
    
    DipoleModelType dipoleModelType() const;
    DipoleModelParameterSet dipoleModelParameterSet() const;
    
    // bSat, bNonSat
    double BG() const;
    double mu02() const;
    double lambdaG() const;
    double Ag() const;
    double C() const;
    double rMax() const;  // STU only

    // bCGC
    double kappa() const;
    double N0() const;
    double x0() const;
    double lambda() const;
    double gammas() const;
    double Bcgc() const;
    
    double quarkMass(unsigned int) const;
    
    double boostedGaussianR2(int vmID);
    double boostedGaussianNL(int vmID);
    double boostedGaussianNT(int vmID);
    double boostedGaussianQuarkMass(int vmID);
    
    bool list(ostream& = cout); 
    
private:
    void setupParameters();
    void setup_bSat();
    void setup_bNonSat();
    void setup_bCGC();
    void setup_boostedGaussiansWaveFunction();
    
private:
    DipoleModelType         mDipoleModelType;
    string                  mDipoleModelParameterSetName;
    DipoleModelParameterSet mDipoleModelParameterSet;
    
    double mQuarkMass[6];

    // bSat, bNonSat
    double mBG;
    double mMu02;
    double mLambdaG;
    double mAg;
    double mC;
    double mRMax; // STU only

    // bCGC
    double mKappa;
    double mN0;
    double mX0;
    double mLambda;
    double mGammas;
    double mBcgc;
    
    // boosted Gaussian wave function parameters
    double mBoostedGaussianR2_rho;
    double mBoostedGaussianNL_rho;
    double mBoostedGaussianNT_rho;
    double mBoostedGaussianQuarkMass_rho;
    double mBoostedGaussianR2_phi;
    double mBoostedGaussianNL_phi;
    double mBoostedGaussianNT_phi;
    double mBoostedGaussianQuarkMass_phi;
    double mBoostedGaussianR2_jpsi;
    double mBoostedGaussianNL_jpsi;
    double mBoostedGaussianNT_jpsi;
    double mBoostedGaussianQuarkMass_jpsi;
    double mBoostedGaussianR2_ups;
    double mBoostedGaussianNL_ups;
    double mBoostedGaussianNT_ups;
    double mBoostedGaussianQuarkMass_ups;

    // hold the custom parameter (internal only)
    vector<double> mCustomParameters;
};

inline double DipoleModelParameters::BG() const {return mBG;}
inline double DipoleModelParameters::mu02() const {return mMu02;}
inline double DipoleModelParameters::C() const {return mC;}
inline double DipoleModelParameters::rMax() const {return mRMax;}
inline double DipoleModelParameters::lambdaG() const {return mLambdaG;}
inline double DipoleModelParameters::Ag() const {return mAg;}
inline double DipoleModelParameters::kappa() const {return mKappa;}
inline double DipoleModelParameters::N0() const {return mN0;}
inline double DipoleModelParameters::x0() const {return mX0;}
inline double DipoleModelParameters::lambda() const {return mLambda;}
inline double DipoleModelParameters::gammas() const {return mGammas;}
inline double DipoleModelParameters::Bcgc() const {return mBcgc;}
inline double DipoleModelParameters::quarkMass(unsigned int i) const
{
    if (i < 6)
        return mQuarkMass[i];
    else
        return 0;
}
#endif
