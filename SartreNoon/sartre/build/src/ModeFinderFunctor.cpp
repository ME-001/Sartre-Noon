//==============================================================================
//  ModeFinderFunctor.cpp
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
//  $Date: 2019-03-08 14:12:33 -0500 (Fri, 08 Mar 2019) $
//  $Author: ullrich $
//==============================================================================
#include "ModeFinderFunctor.h"    
#include "CrossSection.h"    
#include "Kinematics.h"    
#include <iostream>

using namespace std;    

#define PR(x) cout << #x << " = " << (x) << endl;

ModeFinderFunctor::ModeFinderFunctor()    
{    
    mCrossSection = 0;    
    mQ2 = mVmMass = mMinT = mMaxT = 0;    
}    
    
ModeFinderFunctor::ModeFinderFunctor(CrossSection* cs, double Q2, double vmMass, double tmin, double tmax)    
{    
    mCrossSection = cs;    
    mQ2 = Q2;    
    mVmMass = vmMass;    
    mMinT = tmin;  
    mMaxT = tmax;   
}    
      
void ModeFinderFunctor::setVmMass(double val) {mVmMass = val;}    
    
void ModeFinderFunctor::setQ2(double val) {mQ2 = val;}    
   
void ModeFinderFunctor::setMinT(double val) {mMinT = val;}    
    
void ModeFinderFunctor::setMaxT(double val) {mMaxT = val;}    
  
double ModeFinderFunctor::DoEval(double W2) const    
{    
    if (mCrossSection) {    
        double t = Kinematics::tmax(0, mQ2, W2, mVmMass); // first arg (t) set to 0 here    
        if (t > mMaxT) t = mMaxT; // don't exceed given (table) maximum (smallest |t|)  
        if (t < mMinT) return 0;  // lower table limits (treat different than max limit)  
        double result = (*mCrossSection)(t, mQ2, W2); // t, Q2, W2    
        return -result; // minimum=maximum    
    }    
    else {     
        return 0;    
    }    
}    
    
ROOT::Math::IGenFunction* ModeFinderFunctor::Clone() const    
{    
    return new ModeFinderFunctor(mCrossSection, mQ2, mVmMass, mMinT, mMaxT);    
}    

UPCModeFinderFunctor::UPCModeFinderFunctor()
{
    mCrossSection = 0;
    mVmMass = mHBeamEnergy = mEBeamEnergy = 0;
}

UPCModeFinderFunctor::UPCModeFinderFunctor(CrossSection* cs, double vmMass, double hEnergy, double eEnergy)
{
    mCrossSection = cs;
    mVmMass = vmMass;
    mHBeamEnergy = hEnergy;
    mEBeamEnergy = eEnergy;
}
    
double UPCModeFinderFunctor::DoEval(double val) const
{
    if (mCrossSection) {
        double xpom = exp(val);  // x comes as log(x)
        double t = Kinematics::tmax(xpom);
        bool ok = Kinematics::validUPC(mHBeamEnergy, mEBeamEnergy, t, xpom, mVmMass, false);
        if (ok) {
            double result = (*mCrossSection)(t, xpom);
            return -result; // minimum=maximum
        }
        else {
            return 0;
        }
    }
    else {
        return 0;
    }
}

ROOT::Math::IGenFunction* UPCModeFinderFunctor::Clone() const
{
    return new UPCModeFinderFunctor(mCrossSection, mVmMass, mHBeamEnergy, mEBeamEnergy);
}
