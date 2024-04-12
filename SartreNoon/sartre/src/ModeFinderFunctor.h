//==============================================================================
//  ModeFinderFunctor.h
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
//         
//  Functor class.        
//         
//  Used by BrentMinimizer1D algotithm to find maxima in cross-section         
//  in W2 for a given t and Q2.         
//==============================================================================        
#ifndef ModeFinderFunctor_h         
#define ModeFinderFunctor_h         
#include "Math/IFunction.h"         
                  
class CrossSection;         
         
class ModeFinderFunctor : public ROOT::Math::IGenFunction {         
public:         
    ModeFinderFunctor();         
    ModeFinderFunctor(CrossSection*, double Q2, double vmMass, double tmin, double tmax);         
             
    double DoEval(double) const;         
    ROOT::Math::IGenFunction* Clone() const;          
    void setQ2(double);         
    void setVmMass(double);         
    void setMaxT(double);         
    void setMinT(double);         
             
private:         
    CrossSection *mCrossSection;         
    double mQ2;         
    double mVmMass;  
    double mMinT;  
    double mMaxT;    
};

class UPCModeFinderFunctor : public ROOT::Math::IGenFunction {
public:
    UPCModeFinderFunctor();
    UPCModeFinderFunctor(CrossSection*, double vmMass, double hBeamEnergy, double eBeamEnergy);
    
    double DoEval(double) const;  // xpom
    ROOT::Math::IGenFunction* Clone() const;

private:
    CrossSection *mCrossSection;
    double mVmMass;
    double mHBeamEnergy;
    double mEBeamEnergy;
};
#endif         
         
