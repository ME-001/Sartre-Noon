//==============================================================================
//  Kinematics.h
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
//
//------------------------------------------------------------------------------
//
//  Note that Kinematics is only meant for ep and eA diffactive events.
//  It is not applicable for UPC kinematics. Use UpcKinematics instead.
//
//==============================================================================
#ifndef Kinematics_h         
#define Kinematics_h         
#include "Constants.h"         
#include "Settings.h"         
#include "TLorentzVector.h"         
#include "Math/IFunction.h"    
         
class Kinematics {         
public:         
    static double         xprobe(double Q2, double W2, double vmM);         
    static TLorentzVector electronBeam(double eE, bool upc = false);
    static TLorentzVector hadronBeam(double eH);         
    static double         s(const TLorentzVector& e, const TLorentzVector& p);         
    static double         s(double eE, double hE, bool upc = false);
    static double         x(double Q2, double W2);             
    static double         y(double Q2, double x, double s);             
    static double         ymin(double s, double vmM);             
    static double         ymax(double s, double vmM);             
    static double         W2(double Q2, double x);         
    static double         W2(double Q2, double xprobe, double vmM);         
    static double         W2min(double vmM);         
    static double         W2max(double s);         
    static double         Q2min(double y);         
    static double         Q2max(double s);         
    static double         xpomeron(double t, double Q2, double W2, double vmM);         
    static double         tmax(double xpom); // often referred to as tmin implying min |T|         
    static double         tmax(double t, double Q2, double W2, double vmM);         
    static double         tmin(double hE);   // the smallest t value possible          
    
    static double         Egamma(double xpom, double t, double vmM, double hBeamEnergy, double eBeamEnergy,  double MY2minusM2 = 0);

    static bool           validUPC(double hBeamEnergy, double eBeamEnergy, double t,
				                   double xpom, double vmMass, bool verbose = false);   // UPC only
    static bool           valid(double s, double t, double Q2, double W2, double vmMass,          
                                bool useTrueXp = false,          
                                bool verbose = false);         
    
    static bool           error();

    static double         xpomMin(double massVM, double t, TLorentzVector hBeam, TLorentzVector eBeam,  double MY2minusM2 = 0);   // UPC only

private:
    static double         xpomMin(double s, double vmM, double t);   // UPC only

private:         
    static bool   mError;

    static double mXpomMin;
    static double mTold;
    static bool   mXpomMinIsEvaluated;
};         


//-------------------------------------------------------------------------------    
//    
//  Helper class needed to find root in
//  Kinematics::validUPC()
//
//-------------------------------------------------------------------------------    

class LowerXpomeronFormula : public ROOT::Math::IBaseFunctionOneDim
{
public:    
    double DoEval(double) const;
    ROOT::Math::IBaseFunctionOneDim* Clone() const;
    void calculateValidRange(double&, double&);
public:        
    double mT;
    double mVmMass2;
    double mXpomMin;
    double mMY2;
    TLorentzVector mElectronBeam;
    TLorentzVector mHadronBeam;
    double mMY2minusM2;
};

#endif         
