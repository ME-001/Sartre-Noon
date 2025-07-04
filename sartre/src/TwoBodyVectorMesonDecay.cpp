//==============================================================================
//  TwoBodyVectorMesonDecay.cpp
//
//  Copyright (C) 2019 Tobias Toll and Thomas Ullrich
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
//  $Date:$
//  $Author:$
//==============================================================================
#include "TwoBodyVectorMesonDecay.h"
#include "TGenPhaseSpace.h"
#include "TVector3.h"
#include <limits>
#include <utility>

#define PR(x) cout << #x << " = " << (x) << endl;

TwoBodyVectorMesonDecay::TwoBodyVectorMesonDecay()
{
    mSettings = EventGeneratorSettings::instance();
    mRandom = EventGeneratorSettings::randomGenerator();
}

pair<TLorentzVector, TLorentzVector>
TwoBodyVectorMesonDecay::decayVectorMeson(TLorentzVector& vm, int daughterID)
{
    //
    //  Simple 2-body decay flat in phase space
    //
    TGenPhaseSpace decay;
    
    pair<TLorentzVector, TLorentzVector> daughters;
    double daughterMass = mSettings->lookupPDG(daughterID)->Mass();
    double daughterMasses[2] = {daughterMass, daughterMass};
    
    if (decay.SetDecay(vm, 2, daughterMasses)) {
        double weight = decay.Generate(); // weight is always 1 here
        if ((weight-1) > numeric_limits<float>::epsilon()) {
            cout << "TwoBodyVectorMesonDecay::decayVectorMeson(): Warning weight != 1, weight = " << weight << endl;
        }
        daughters.first = *decay.GetDecay(0);
        daughters.second = *decay.GetDecay(1);
    }
    else {
        cout << "TwoBodyVectorMesonDecay::decayVectorMeson(): Warning, kinematics of vector meson does not allow decay!" << endl;
    }
    
    return daughters;
}

pair<TLorentzVector, TLorentzVector>
TwoBodyVectorMesonDecay::decayVectorMeson(TLorentzVector& vm, Event& event, int daughterID)
{
    //
    //   Decay correctly treated with polarization of the virtual photon taken into account
    //   (SCHC approximation).
    //   This code was developed with the help of Athira Vijayakumar and Barak Schmookler
    //   from Stony Brook University.
    //
    //   Essentially we create the daughters in the rest frame of the vector meson.
    //   Angels are generated randomly, flat in phi, and cos(theta) according to the
    //   respective distribution found in literature (see cosTheta()). At the end we
    //   boost back in the lab system.
    //
    
    pair<TLorentzVector, TLorentzVector> daughters;
    
    double daughterMass = mSettings->lookupPDG(daughterID)->Mass();
    
    //
    //  Get matrix element from the cross-section ratio sig_L/sig_T.
    //  An alternative is to not use the ratio from Sartre directly but calculate
    //  it using Eq. 16 in PLB 449, 328. However, this is (i) model dependent and
    //  (ii) would introduce value that is not always identical to the one from
    //  Sartre. Tests by Athira showed that the two are very close so we use the
    //  actual value generated by Sartre, which is stored in the Event structure.
    //  Note also that the PLB version is t-integrated while Sartre's one is not.
    //
    double polarizationParameter = (1 - event.y)/(1 - event.y + (event.y*event.y)/2);
    double a = event.crossSectionRatioLT * polarizationParameter;
    double r0400 = a/(1+a);
    
    //
    //  Generate random angles with the proper distributions
    //
    TRandom3 *random = EventGeneratorSettings::randomGenerator();
    double phi = random->Uniform(2*M_PI);
    double costh = cosTheta(r0400, daughterID);
    if (fabs(costh) > 1) {
        cout << "TwoBodyVectorMesonDecay::decayVectorMeson(): Error, falling back to simple decay scheme." << endl;
        return decayVectorMeson(vm, daughterID);
    }
    double theta = acos(costh);
    
    //
    //  Daughters in center-of-mass of vector meson,
    //  relying on both daughters having the same mass (and spin).
    //
    double psquare = vm.M2()/4 - daughterMass*daughterMass;
    if (psquare < 0) {
        cout << "TwoBodyVectorMesonDecay::decayVectorMeson(): Error, vector meson mass too small for decay into given "
                "daughters. Will result in 'nan' values." << endl;
    }
    double p = sqrt(psquare);
    double E = sqrt(p*p + daughterMass*daughterMass);
    TVector3 p3(sin(theta)*cos(phi)*p, sin(theta)*sin(phi)*p, cos(theta)*p);
    daughters.first = TLorentzVector(p3, E);
    daughters.second = TLorentzVector(-p3, E);
    
    //
    //  Boost to lab
    //
    TVector3 labVec = vm.Vect()*(1/vm.Energy());
    daughters.first.Boost(labVec);
    daughters.second.Boost(labVec);
    
    return daughters;
}

double TwoBodyVectorMesonDecay::cosTheta(double r, int daughterID)
{
    //
    //  cos(theta) is generated according to the distributions derived by H1 et al.
    //  See: Eur. Phys. J. C 6, 603 (1999) and Eur. Phys. J. C 13, 371 (2000)
    //
    //  To generate a random number we invert the cumulative of this distribution.
    //  It ends up with a (reduced) cubic equation that is solved analytically.
    //  We cover several cases since the solution scheme depends on the matrix
    //  element r (R^04_00).
    //
    //  If something goes wrong we return 42 instead of -1 < cos(theta) < 1.
    //
    double p, q;
    
    double rnd = mRandom->Uniform(1);
    
    if (abs(daughterID) == 211 || abs(daughterID) == 321) {     //  spin 0 decay particle
        p = (3-3*r)/(3*r-1);
        q = 1+(3-3*r-4*rnd)/(3*r-1);
    }
    else if (abs(daughterID) == 11 || abs(daughterID) == 13) {  // spin 1/2 decay particle
        p = (3+3*r)/(1-3*r);
        q = (4-8*rnd)/(1-3*r);
    }
    else {
        cout << "TwoBodyVectorMesonDecay::cosTheta(): Cannot handle particle ID=" << daughterID << "." << endl;
        return 42;
    }
    
    double R = sqrt(fabs(p)/3.)*(q < 0 ? -1 : 1);
    double D = pow(p/3, 3) + q*q/4;
    
    double result, phi;
    
    if (p > 0) {
        phi = asinh(q/(2*R*R*R));
        result = -2*R*sinh(phi/3);
    }
    else if (p < 0) {
        if (D <= 0) {
            phi = acos(q/(2*R*R*R));
            //result = -2*R*cos(phi/3);
            //result = -2*R*cos(phi/3 + 2*M_PI/3);
            result = -2*R*cos(phi/3 + 4*M_PI/3);
        }
        else {
            phi = acosh(q/(2*R*R*R));
            result = -2*R*cosh(phi/3);
        }
    }
    else {
        if (-q >= 0)
            result = pow(-q, 1./3.);
        else {
            cout << "TwoBodyVectorMesonDecay::cosTheta(): Warning, cannot generate cos(theta)." << endl;
            return 42;
        }
    }
    
    return result;
}
