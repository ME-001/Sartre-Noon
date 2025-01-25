//==============================================================================
//  CrossSection.cpp
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
//  $Date: 2019-12-18 11:55:18 -0500 (Wed, 18 Dec 2019) $
//  $Author: ullrich $
//==============================================================================
#include "CrossSection.h"    
#include "EventGeneratorSettings.h"    
#include "Kinematics.h"    
#include "TableCollection.h"    
#include "Table.h"    
#include "Constants.h"    
#include "TH2F.h"
#include "TFile.h"
#include "DglapEvolution.h"
#include <cstdio>
#include <cmath>
#include <limits>

#define PR(x) cout << #x << " = " << (x) << endl;    
  
CrossSection::CrossSection(TableCollection* tc, TableCollection* ptc)    
{    
    mSettings = EventGeneratorSettings::instance();    
    mRandom = mSettings->randomGenerator();    
    mS = Kinematics::s(mSettings->eBeam(), mSettings->hBeam());    
    mPhotonFlux.setS(mS);    
    mTableCollection = tc;    
    mProtonTableCollection = ptc;    
    TParticlePDG *vectorMesonPDG = mSettings->lookupPDG(mSettings->vectorMesonId());        
    mVmMass = vectorMesonPDG->Mass();    
    mCheckKinematics = true;
    mCrossSectionRatioLT = 0;
}    
    
CrossSection::~CrossSection() { /* no op */ }    
  
void CrossSection::setTableCollection(TableCollection* tc) {mTableCollection = tc;}    

void CrossSection::setProtonTableCollection(TableCollection* ptc) {mProtonTableCollection = ptc;}    
    
void CrossSection::setCheckKinematics(bool val) {mCheckKinematics = val;}  
  
double CrossSection::operator()(double t, double Q2, double W2)
{
    return dsigdtdQ2dW2_total_checked(t, Q2, W2);
}

// UPC Version
double CrossSection::operator()(double t, double xpom)
{
    return dsigdtdxp_total_checked(t, xpom);
}

double CrossSection::operator()(const double* array)
{
    if (mSettings->UPC())
        return dsigdtdxp_total_checked(array[0], array[1]);
    else
        return dsigdtdQ2dW2_total_checked(array[0], array[1], array[2]);
}

//
//   PDF passed to UNURAN
//   Array holds: t, log(Q2), W2 for e+p/A running
//                t, log(xpom)  for UPC
//
double CrossSection::unuranPDF(const double* array)    // array is t, log(Q2), W2
{                                                      // or t and log(xpom)
    double result = 0;
    if (mSettings->UPC()) {
        double xpom = exp(array[1]);
        result = dsigdtdxp_total_checked(array[0], xpom);  // t, xpom
        result *= xpom; // Jacobian
    }
    else {
        double Q2 = exp(array[1]);
        result = dsigdtdQ2dW2_total_checked(array[0], Q2, array[2]);  // t, Q2, W2
        result *= Q2;   // Jacobian
    }
    return log(result);
}

double CrossSection::dsigdtdQ2dW2_total_checked(double t, double Q2, double W2)
{
    double result = 0;
        
    //    
    //  Check if in kinematically allowed region    
    //    
    if (mCheckKinematics && !Kinematics::valid(mS, t, Q2, W2, mVmMass, false, (mSettings->verboseLevel() > 2) )) {    
        if (mSettings->verboseLevel() > 2)     
            cout << "CrossSection::dsigdtdQ2dW2_total_checked(): warning, t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2)
                << " is outside of kinematically allowed region. Return 0." << endl;
        return result;    
    }    

    //    
    //  Total cross-section dsig2/(dQ2 dW2 dt)    
    //  This is the probability density function needed for UNU.RAN    
    //    
    double csT = dsigdtdQ2dW2_total(t, Q2, W2, transverse);    
    double csL = dsigdtdQ2dW2_total(t, Q2, W2, longitudinal);    
    result = csT + csL;
    mCrossSectionRatioLT = csL/csT;

    //
    //  Polarization    
    //    
    if (mRandom->Uniform(result) <= csT)    
        mPolarization = transverse;    
    else    
        mPolarization = longitudinal;    
        
    //    
    // Diffractive Mode    
    //    
    double sampleRange = (mPolarization == transverse ? csT : csL);    
    double sampleDivider = dsigdtdQ2dW2_coherent(t, Q2, W2, mPolarization);    
    if (mRandom->Uniform(sampleRange) <= sampleDivider)    
        mDiffractiveMode = coherent;    
    else    
        mDiffractiveMode = incoherent;    
     
    //    
    // Print-out at high verbose levels    
    //    
    if (mSettings->verboseLevel() > 10) {     // Spinal Tap ;-)
        cout << "CrossSection::dsigdtdQ2dW2_total_checked(): " << result;
        cout << " at t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2);    
        cout << " (" << (mPolarization == transverse ? "transverse" : "longitudinal");    
        cout << " ," << (mDiffractiveMode == coherent ? "coherent" : "incoherent");    
        cout << ')' << endl;    
    }    
    
    //    
    // Check validity of return value    
    //    
    if (std::isnan(result)) {    
        cout << "CrossSection::dsigdtdQ2dW2_total_checked(): Error, return value = NaN at" << endl;
        cout << "                                            t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2) << endl;
        result = 0;            
    }    
    if (std::isinf(result)) {    
        cout << "CrossSection::dsigdtdQ2dW2_total_checked(): Error, return value = inf at" << endl;
        cout << "                                            t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2) << endl;
        result = 0;    
    }    
    if (result < 0) {    
        cout << "CrossSection::dsigdtdQ2dW2_total_checked(): Error, negative cross-section at" << endl;
        cout << "                                            t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2) << endl;
        result = 0;    
    }    
        
    return result;    
}    

double CrossSection::dsigdtdxp_total_checked(double t, double xpom)
{
    double result = 0;
    
    //
    //  Check if in kinematically allowed region
    //
    if (mCheckKinematics && !Kinematics::validUPC(mSettings->hadronBeamEnergy(),
                                                  mSettings->electronBeamEnergy(),
                                                  t, xpom, mVmMass,
                                                  (mSettings->verboseLevel() > 2) )) {
        if (mSettings->verboseLevel() > 2)
            cout << "CrossSection::dsigdtdxp_total_checked(): warning, t=" << t << ", xpom=" << xpom
                 << " is outside of kinematically allowed region. Return 0." << endl;
        return result;
    }
    
    //
    //  Total cross-section dsig2/(dt dxp)
    //  This is the probability density function needed for UNU.RAN
    //
    result = dsigdtdxp_total(t, xpom);
    
    //
    //  Polarization
    //
    mPolarization = transverse; // always
    
    //
    // Diffractive Mode
    //
    double sampleRange = result;
    double sampleDivider = dsigdtdxp_coherent(t, xpom);
    if (mRandom->Uniform(sampleRange) <= sampleDivider)
        mDiffractiveMode = coherent;
    else
        mDiffractiveMode = incoherent;
    
    //
    // Print-out at high verbose levels
    //
    if (mSettings->verboseLevel() > 10) {
        cout << "CrossSection::dsigdtdxp_total_checked(): " << result;
        cout << " at t=" << t << ", xp=" << xpom;
        cout << " (" << (mDiffractiveMode == coherent ? "coherent" : "incoherent");
        cout << ')' << endl;
    }
    
    //
    // Check validity of return value
    //
    if (std::isnan(result)) {
        cout << "CrossSection::dsigdtdxp_total_checked(): Error, return value = NaN at" << endl;
        cout << "                                        t=" << t << ", xp=" << xpom << endl;
        result = 0;
    }
    if (std::isinf(result)) {
        cout << "CrossSection::dsigdtdxp_total_checked(): Error, return value = inf at" << endl;
        cout << "                                        t=" << t << ", xp=" << xpom << endl;
        result = 0;
    }
    if (result < 0) {
        cout << "CrossSection::dsigdtdxp_total_checked(): Error, negative cross-section at" << endl;
        cout << "                                        t=" << t << ", xp=" << xpom << endl;
        result = 0;
    }
    
    return result;
}

double CrossSection::dsigdt_total(double t, double Q2, double W2, GammaPolarization pol) const
{
    double result = 0;
    if (mSettings->tableSetType() == coherent_and_incoherent) {
        result = dsigdt_coherent(t, Q2, W2, pol) + dsigdt_incoherent(t, Q2, W2, pol);
    }
    else if (mSettings->tableSetType() == total_and_coherent) {
        result = mTableCollection->get(Q2,  W2,  t, pol, mean_A2);  // units now are fm^4
        result /= (16*M_PI);
        result /= hbarc2;                  // in fm^2/GeV^2
        result *= 1e7;                     // in nb/GeV^2
        
        double lambda = 0;
        if (mSettings->correctForRealAmplitude()){
            lambda = logDerivateOfAmplitude(t, Q2, W2, pol);
            result *= realAmplitudeCorrection(lambda);
        }
        if (mSettings->correctSkewedness()){
            lambda = logDerivateOfGluonDensity(t, Q2, W2, pol);
            result *= skewednessCorrection(lambda);
        }
    }
    return result;
}

//
//   UPC Version
//
double CrossSection::dsigdt_total(double t, double xpom) const
{
    double result = 0;
    if (mSettings->tableSetType() == coherent_and_incoherent) {
        result = dsigdt_coherent(t, xpom) + dsigdt_incoherent(t, xpom);
    }
    else if (mSettings->tableSetType() == total_and_coherent) {
        result = mTableCollection->get(xpom, t, mean_A2);  // units now are fm^4
        result /= (16*M_PI);
        result /= hbarc2;                  // in fm^2/GeV^2
        result *= 1e7;                     // in nb/GeV^2
        
        double lambda = 0;
        if (mSettings->correctForRealAmplitude()) {
            lambda = logDerivateOfAmplitude(t, xpom);
            result *= realAmplitudeCorrection(lambda);
        }
        if (mSettings->correctSkewedness()) {
            lambda = logDerivateOfGluonDensity(t, xpom);
            result *= skewednessCorrection(lambda);
        }
    }
    return result;
}

double CrossSection::dsigdt_coherent(double t, double Q2, double W2, GammaPolarization pol) const
{    
    double val = mTableCollection->get(Q2,  W2,  t, pol, mean_A);  // fm^2
    double result = val*val/(16*M_PI); // units now are fm^4
    result /= hbarc2;                  // in fm^2/GeV^2
    result *= 1e7;                     // in nb/GeV^2
    
    double lambda = 0;
    if (mSettings->correctForRealAmplitude()) {
        lambda = logDerivateOfAmplitude(t, Q2, W2, pol);
        result *= realAmplitudeCorrection(lambda);
    }
    if (mSettings->correctSkewedness()) {
        lambda = logDerivateOfGluonDensity(t, Q2, W2, pol);
        result *= skewednessCorrection(lambda);
    }
    
    return result;    
}    

//
//   UPC Version
//
double CrossSection::dsigdt_coherent(double t, double xpom) const
{
    double val = mTableCollection->get(xpom, t, mean_A);  // fm^2
    double result = val*val/(16*M_PI); // units now are fm^4
    result /= hbarc2;                  // in fm^2/GeV^2
    result *= 1e7;                     // in nb/GeV^2
    
    double lambda = 0;
    if (mSettings->correctForRealAmplitude()) {
        lambda = logDerivateOfAmplitude(t, xpom);
        result *= realAmplitudeCorrection(lambda);
    }
    if (mSettings->correctSkewedness()) {
        lambda = logDerivateOfGluonDensity(t, xpom);
        result *= skewednessCorrection(lambda);
    }
    
    return result;
}


double CrossSection::dsigdtdQ2dW2_total(double t, double Q2, double W2, GammaPolarization pol) const
{
    double result = dsigdt_total(t, Q2, W2, pol);
    if (mSettings->applyPhotonFlux()) result *= mPhotonFlux(Q2,W2,pol);
    
    return result;
}

//
//   UPC version
//
double CrossSection::dsigdtdxp_total(double t, double xpom) const
{
    double result = dsigdt_total(t, xpom);
    if (mSettings->applyPhotonFlux()) {
      result *= UPCPhotonFlux(t, xpom);
    }
    
    return result;
}

double CrossSection::dsigdt_incoherent(double t, double Q2, double W2, GammaPolarization pol) const
{
    double result = mTableCollection->get(Q2,  W2,  t, pol, variance_A);  // fm^4
    result /= (16*M_PI);
    result /= hbarc2;                  // in fm^2/GeV^2
    result *= 1e7;                     // in nb/GeV^2
    
    double lambda = 0;
    if (mSettings->correctForRealAmplitude()){ 
      lambda = logDerivateOfAmplitude(t, Q2, W2, pol);
      result *= realAmplitudeCorrection(lambda);
    }
    if (mSettings->correctSkewedness()){
      lambda = logDerivateOfGluonDensity(t, Q2, W2, pol);
      result *= skewednessCorrection(lambda);
    }

    return result;
}

//
//   UPC Version
//
double CrossSection::dsigdt_incoherent(double t, double xpom) const
{
    double result = mTableCollection->get(xpom, t, variance_A);  // fm^4
    result /= (16*M_PI);
    result /= hbarc2;                  // in fm^2/GeV^2
    result *= 1e7;                     // in nb/GeV^2
    
    double lambda = 0;
    if (mSettings->correctForRealAmplitude()) {
        lambda = logDerivateOfAmplitude(t, xpom);
        result *= realAmplitudeCorrection(lambda);
    }
    if (mSettings->correctSkewedness()) {
        lambda = logDerivateOfGluonDensity(t, xpom);
        result *= skewednessCorrection(lambda);
    }
    
    return result;
}

double CrossSection::dsigdtdQ2dW2_coherent(double t, double Q2, double W2, GammaPolarization pol) const    
{    
    double result = dsigdt_coherent(t, Q2, W2, pol);    
    if (mSettings->applyPhotonFlux())  result *= mPhotonFlux(Q2,W2,pol);    
        
    return result;    
}    

//
//   UPC Version
//
double CrossSection::dsigdtdxp_coherent(double t, double xpom) const
{
    double result = dsigdt_coherent(t, xpom);
    if (mSettings->applyPhotonFlux()) {
        result *= UPCPhotonFlux(t, xpom);
    }

    return result;
}

double CrossSection::UPCPhotonFlux(double t, double xpom) const{    
    double Egamma = Kinematics::Egamma(xpom, t, mVmMass,
			       mSettings->hadronBeamEnergy(), mSettings->electronBeamEnergy());
    double result = mPhotonFlux(Egamma);
    //
    // Jacobian = dEgamma/dN, such that dsig/dtdxp = dsigd/tdEgam*dEgam/dxp
    //
    double h = xpom*1e-3;
    double Egamma_p = Kinematics::Egamma(xpom+h, t, mVmMass,
				 mSettings->hadronBeamEnergy(), mSettings->electronBeamEnergy());
    double Egamma_m = Kinematics::Egamma(xpom-h, t, mVmMass,
				 mSettings->hadronBeamEnergy(), mSettings->electronBeamEnergy());
    double jacobian = (Egamma_p-Egamma_m)/(2*h);
    result *= fabs(jacobian);  
    return result;
}

GammaPolarization CrossSection::polarizationOfLastCall() const {return mPolarization;}    
    
DiffractiveMode CrossSection::diffractiveModeOfLastCall() const {return mDiffractiveMode;}

double CrossSection::crossSectionRatioLTOfLastCall() const {return mCrossSectionRatioLT;}


double CrossSection::logDerivateOfAmplitude(double t, double Q2, double W2, GammaPolarization pol) const
{
    double lambda = 0;
    bool   lambdaFromTable = true;
    Table *table = 0;
    
    if (!mProtonTableCollection) {
        cout << "CrossSection::logDerivateOfAmplitude(): no proton table defined to obtain lambda." << endl;
        cout << "                                        Corrections not available. Should be off." << endl;
        return 0;
    }
   
    //
    //  If the lambda table is present we use the more accurate and numerically 
    //  stable table value. Otherwise we calculate it from the <A> table(s).
    //
    //  Note (TU): if the lambda value from a table is not valid and not > 0,
    //             get() returns lambda=0 and table=0. This enforces a renewed
    //             calculation.
    //
    lambda = mProtonTableCollection->get(Q2,  W2,  t, pol, lambda_real, table);
    
    if (!table) {  // no lambda value from correct table, use fallback solution
        
        lambdaFromTable = false;
        
        double value = mProtonTableCollection->get(Q2,  W2,  t, pol, mean_A, table); // use obtained table from here on
        
        if (value <= 0) {
            cout << "CrossSection::logDerivateOfAmplitude(): got invalid value from table, value=" << value << '.' << endl;
            cout << "                                        t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2) << ", pol="
                 << (pol == transverse ? 'T' : 'L') << endl;
            return 0;
        }

        if (!table) {
            cout << "CrossSection::logDerivateOfAmplitude(): got invalid pointer to lookup table." << endl;
            return 0;
        }
        
        //
        //  Note: the derivate taken from values in the table is a delicate issue.
        //  The standard interpolation method(s) used in Table::get() are
        //  at times not accurate and causes ripples on lambda(W2).
        //  However, interpolation methods or robust fitting is by far too
        //  expensive in terms of CPU time per point.
        //
        
        //
        //  Calculate the derivative using simple method.
        //
        double derivative;
        double hplus, hminus;
        double dW2 = table->binWidthW2();
        hplus = hminus = 0.5*dW2; // half bin width is found to be the best choice after some testing
        hplus = min(hplus, table->maxW2()-W2);
        hminus = min(hminus, W2-table->minW2());
        hminus -= numeric_limits<float>::epsilon();
        hplus  -= numeric_limits<float>::epsilon();
        if (hminus < 0) hminus = 0;      
        if (hplus < 0) hplus = 0;    
        if (hplus == 0 && hminus == 0) {    
            cout << "CrossSection::logDerivateOfAmplitude(): Warning, cannot find derivative." << endl;    
            return 0;    
        }   
      
        double a =  table->get(Q2, W2+hplus, t);
        if (a <= 0) {
            cout << "CrossSection::logDerivateOfAmplitude(): got invalid value from table, value=" << a << '.' << endl;
            cout << "                                        t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2+hplus) << ", pol="
            << (pol == transverse ? 'T' : 'L') << endl;
            return 0;
        }
        double b  = table->get(Q2, W2-hminus, t);
        if (b <= 0) {
            cout << "CrossSection::logDerivateOfAmplitude(): got invalid value from table, value=" << b << '.' << endl;
            cout << "                                        t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2-hminus) << ", pol="
            << (pol == transverse ? 'T' : 'L') << endl;
            return 0;
        }
        derivative = (a-b)/(hplus+hminus);
        
        //
        //  Finally calculate lambda
        //
        double jacobian = (W2-protonMass2+Q2)/value;
        lambda = jacobian*derivative;
        
    } // end fall back solution

    if (mSettings->verboseLevel() > 3) {
        cout << "CrossSection::logDerivateOfAmplitude(): ";
        if (lambdaFromTable) 
            cout << "Info, lambda taken from table." << endl;
        else
            cout << "Info, lambda derived numerically from proton amplitude table" << endl;
    }
    
    //
    //  Check lambda value.
    //  At a lambda of ~0.6 both corrections have equal value
    //  of around 2.9. This will yield excessive large (unphysical)
    //  corrections. Large values are typically caused by fluctuations
    //  and glitches in the tables and should be rare.
    //

    double maxLambda = mSettings->maxLambdaUsedInCorrections();
    if (fabs(lambda) > maxLambda) {
        if (mSettings->verboseLevel() > 2) {
            cout << "CrossSection::logDerivateOfAmplitude(): ";
            cout << "Warning, lambda is excessively large (" << lambda << ") at " ;
            cout << "Q2=" << Q2 << ", W2=" << W2 << ", t=" << t << endl;
            cout << "Set to " << (lambda > 0 ? maxLambda : -maxLambda) << "." << endl;
        }
        lambda = lambda > 0 ? maxLambda : -maxLambda;
    }
    
    if (std::isinf(lambda)) {
        cout << "CrossSection::logDerivateOfAmplitude(): error, lambda = inf for pol=" << (pol == transverse ? 'T' : 'L') << endl;
        cout << "                                        t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2) << endl;
        cout << "                                        Set to 0." << endl;
        return 0;
    }
    
    if (std::isnan(lambda)) {
        cout << "CrossSection::logDerivateOfAmplitude(): error, lambda is NaN for pol=" << (pol == transverse ? 'T' : 'L') << endl;
        cout << "                                        t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2) << endl;
        cout << "                                        Set to 0." << endl;
        return 0;
    }
    return lambda;
} 

//
//   UPC only version
//
double CrossSection::logDerivateOfAmplitude(double t, double xpom) const
{
    double lambda = 0;
    bool   lambdaFromTable = true;
    Table *table = 0;
    
    if (!mProtonTableCollection) {
        cout << "CrossSection::logDerivateOfAmplitude(): no proton table defined to obtain lambda." << endl;
        cout << "                                        Corrections not available. Should be off." << endl;
        return 0;
    }
    
    //
    //   Usual numeric issues at boundaries.
    //   Subtracting an eps in log(xpom) does the trick.
    //
    if (xpom > mProtonTableCollection->maxX()) {
        xpom = exp(log(xpom)-numeric_limits<float>::epsilon());
    }
    if (xpom < mProtonTableCollection->minX()) {
        xpom = exp(log(xpom)+numeric_limits<float>::epsilon());
    }
    
    //
    //  If the lambda table is present we use the more accurate and numerically
    //  stable table value. Otherwise we calculate it from the <A> table(s).
    //
    lambda = mProtonTableCollection->get(xpom, t, lambda_real, table);
    
    if (!table) {  // no lambda value from correct table, use fallback solution
        
        lambdaFromTable = false;
        
        (void) mProtonTableCollection->get(xpom, t, mean_A, table); // use obtained table from here on
        
        if (!table) {
            cout << "CrossSection::logDerivateOfAmplitude(): got invalid pointer to lookup table." << endl;
            return 0;
        }
        
        //
        //  Here's the tricky part (see comments in non-UPC version above)
        //
        double theLogxpom = log(xpom);
        double dlogxpom = table->binWidthX();  // assuming table is in log x ???????? FIX later
        double maxLogxpom = log(table->maxX());
        double derivative;
        double hplus, hminus;
        hplus = hminus = 0.5*dlogxpom;
        hplus  = min(hplus, fabs(maxLogxpom-theLogxpom));
        hminus = min(hminus, fabs(theLogxpom-log(table->minX())));
        hminus -= numeric_limits<float>::epsilon();
        hplus  -= numeric_limits<float>::epsilon();
        if (hminus < 0) hminus = 0;
        if (hplus < 0) hplus = 0;
        if (hplus == 0 && hminus == 0) {
            cout << "CrossSection::logDerivateOfAmplitude(): Warning, cannot find derivative." << endl;
            return 0;
        }
        
        double a =  table->get(exp(theLogxpom+hplus), t);
        if (a <= 0) {
            cout << "CrossSection::logDerivateOfAmplitude(): got invalid value from table, value=" << a << '.' << endl;
            cout << "                                        t=" << t << ", W=" << sqrt(exp(theLogxpom+hplus)) << endl;
            return 0;
        }
        double b  = table->get(exp(theLogxpom-hminus), t);
        if (a <= 0) {
            cout << "CrossSection::logDerivateOfAmplitude(): got invalid value from table, value=" << b << '.' << endl;
            cout << "                                        t=" << t << ", W=" << sqrt(exp(theLogxpom-hminus)) << endl;
            return 0;
        }
        derivative = log(a/b)/(hplus+hminus);
        
        //
        //  Finally calculate lambda
        //  Directly dlog(A)/-dlog(xpom) here.
        //
        lambda = -derivative;
    }
    
    if (mSettings->verboseLevel() > 3) {
        cout << "CrossSection::logDerivateOfAmplitude(): ";
        if (lambdaFromTable)
        cout << "Info, lambda taken from table." << endl;
        else
        cout << "Info, lambda derived numerically from proton amplitude table" << endl;
        cout << "                                t=" << t << ", xpom=" << xpom << endl;
    }
    
    //
    //  Check lambda value.
    //  At a lambda of ~0.6 both corrections have equal value
    //  of around 2.9. This will yield excessive large (unphysical)
    //  corrections. Large values are typically caused by fluctuations
    //  and glitches in the tables and should be rare.
    //
    
    double maxLambda = mSettings->maxLambdaUsedInCorrections();
    if (fabs(lambda) > maxLambda) {
        if (mSettings->verboseLevel() > 2) {
            cout << "CrossSection::logDerivateOfAmplitude(): ";
            cout << "Warning, lambda is excessively large (" << lambda << ") at " ;
            cout << "xpom=" << xpom << ", t=" << t << endl;
            cout << "Set to " << (lambda > 0 ? maxLambda : -maxLambda) << "." << endl;
        }
        lambda = lambda > 0 ? maxLambda : -maxLambda;
    }
    
    if (std::isinf(lambda)) {
        cout << "CrossSection::logDerivateOfAmplitude(): error, lambda = infinity for"  << endl;
        cout << "                                        t=" << t << ", xpom=" << xpom << endl;
        cout << "                                        Set to 0." << endl;
        return 0;
    }
    
    if (std::isnan(lambda)) {
        cout << "CrossSection::logDerivateOfAmplitude(): error, lambda is NaN for" << endl;
        cout << "                                        t=" << t << ", xpom=" << xpom << endl;
        cout << "                                        Set to 0." << endl;
        return 0;
    }
    
    return lambda;
}

double CrossSection::logDerivateOfGluonDensity(double t, double Q2, double W2, GammaPolarization pol) const
{
    double lambda = 0;
    bool   lambdaFromTable = true;
    Table *table = 0;
    
    if (!mProtonTableCollection) {
        cout << "CrossSection::logDerivateOfGluonDensity(): no proton table defined to obtain lambda." << endl;
        cout << "                                           Corrections not available. Should be off." << endl;
        return 0;
    }
   
    //
    //  If the lambda table is present we use the more accurate and numerically 
    //  stable table value. Otherwise we calculate it from the <A> table(s).
    //    
    lambda = mProtonTableCollection->get(Q2,  W2,  t, pol, lambda_skew, table);
    
    if (!table) {
        //
        // no lambda value from correct table, use fallback solution
        // lambda_skew ~ lambda_real
        //
        lambda = logDerivateOfAmplitude(t, Q2, W2, pol) ; // This is an approximation in this case
    } // end fall back solution

    if (mSettings->verboseLevel() > 3) {
        cout << "CrossSection::logDerivateOfGluonDensity(): ";
        if (lambdaFromTable) 
            cout << "Info, lambda taken from table." << endl;
        else
            cout << "Info, lambda taken from logDerivateOfAmplitude as an approximation" << endl;
    }
    
    //
    //  Checking lambda.
    //  At a lambda of ~0.6 both corrections have equal value
    //  of around 2.9. This will yield excessive large (unphysical)
    //  corrections. Large values are typically caused by fluctuations
    //  and glitches in the tables and should be rare.
    //
    
    double maxLambda = mSettings->maxLambdaUsedInCorrections();
    if (fabs(lambda) > maxLambda) {
        if (mSettings->verboseLevel() > 2) {
            cout << "CrossSection::logDerivateOfGluonDensity(): ";
            cout << "Warning, lambda is excessively large (" << lambda << ") at " ;
            cout << "Q2=" << Q2 << ", W2=" << W2 << ", t=" << t << endl;
            cout << "Set to " << (lambda > 0 ? maxLambda : -maxLambda) << "." << endl;
        }
        lambda = lambda > 0 ? maxLambda : -maxLambda;
    }
    
    if (std::isinf(lambda)) {
        cout << "CrossSection::logDerivateOfGluonDensity(): error, lambda = inf for pol=" << (pol == transverse ? 'T' : 'L') << endl;
        cout << "                                        t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2) << endl;
        cout << "                                        Set to 0." << endl;
        return 0;
    }
    
    if (std::isnan(lambda)) {
        cout << "CrossSection::logDerivateOfGluonDensity(): error, lambda is NaN for pol=" << (pol == transverse ? 'T' : 'L') << endl;
        cout << "                                        t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2) << endl;
        cout << "                                        Set to 0." << endl;
        return 0;
    }
    
    return lambda;
} 

//
//   UPC version
//
double CrossSection::logDerivateOfGluonDensity(double t, double xpom) const
{
    double lambda = 0;
    bool   lambdaFromTable = true;
    Table *table = 0;
    
    if (!mProtonTableCollection) {
        cout << "CrossSection::logDerivateOfGluonDensity(): no proton table defined to obtain lambda." << endl;
        cout << "                                           Corrections not available. Should be off." << endl;
        return 0;
    }
    
    //
    //   Usual numeric issues at boundaries.
    //   Subtracting an eps in log(xpom) does the trick.
    //
    if (xpom > mProtonTableCollection->maxX()) {
        xpom = exp(log(xpom)-numeric_limits<float>::epsilon());
    }
    if (xpom < mProtonTableCollection->minX()) {
        xpom = exp(log(xpom)+numeric_limits<float>::epsilon());
    }
    
    //
    //  If the lambda table is present we use the more accurate and numerically
    //  stable table value. Otherwise we calculate it from the <A> table(s).
    //
    lambda = mProtonTableCollection->get(xpom, t, lambda_skew, table);
    
    if (!table) {
        //
        // no lambda value from correct table, use fallback solution
        // lambda_skew ~ lambda_real
        //
        
        lambdaFromTable = false;
        
        lambda = logDerivateOfAmplitude(t, xpom);
    }
    
    if (mSettings->verboseLevel() > 3) {
        cout << "CrossSection::logDerivateOfGluonDensity(): ";
        if (lambdaFromTable)
        cout << "Info, lambda taken from table." << endl;
        else
        cout << "Info, lambda taken from logDerivateOfAmplitude as approximation." << endl;
    }
    
    //
    //  Check lambda value.
    //  At a lambda of ~0.6 both corrections have equal value
    //  of around 2.9. This will yield excessive large (unphysical)
    //  corrections. Large values are typically caused by fluctuations
    //  and glitches in the tables and should be rare.
    //
    
    double maxLambda = mSettings->maxLambdaUsedInCorrections();
    if (fabs(lambda) > maxLambda) {
        if (mSettings->verboseLevel() > 2) {
            cout << "CrossSection::logDerivateOfGluonDensity(): ";
            cout << "Warning, lambda is excessively large (" << lambda << ") at " ;
            cout << "xpom=" << xpom << ", t=" << t << endl;
            cout << "Set to " << (lambda > 0 ? maxLambda : -maxLambda) << "." << endl;
        }
        lambda = lambda > 0 ? maxLambda : -maxLambda;
    }
    
    if (std::isinf(lambda)) {
        cout << "CrossSection::logDerivateOfGluonDensity(): error, lambda = infinity for"  << endl;
        cout << "                                        t=" << t << ", xpom=" << xpom << endl;
        cout << "                                        Set to 0." << endl;
        return 0;
    }
    
    if (std::isnan(lambda)) {
        cout << "CrossSection::logDerivateOfGluonDensity(): error, lambda is NaN for" << endl;
        cout << "                                        t=" << t << ", xpom=" << xpom << endl;
        cout << "                                        Set to 0." << endl;
        return 0;
    }
    
    return lambda;
}

double CrossSection::realAmplitudeCorrection(double lambda) const     
{          
    //    
    // Correction factor for real amplitude contribution    
    //    
    double beta = tan(lambda*M_PI/2.);    
    double correction = 1 + beta*beta;    
    
    return correction;    
}    

double CrossSection::realAmplitudeCorrection(double t, double Q2, double W2, GammaPolarization pol) const     
{          
    double lambda = logDerivateOfAmplitude(t, Q2, W2, pol);   
    double correction = realAmplitudeCorrection(lambda);    
    // correction *= exp(-10*Kinematics::x(Q2, W2));  // damped   
    return correction;    
}    

//
//    UPC version
//
double CrossSection::realAmplitudeCorrection(double t, double xpom) const
{
    double lambda = logDerivateOfAmplitude(t, xpom);
    double correction = realAmplitudeCorrection(lambda);
    return correction;
}

double CrossSection::skewednessCorrection(double lambda) const     
{    
    //    
    // Skewedness correction   
    //   
    double R = pow(2.,2*lambda+3)*TMath::Gamma(lambda+2.5)/(sqrt(M_PI)*TMath::Gamma(lambda+4));    
    double correction = R*R;    
    
    return correction;    
}    

double CrossSection::skewednessCorrection(double t, double Q2, double W2, GammaPolarization pol) const
{
    double lambda = logDerivateOfAmplitude(t, Q2, W2, pol);
    double correction = skewednessCorrection(lambda);
    // correction *= exp(-10*Kinematics::x(Q2, W2));  // damped
    return correction;
}

//
//    UPC version
//
double CrossSection::skewednessCorrection(double t, double xpom) const
{    
    double lambda = logDerivateOfAmplitude(t, xpom);
    double correction = skewednessCorrection(lambda);        
    return correction;    
}    

