 //==============================================================================
//  CrossSection.h
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
//         
//  Functor class.        
//         
//  operator() returns d^3sig/(dt dQ2 dW2) in nb/GeV^6.         
//         
//===============================================================================         
#ifndef CrossSection_h         
#define CrossSection_h         
#include "Enumerations.h"         
#include "PhotonFlux.h"         
         
class TRandom3;         
class EventGeneratorSettings;         
class TableCollection;         
         
class CrossSection {         
public:         
    CrossSection(TableCollection* = 0, TableCollection* = 0);         
    ~CrossSection();         
             
    double operator()(double t, double Q2, double W2);           
    double operator()(double t, double xpom);  // UPC version
    double operator()(const double*);          // array of t, Q2, W2 or t, xpom for UPC
             
    double unuranPDF(const double*);          // for UNU.RAN using log(Q2) (or log(xpom) for UPC)
                                              // and returning log of cross-section
    
    void setTableCollection(TableCollection*);         
    void setProtonTableCollection(TableCollection*);         
    GammaPolarization polarizationOfLastCall() const;         
    DiffractiveMode diffractiveModeOfLastCall() const;         
    double crossSectionRatioLTOfLastCall() const;

    void setCheckKinematics(bool);  
    double dsigdtdQ2dW2_total(double t, double Q2, double W2, GammaPolarization) const;

protected:
    friend class Sartre;  // mostly for debugging and QA
    
    double dsigdtdQ2dW2_total_checked(double t, double Q2, double W2);
    double dsigdtdxp_total_checked(double t, double xpom);

    double dsigdt_total(double t, double Q2, double W2, GammaPolarization) const;  // modified
    double dsigdt_coherent(double t, double Q2, double W2, GammaPolarization) const;         
    double dsigdt_incoherent(double t, double Q2, double W2, GammaPolarization) const;  // new
    double dsigdtdQ2dW2_coherent(double t, double Q2, double W2, GammaPolarization) const;
    
    double dsigdt_total(double t, double xpom) const;  // UPC version
    double dsigdt_coherent(double t, double xpom) const;  // UPC version
    double dsigdt_incoherent(double t, double xpom) const;  // UPC version
    double dsigdtdxp_total(double t, double xpom) const;  // UPC
    double dsigdtdxp_coherent(double t, double xpom) const;  // UPC

    double logDerivateOfAmplitude(double t, double Q2, double W2, GammaPolarization) const;  
    double logDerivateOfGluonDensity(double t, double Q2, double W2, GammaPolarization) const;  
    double realAmplitudeCorrection(double t, double Q2, double W2, GammaPolarization pol) const;        
    double skewednessCorrection(double t, double Q2, double W2, GammaPolarization pol) const;

    double logDerivateOfAmplitude(double t, double xpom) const; // UPC version
    double logDerivateOfGluonDensity(double t, double xpom) const; // UPC version
    double realAmplitudeCorrection(double t, double xpom) const; // UPC version
    double skewednessCorrection(double t, double xpom) const; // UPC version

    double skewednessCorrection(double lambda) const;
    double realAmplitudeCorrection(double lambda) const;

    double UPCPhotonFlux(double t, double xpom) const;
    
private:         
    TRandom3*          mRandom;         
    GammaPolarization  mPolarization;         
    DiffractiveMode    mDiffractiveMode;         
    PhotonFlux         mPhotonFlux;         
    TableCollection*   mTableCollection;         
    TableCollection*   mProtonTableCollection;         
    EventGeneratorSettings* mSettings;         
             
    double mS;         
    double mVmMass;
    double mCrossSectionRatioLT;
    
    bool mCheckKinematics;  
};         
#endif         
