//==============================================================================
//  PhotonFlux.cpp
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
//  Author: Thomas Ullrich, Tobias Toll
//  Last update: 
//  $Date: 2021-07-02 14:28:12 -0400 (Fri, 02 Jul 2021) $
//  $Author: ullrich $
//==============================================================================
#include "PhotonFlux.h"    
#include "Constants.h"    
#include "Kinematics.h"    
#include "EventGeneratorSettings.h"  
#include "Nucleus.h"
#include <cmath>    
#include "TF1.h"
#include "Math/Functor.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

#define PR(x) cout << #x << " = " << (x) << endl;    

PhotonFlux::PhotonFlux()
{
    mS = 0;
    mIsUPC = false;
    mSettings = EventGeneratorSettings::instance();
    
    if (mSettings->UPC()){
        mNucleusUPC = new Nucleus(mSettings->UPCA());
        mNucleus = new Nucleus(mSettings->A());
        if (mSettings->UPCA()>1 && mSettings->A()>1) {
            cout << "PhotonFlux::PhotonFlux(): Calculating TAA table ... ";
            setupTAAlookupTable();
            cout << "done."  << endl;
        }
        mEBeamEnergy = mSettings->electronBeamEnergy(); // the beam that shines
        mIsUPC = true;
    }
    
    mSIsSet = false;
}    

PhotonFlux::PhotonFlux(double s)
{
    mS = s;
    mIsUPC = false;
    mSettings = EventGeneratorSettings::instance();
    
    if (mSettings->UPC()) {
        mNucleusUPC = new Nucleus(mSettings->UPCA());
        mNucleus = new Nucleus(mSettings->A());
        if (mSettings->UPCA()>1 && mSettings->A()>1){
            cout << "PhotonFlux::PhotonFlux(): Calculating TAA table ... ";
            setupTAAlookupTable();
            cout << "done" << endl;
        }
        mEBeamEnergy = mSettings->electronBeamEnergy(); //The beam that shines
        mIsUPC = true;
        cout<<"PhotonFlux::PhotonFlux(): Warning Sartre is not yet ready for UPC..." << endl;
    }
    
    mSIsSet = true;
}    

void PhotonFlux::setS(double s)
{
    mS = s;
    mSIsSet = true;
}

double PhotonFlux::operator()(double Q2, double W2, GammaPolarization p) const  
{
    if (!mSIsSet) cout<<"Warning in PhotonFlux::operator() s is not set!"<<endl;
    if (p == transverse)
        return fluxTransverse(Q2, W2);
    else
        return fluxLongitudinal(Q2, W2);
}

double PhotonFlux::operator()(double Egamma) const
{
    if (!mSIsSet) cout << "PhotonFlux::operator(): Warning, s is not set!" << endl;
    return nuclearPhotonFlux(Egamma);
}

double PhotonFlux::fluxTransverse(double Q2, double W2) const    
{    
    //
    //  Transverse photon flux
    //
    double y = Kinematics::y(Q2, Kinematics::x(Q2, W2), mS);
    
    if (y<0 || y>1 || Kinematics::error()) return 0;
    if (Q2 > Kinematics::Q2max(mS)) return 0;
    
    double result = alpha_em/(2*M_PI*Q2*mS*y);
    result *= (1 + (1-y)*(1-y)) - (2*(1-y)*Kinematics::Q2min(y)/Q2);
    
    return result;
}    

double PhotonFlux::fluxLongitudinal(double Q2, double W2) const    
{    
    //
    //  Longitudinal photon flux
    //
    double y = Kinematics::y(Q2, Kinematics::x(Q2, W2), mS);
    
    if (y<0 || y>1 || Kinematics::error()) return 0;
    if (Q2 > Kinematics::Q2max(mS)) return 0;
    
    double result = alpha_em/(2*M_PI*Q2*mS*y);
    result *= 2*(1-y);
    
    return result;
}    


double PhotonFlux::nuclearPhotonFlux(double Egamma) const
{
    if (Egamma < 0) {
        if (mSettings->verbose() && mSettings->verboseLevel() > 4)
            cout << "PhotonFlux::nuclearPhotonFlux(): Warning, negative Egamma as argument.\n"
                 << "                                 Return 0 for photon flux." << endl;
        return 0;
    }
    double bmin = 0.;
    double bmax = 1e10; //Could be arbitrary large (TF1 doesn't care)
    TF1 fluxFunction("fluxFunction", this, &PhotonFlux::unintegratedNuclearPhotonFlux, bmin, bmax, 1);
    fluxFunction.SetParameter(0, Egamma);
    ROOT::Math::WrappedTF1 wrappedFluxFunction(fluxFunction);
    ROOT::Math::GaussIntegrator fluxIntegrator;
    fluxIntegrator.SetFunction(wrappedFluxFunction);
    fluxIntegrator.SetAbsTolerance(0.);
    fluxIntegrator.SetRelTolerance(1e-5);
    
    return fluxIntegrator.IntegralUp(bmin)/hbarc;   //dN/dEgamma GeV^-1
}

double PhotonFlux::unintegratedNuclearPhotonFlux(double* var, double* par) const
{
    //
    // https://arxiv.org/abs/1607.03838v1
    //
    double b = var[0]; //fm
    double eGamma = par[0]; //GeV
    
    int Apom = mNucleus->A();
    int AUPC = mNucleusUPC->A();
    
    double Z = mNucleusUPC->Z();
    double ebeamMass = mNucleusUPC->nuclearMass()/mNucleusUPC->A(); //Mass per nucleon //GeV
    double beamLorentzGamma = mEBeamEnergy/ebeamMass;
    double beamLorentzGamma2 = beamLorentzGamma*beamLorentzGamma;
    
    double xi = eGamma*b/hbarc/beamLorentzGamma; //GeV0
    
    double K0 = TMath::BesselK0(xi);
    double K1 = TMath::BesselK1(xi);
    
    double prefactor = Z*Z*alpha_em/(M_PI*M_PI)*eGamma/beamLorentzGamma2; //GeV1
    double N = prefactor*(K1*K1+K0*K0/beamLorentzGamma2); //GeV1
    
    //Spectrum is calculated under the assumption that there is no
    //hadronic interaction between the beam particles.
    //Therefore it has to be multiplied by the probability for this:
    double Pnohad = 1;
    if (Apom>1 && AUPC>1) { //nucleus-nucleus interaction
        Pnohad = exp(-sigma_nn(mS)*mTAA_of_b->Interpolate(b)); //GeV0
    }
    else if (Apom == 1 && AUPC > 1) { //proton-nucleus
        Pnohad=exp(-sigma_nn(mS)*mNucleusUPC->T(b));
    }
    else if (Apom > 1 && AUPC == 1) { //nucleus-proton
        Pnohad = exp(-sigma_nn(mS)*mNucleus->T(b));
    }
    else { //proton-proton
        double b02 = 19.8; //GeV^-2
        double scamp = 1-exp(-b*b/hbarc2/2./b02);
        Pnohad = scamp*scamp;
    }
    
    double result = 2*M_PI*b/hbarc*N*Pnohad; //GeV0
    return result;
}


void PhotonFlux::setupTAAlookupTable()
{
    double rbhigh = upperIntegrationLimit*(mNucleus->radius()+mNucleusUPC->radius());
    mTAA_of_b = new TH1D("mTAA_of_b", "mTAA_of_b", 1000, 0, rbhigh);
    for (unsigned int i=1; i<=1000; i++) {
        double b = mTAA_of_b->GetBinCenter(i);
        mTAA_of_b->SetBinContent(i, TAA(b));
    }
}

double PhotonFlux::sigma_nn(double s) const {
    //
    // For parameters see
    // http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20080014212.pdf table 2
    //
    // Following KNSGB: http://arxiv.org/abs/1607.03838v1 (eq.6)
    //
    double Z = 33.73; //mb
    double B = 0.2838; //mb
    double s0 = 1.;//GeV2
    double Y1 = 13.67; //mb
    double s1 = 1.; //GeV2
    double eta1 = 0.412;
    double Y2 = 7.77; //mb
    double eta2 = 0.5626;
    
    double result = Z+B*log(s/s0)*log(s/s0)+Y1*pow(s1/s, eta1)-Y2*pow(s1/s, eta2); //mb = 1e-3b = 1e-1 fm2
    
    result *= 1e-1; //fm2 (1 barn = 1e2 fm2)
    
    return result/hbarc2; //GeV^-2
}

double PhotonFlux::TAA(double b)
{
    //
    // Integral over dr^2=2*pi*r*dr
    //
    
    double rbhigh=upperIntegrationLimit*(mNucleus->radius()+mNucleusUPC->radius());
    
    double lo[2] = {0., 0.}; //r, theta
    double hi[2] = {rbhigh, 2*M_PI};
    ROOT::Math::Functor wf(this, &PhotonFlux::TAAForIntegration, 2);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    ig.SetFunction(wf);
    ig.SetAbsTolerance(0.);
    ig.SetRelTolerance(1e-4);
    mB=b;
    
    double result = ig.Integral(lo, hi); //GeV^3*fm
    result /= hbarc; //GeV^2
    return result;
}

double PhotonFlux::TAAForIntegration(const double* var) const
{
    double r = var[0]; //fm
    double phi = var[1];
    double b = mB; //fm
    
    //Put A1 in the origin, and A2 on the y-axis at distance b:
    double arg1 = r;
    double arg2 = sqrt(r*r+b*b-2*b*r*sin(phi)); //fm
    double T1 = mNucleusUPC->T(arg1); //GeV^2
    double T2 = mNucleus->T(arg2); //GeV^2
    
    return r/hbarc*T1*T2; //GeV^3
}


