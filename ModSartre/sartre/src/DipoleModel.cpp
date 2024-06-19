//==============================================================================
//  DipoleModel.cpp
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
//  Author: Tobias Toll
//  Last update: 
//  $Date: 2019-04-01 14:39:33 -0400 (Mon, 01 Apr 2019) $
//  $Author: ullrich $
//==============================================================================
#include <fstream>  
#include <iostream>  
#include <sstream>  
#include <cmath>  
#include "DipoleModel.h"  
#include "TableGeneratorSettings.h"  
#include "DglapEvolution.h"  
#include "Constants.h"  
#include "TFile.h"  
#include "TVector3.h"  
#include "TMath.h"  
#include "TH2F.h"  
#include "TF1.h"  
  
#define PR(x) cout << #x << " = " << (x) << endl;  
  
using namespace std;  
  
DipoleModel::DipoleModel()  
{  
    mConfigurationExists=false;      
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();  
    unsigned int A = settings->A();  
    mNucleus.init(A);  
    mIsInitialized=true;  
    mParameters = 0;
}  
  
DipoleModel::~DipoleModel()
{
    delete mParameters;
}
  
const TableGeneratorNucleus* DipoleModel::nucleus() const { return &mNucleus; }  
  
bool DipoleModel::configurationExists() const { return mConfigurationExists; }  
  
double DipoleModel::bDependence(double) { return 0; }  

double DipoleModel::bDependence(double, double) { return 0; }

double DipoleModel::dsigmadb2ep(double, double, double) { return 0;}  
  
  
//***********bSat:*****************************************************  
DipoleModel_bSat::DipoleModel_bSat()  
{  
    mBDependence = 0;  
    mSigma_ep_LookupTable = 0;
    //
    //  Set the parameters. Note that we enforce here the bSat model
    //  independent of what the settings say.
    //
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    mParameters = new DipoleModelParameters(bSat, settings->dipoleModelParameterSet());

}  
  
DipoleModel_bSat::~DipoleModel_bSat()  
{  
    delete mSigma_ep_LookupTable;  
    delete mBDependence;
}  
  
DipoleModel_bSat& DipoleModel_bSat::operator=(const DipoleModel_bSat& dp)  
{  
    if (this != &dp) {  
        delete mBDependence;  
        delete mSigma_ep_LookupTable;  
	
        DipoleModel::operator=(dp);  
        mBDependence = new TH2F(*(dp.mBDependence));  
        mBDependence->SetDirectory(0);  
        mSigma_ep_LookupTable = new TH1F(*(dp.mSigma_ep_LookupTable));  
        mSigma_ep_LookupTable->SetDirectory(0);
	  
    }  
    return *this;      
}  
  
DipoleModel_bSat::DipoleModel_bSat(const DipoleModel_bSat& dp) : DipoleModel(dp)  
{  
    mBDependence = new TH2F(*(dp.mBDependence));  
    mBDependence->SetDirectory(0);
}  
  
void DipoleModel_bSat::createConfiguration(int iConfiguration)  
{  
    if (!mIsInitialized) {  
        cout << "DipoleModel_bSat::createConfiguration(): DipoleModel class has not been initialized! Stopping." << endl;  
        exit(1);  
    }  
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();  
    unsigned int A = mNucleus.A();  
    string path=settings->bSatLookupPath();  
    ostringstream filename;  
    filename.str("");  
    filename << path << "/bSat_bDependence_A" << A <<".root";  
    ifstream ifs(filename.str().c_str());  
    if (!ifs) {  
        cout << "DipoleModel_bSat::createConfiguration(): File does not exist: " << filename.str().c_str() << endl;  
        cout << "Stopping." << endl;  
        exit(1);  
    }  
    TFile* lufile= new TFile(filename.str().c_str());  
    ostringstream histoName;  
    histoName.str( "" );  
    histoName << "Configuration_" << iConfiguration;  
    lufile->GetObject( histoName.str().c_str(), mBDependence );  
    mBDependence->SetDirectory(0);  
    lufile->Close();  
    mConfigurationExists=true;  
}  
  
double DipoleModel_bSat::dsigmadb2(double r, double b, double phi, double xprobe)  
{  
    double bDep=bDependence(b, phi);  
    double muQ2 = mParameters->C()/(r*r/hbarc2) + mParameters->mu02();
    double asxg = DglapEvolution::instance().alphaSxG(xprobe, muQ2);
    double omega = ((M_PI*M_PI)/Nc)*(r*r/hbarc2)*asxg*bDep;
    double result = 2.*(1. - exp(-omega/2));  

    return result;  
}  
  
double DipoleModel_bSat::bDependence(double b, double phi)  
{  
    double result = mBDependence->Interpolate(b, phi);  
    return result;  
}  
  
double DipoleModel_bSat::dsigmadb2ep(double r, double b, double xprobe)  
{  
    const double BG = mParameters->BG(); // GeV^-2  
    double arg = b*b/(2*BG);  
    arg /= hbarc2;  
    double bDep= 1/(2*M_PI*BG) * exp(-arg);  
    double Mu02 = mParameters->mu02(); // GeV^2  
    double muQ2 = mParameters->C()/(r*r/hbarc2) + Mu02;
    double asxg = DglapEvolution::instance().alphaSxG(xprobe, muQ2);
    double omega = ((M_PI*M_PI)/Nc)*(r*r/hbarc2)*asxg*bDep;
    double result = 2.*(1. - exp(-omega/2));  

    return result;  
}  

double DipoleModel_bSat::coherentDsigmadb2(double r,  double b, double xprobe){  
    xprobe*=1.;  
    double sigmap = mSigma_ep_LookupTable->Interpolate(r);  
    int A=nucleus()->A();  
    double TA=nucleus()->T(b)/A;   
    double result = 2 * ( 1 - pow(1 - TA/2.*sigmap, A) );  
    return result;  
}  
  
void DipoleModel_bSat::createSigma_ep_LookupTable(double xprobe)  
{  
    double rbRange=3.*nucleus()->radius();  
    TF1* dsigmaForIntegration = new TF1("dsigmaForIntegration", this, &DipoleModel_bSat::dsigmadb2epForIntegration, 0., rbRange, 2);  
    mSigma_ep_LookupTable = new TH1F("", "", 1000, 0, rbRange);  
    dsigmaForIntegration->SetNpx(1000);  
    for (int iR=1; iR<=1000; iR++) {  
        double r=mSigma_ep_LookupTable->GetBinCenter(iR);  
        dsigmaForIntegration->SetParameter(0, r);  
        dsigmaForIntegration->SetParameter(1, xprobe);  
        double result=dsigmaForIntegration->Integral(0, rbRange);  
        mSigma_ep_LookupTable->SetBinContent(iR, result);  
    }  
    delete dsigmaForIntegration;  
}  
  
double DipoleModel_bSat::dsigmadb2epForIntegration(double *x, double* par)  
{  
    double r=par[0];  
    double xprobe=par[1];  
    double b = *x;  
    return 2*M_PI*b/hbarc2*dsigmadb2ep(r, b, xprobe);  
}  
  
  
//***********bNonSat:*************************************************  
DipoleModel_bNonSat::DipoleModel_bNonSat()
{
    //
    //  Set the parameters. Note that we need bNonSat to calculate
    //  the skewedness correction for bSat.   
    //
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    mParameters = new DipoleModelParameters(settings->dipoleModelType(), settings->dipoleModelParameterSet());
}

DipoleModel_bNonSat::~DipoleModel_bNonSat(){}  
  
  
double DipoleModel_bNonSat::dsigmadb2ep(double r, double b, double xprobe)  
{  
    const double BG = mParameters->BG(); // GeV^-2  
    double arg = b*b/(2*BG);  
    arg /= hbarc2;  
    double bDep= 1/(2*M_PI*BG) * exp(-arg);  
    double Mu02 = mParameters->mu02(); // GeV^2  
    double muQ2 = mParameters->C()/(r*r/hbarc2) + Mu02;
    double asxg = DglapEvolution::instance().alphaSxG(xprobe, muQ2);
    double omega = ((M_PI*M_PI)/Nc)*(r*r/hbarc2)*asxg*bDep;
      
    return omega;  
}  
  
double DipoleModel_bNonSat::dsigmadb2(double r, double b, double phi, double xprobe)  
{  
    double bDep=bDependence(b, phi);  
    double muQ2 = mParameters->C()/(r*r/hbarc2) + mParameters->mu02();
    double asxg = DglapEvolution::instance().alphaSxG(xprobe, muQ2);
    double omega = ((M_PI*M_PI)/Nc)*(r*r/hbarc2)*asxg*bDep;
  
    return omega;  
}  
  
double DipoleModel_bNonSat::coherentDsigmadb2(double r, double b, double xprobe){  
    int A=nucleus()->A();  
    double TA=nucleus()->T(b)/A;   
    double muQ2 = mParameters->C()/(r*r/hbarc2) + mParameters->mu02();  
    double asxg = DglapEvolution::instance().alphaSxG(xprobe, muQ2);
    double result=A*TA*M_PI*M_PI/Nc*r*r/hbarc2*asxg;  
    return result;  
}  
  
  
//***********bCGC:*****************************************************  

DipoleModel_bCGC::DipoleModel_bCGC()
{
    //
    //  Set the parameters. Note that we enforce here the bNonSat model
    //  independent of what the settings say.
    //
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    mParameters = new DipoleModelParameters(bCGC, settings->dipoleModelParameterSet());
}


void DipoleModel_bCGC::createConfiguration(int iConfiguration)
{  
    if (!mIsInitialized) {  
        cout << "DipoleModel_bCGC::createConfigurationDipoleModel class has not been initialized! Stopping." << endl;  
        exit(1);  
    }  
    //iConfiguration is a dummy for bCGC. Do this to avoid compiler warnings:  
    iConfiguration++;   
    mNucleus.generate();    
    mConfigurationExists=true;  
}  
  
double DipoleModel_bCGC::dsigmadb2(double r, double b, double phi, double x)  
{  
    double result=1;  
    for (unsigned int iA=0; iA<mNucleus.A(); iA++) {  
        double absdeltab=( TVector3(b*cos(phi), b*sin(phi), 0.)  
                          -mNucleus.configuration.at(iA).position() ).Perp();  
        result*=(1.-0.5*dsigmadb2ep(r, absdeltab, x));  
    }  
    return 2.*(1.-result);  
}  
  
double DipoleModel_bCGC::dsigmadb2ep(double r, double b, double xprobe)  
{  
    double Y = log(1/xprobe);  
    double kappa = mParameters->kappa();
    double N0 = mParameters->N0();
    double x0 = mParameters->x0();
    double lambda = mParameters->lambda();
    double gammas = mParameters->gammas();
    double A = -N0*N0*gammas*gammas/((1-N0)*(1-N0)*log(1-N0));  
    double B = 0.5*pow(1-N0,-(1-N0)/(N0*gammas));  
    double Qs = pow(x0/xprobe,lambda/2)*sqrt(DipoleModel_bCGC::bDependence(b));  
    double rQs = r*Qs/hbarc;  
    double result=0;  
    if (rQs <= 2)   
        result = 2*N0*pow(0.5*rQs, 2*(gammas+(1/(kappa*lambda*Y))*log(2/rQs)));  
    else   
        result = 2*(1 - exp(-A*log(B*rQs)*log(B*rQs)));  
    return result;  
}  
  
double DipoleModel_bCGC::bDependence(double b)  
{  
    double gammas = mParameters->gammas();
    double Bcgc = mParameters->Bcgc();
    return exp(-0.5*b*b/Bcgc/gammas/hbarc2);  
}  
  
