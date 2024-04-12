//==============================================================================
//  DipoleModel.h
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
#ifndef DipoleModel_h  
#define DipoleModel_h  
#include "AlphaStrong.h"  
#include "TableGeneratorNucleus.h"  
#include "DipoleModelParameters.h"

class TH2F;  
class TH1F;  
  
class DipoleModel {  
public:  
    DipoleModel();  
    virtual ~DipoleModel();  
      
    const TableGeneratorNucleus* nucleus() const;  
    bool  configurationExists() const;  
      
    virtual void   createConfiguration(int)=0;  
    virtual double dsigmadb2(double, double, double, double)=0;  
    virtual double bDependence(double);
    virtual double bDependence(double, double);
    virtual double dsigmadb2ep(double, double, double);
    virtual double coherentDsigmadb2(double, double, double) {return 0;};
    virtual void   createSigma_ep_LookupTable(double) {/* nothing */};
    
protected:  
    TableGeneratorNucleus mNucleus;  
    DipoleModelParameters *mParameters;

    AlphaStrong mAs;
    bool        mConfigurationExists;  
    bool        mIsInitialized;
};
  
class DipoleModel_bSat : public DipoleModel {  
public:  
    DipoleModel_bSat();  
    DipoleModel_bSat(const DipoleModel_bSat&);
    DipoleModel_bSat& operator=(const DipoleModel_bSat&);  
    ~DipoleModel_bSat();

    void   createSigma_ep_LookupTable(double);
    
protected:  
    void   createConfiguration(int);    
    double dsigmadb2(double, double, double, double);  
    double bDependence(double, double);  
    double dsigmadb2ep(double, double, double);  
      
protected:  
    TH2F*  mBDependence;  
  
private:  
    double dsigmadb2epForIntegration(double*, double*);  
    TH1F*  mSigma_ep_LookupTable;  
    double coherentDsigmadb2(double, double, double);  

};  
  
class DipoleModel_bNonSat : public DipoleModel_bSat{  
public:  
    DipoleModel_bNonSat();  
    ~DipoleModel_bNonSat();  
    
private:  
    double dsigmadb2(double, double, double, double);  
    double dsigmadb2ep(double, double, double);  
    double coherentDsigmadb2(double, double, double);  
};  

class DipoleModel_bCGC : public DipoleModel {  
public:
    DipoleModel_bCGC();
    
private:
    void   createConfiguration(int);    
    double dsigmadb2(double, double, double, double);  
    double dsigmadb2ep(double, double, double);  
    double bDependence(double);  
};  

#endif  
