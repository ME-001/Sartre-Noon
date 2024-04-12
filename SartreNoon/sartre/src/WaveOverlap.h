//==============================================================================
//  WaveOverlap.h
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
//  $Date: 2019-03-08 14:12:33 -0500 (Fri, 08 Mar 2019) $
//  $Author: ullrich $
//==============================================================================
#ifndef WaveOverlap_h  
#define WaveOverlap_h  

class DipoleModelParameters;

class WaveOverlap {  
public:  
    WaveOverlap();  
    virtual ~WaveOverlap();  
    virtual void setWaveOverlapFunctionParameters(int);  
    virtual double T(double, double, double)=0;  
    virtual double L(double, double, double)=0;  
    virtual void setProcess(int);
    virtual void testBoostedGaussianParameters(int);

protected:
    DipoleModelParameters *mParameters;
};
  
class WaveOverlapVM : public WaveOverlap {  
public:  
    WaveOverlapVM();  
    void setWaveOverlapFunctionParameters(int);  
    void setProcess(int);  
    void testBoostedGaussianParameters(int);
    
private:  
    double T(double, double, double);  
    double L(double, double, double);  
    double transverseWaveFunction(double, double);  
    double longitudinalWaveFunction(double, double);  
    double dDrTransverseWaveFunction(double, double);  
    double laplaceRLongitudinalWaveFunction(double, double);  
    double uiDecayWidth(double*, double*);
    double uiNormL(const double*);
    double uiNormT(const double*);
    
private:  
    double mNT, mRT2;  
    double mMf;  // mass of quarks in vector meson
    double mBoostedGaussianMf; // // mass of quarks in vector meson's boosted Gaussian wave fct.
    double mMf2;  
    double mBoostedGaussianMf2;
    double mEf;
    double mMV;  
    double mNL, mRL2;
    
};
  
class WaveOverlapDVCS : public WaveOverlap {  
private:  
    double T(double, double, double);  
    double L(double, double, double);  
};  
  
inline void WaveOverlap::setWaveOverlapFunctionParameters(int) {/* no op*/}  
inline void WaveOverlap::setProcess(int) {/* no op*/};
inline void WaveOverlap::testBoostedGaussianParameters(int) {/* no op*/};

#endif  
