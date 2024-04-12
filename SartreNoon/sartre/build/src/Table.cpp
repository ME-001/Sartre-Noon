//==============================================================================
//  Table.cpp
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
//  $Date: 2019-05-30 16:57:30 -0400 (Thu, 30 May 2019) $
//  $Author: ullrich $
//==============================================================================
//    
//  Table data is stored in a TH3F.    
//    
//  All information is stored in the histogram title.    
//    
//  Table ID is encoded as follows:    
//    
//  histogram title is a number that is to be interpreted     
//  as an uint64_t with the bits set as follows:    
//    
//  bit 0:     content type: 0 for <A>, 1 for <A2>  (see also bits 33 and 45)
//  bit 1:     polarization: T=1, L=0    
//  bit 2:     t encoding: 0 for |t|, 1 for log(|t|)    
//  bit 3:     W2 encoding: 0 for lin, 1 for log    
//  bit 4:     Q2 encoding: 0 for lin, 1 for log     
//  bit 5-7:   dipole model type (Enumerations.h)    
//  bit 8-15:  mass number A    
//  bit 16-31: vector meson ID (PDG)    
//  bit 32:    content encoding: 0 in lin, 1 in log   
//  bit 33:    content type is lambda_<A> (bits 0 and 45 are zero in this case)
//  bit 34-41: priority
//  bit 42-44: dipole model parameter set (Enumerations.h)
//  bit 45:    content type is variance_A (bit 33 and 45 are zero in this case)
//  bit 46:    upc only table (tables with one Q2 bin only)
//  bit 47:    xpomeron encoding: 0 for lin, 1 for log (UPC only)
//  bit 48:    content type is lambda_skewedness (bits 0, 33, and 45 are zero in this case).
//  bit 49-63: not used
//     
//  log implies ln    
//  bit 0 is the LSB    
//    
//  Internal histogram: x = Q2, y = W2, z = t (or the logs)
//                      x = xpomeron, y = t (or the logs)     // UPC only
//    
//  Actual value is taken at the center of the bin.    
//  min/max functions (e.g. minQ2()) return the value of the first/last    
//  bin center.    
//==============================================================================   
#include "Table.h"    
#include "TFile.h"    
#include "TH2F.h"
#include "TH3F.h"
#include <cstdio>
#include <cmath>    
#include <ctime>   
#include <algorithm>    
#include <string>    
#include <sstream>    
#include <fstream>    
#include <iomanip>    
#include <limits>    
#include <unistd.h>   
#include <sys/types.h>   
#include "GridSpline.h"   
    
#define PR(x) cout << #x << " = " << (x) << endl;    
    
Table::Table()     
{    
    mID = 0;    
    m3DHisto = 0;    
    m2DHisto = 0;
    mFillCounter = 0;
    mBackupFrequence = 0;  
}    
    
Table& Table::operator=(const Table& tab)  
{  
    if (this != &tab) {  
        delete m3DHisto;  
        delete m2DHisto;

        mID = tab.mID;         
        if (tab.m3DHisto) {
            m3DHisto = new TH3F(*(tab.m3DHisto));
            m3DHisto->SetDirectory(0);
        }
        if (tab.m2DHisto) {
            m2DHisto = new TH2F(*(tab.m2DHisto));
            m2DHisto->SetDirectory(0);
        }
        mFilename = tab.mID;
        mFillCounter = tab.mID;
        mBackupFrequence = tab.mID;
        mBackupPrefix = tab.mID;
        mLastBackupFilename = tab.mID;
    }  
    return *this;  
}  
  
Table::Table(const Table& tab)  
{  
    mID = tab.mID;         
    if (tab.m3DHisto) {
        m3DHisto = new TH3F(*(tab.m3DHisto));
        m3DHisto->SetDirectory(0);
    }
    if (tab.m2DHisto) {
        m2DHisto = new TH2F(*(tab.m2DHisto));
        m2DHisto->SetDirectory(0);
    }
    mFilename = tab.mID;
    mFillCounter = tab.mID;
    mBackupFrequence = tab.mID;
    mBackupPrefix = tab.mID;
    mLastBackupFilename = tab.mID;
}  
  
Table::~Table()     
{    
    delete m2DHisto;
    delete m3DHisto;
}
    
unsigned int Table::create(int nbinsQ2, double Q2min, double Q2max,    
                           int nbinsW2, double W2min, double W2max,    
                           int nbinsT, double tmin, double tmax,    
                           bool logQ2, bool logW2, bool logt, bool logC,    
                           AmplitudeMoment mom, GammaPolarization pol,     
                           unsigned int A, int vm,    
                           DipoleModelType model, DipoleModelParameterSet pset,
                           const char* filename, unsigned char priority)
{    
    if (m3DHisto != 0) {    
        cout << "Table:create(): cannot create, table already exists." << endl;    
        return 0;    
    }    
    
    if (nbinsQ2 < 2) {
        cout << "Table:create(): number of bins in Q2 must be at least 2." << endl;
        return 0;
    }
    if (nbinsW2 < 2) {
        cout << "Table:create(): number of bins in W2 must be at least 2." << endl;
        return 0;
    }
    if (nbinsT < 2) {
        cout << "Table:create(): number of bins in t must be at least 2." << endl;
        return 0;
    }
    
    if (filename) mFilename = string(filename); // name of file where table is written to
      
    //    
    //  Create table ID first    
    //    
    mID = (vm << 16);    
    mID |= (A << 8);    
    mID |= (model << 5);    
    mID |= (static_cast<uint64_t>(pset) << 42);
    if (mom == mean_A2) mID |= 1;
    if (pol == transverse) mID |= (1 << 1);    
    if (logt)        mID |= (1 << 2);    
    if (logW2)       mID |= (1 << 3);    
    if (logQ2)       mID |= (1 << 4);    
    if (logC)        mID |= (static_cast<uint64_t>(1) << 32);   
    if (mom == lambda_real) mID |= (static_cast<uint64_t>(1) << 33);
    if (mom == lambda_skew) mID |= (static_cast<uint64_t>(1) << 48);
    mID |= (static_cast<uint64_t>(priority) << 34);
    if (mom == variance_A) mID |= (static_cast<uint64_t>(1) << 45);

    ostringstream titlestream;   
    titlestream << mID;   
    string title = titlestream.str();   
        
    //    
    //  Binning      
    //    
    //  Interpolate() only works up to the bin center    
    //  of the first and last bin. To guarantee that     
    //  the full range is accessible we shift the min    
    //  and max of each axis.    
    //    
    tmin = fabs(tmin);    
    tmax = fabs(tmax);    
    if (tmin > tmax) swap(tmin, tmax);     
    if (logQ2) Q2min = log(Q2min);    
    if (logQ2) Q2max = log(Q2max);    
    if (logW2) W2min = log(W2min);    
    if (logW2) W2max = log(W2max);    
    if (logt)  tmin =  log(tmin);    
    if (logt)  tmax =  log(tmax);    
    double binWidthQ2 = (Q2max - Q2min)/(nbinsQ2-1);    
    double binWidthW2 = (W2max - W2min)/(nbinsW2-1);    
    double binWidthT =  (tmax - tmin)/(nbinsT-1);        
        
    //    
    //  Book 3D histo to hold table    
    //    
    m3DHisto = new TH3F("table", title.c_str(),
                        nbinsQ2, Q2min-binWidthQ2/2., Q2max+binWidthQ2/2.,    
                        nbinsW2, W2min-binWidthW2/2., W2max+binWidthW2/2.,    
                        nbinsT,  tmin-binWidthT/2.,   tmax+binWidthT/2.);    
    m3DHisto->SetDirectory(0);    
        
    return nbinsQ2*nbinsW2*nbinsT; // return total number of bins    
}

//
//   Overloaded function for UPC only
//
unsigned int Table::create(int nbinsX, double Xmin, double Xmax,
                           int nbinsT, double tmin, double tmax,
                           bool logX, bool logt, bool logC,
                           AmplitudeMoment mom, GammaPolarization pol,
                           unsigned int A, int vm,
                           DipoleModelType model, DipoleModelParameterSet pset,
                           const char* filename, unsigned char priority)
{
    if (m2DHisto != 0) {
        cout << "Table:create(): cannot create, table already exists." << endl;
        return 0;
    }
    
    if (nbinsX < 2) {
        cout << "Table:create(): number of bins in x must be at least 2." << endl;
        return 0;
    }
    if (nbinsT < 2) {
        cout << "Table:create(): number of bins in t must be at least 2." << endl;
        return 0;
    }

    if (filename) mFilename = string(filename); // name of file where table is written to
    
    //
    //  Create table ID first
    //
    mID = (vm << 16);
    mID |= (A << 8);
    mID |= (model << 5);
    mID |= (static_cast<uint64_t>(pset) << 42);
    if (mom == mean_A2) mID |= 1;
    if (pol == transverse) mID |= (1 << 1);
    if (logt)        mID |= (1 << 2);
    if (logC)        mID |= (static_cast<uint64_t>(1) << 32);
    if (mom == lambda_real) mID |= (static_cast<uint64_t>(1) << 33);
    if (mom == lambda_skew) mID |= (static_cast<uint64_t>(1) << 48);
    mID |= (static_cast<uint64_t>(priority) << 34);
    if (mom == variance_A) mID |= (static_cast<uint64_t>(1) << 45);
    mID |= (static_cast<uint64_t>(1) << 46);  // flag UPC
    if (logX) mID |= (static_cast<uint64_t>(1) << 47);

    ostringstream titlestream;
    titlestream << mID;
    string title = titlestream.str();
    
    //
    //  Binning
    //
    //  Interpolate() only works up to the bin center
    //  of the first and last bin. To guarantee that
    //  the full range is accessible we shift the min
    //  and max of each axis.
    //
    tmin = fabs(tmin);
    tmax = fabs(tmax);
    if (tmin > tmax) swap(tmin, tmax);
    if (logX) Xmin = log(Xmin);
    if (logX) Xmax = log(Xmax);
    if (logt)  tmin =  log(tmin);
    if (logt)  tmax =  log(tmax);
    double binWidthX = (Xmax - Xmin)/(nbinsX-1);
    double binWidthT = (tmax - tmin)/(nbinsT-1);
    
    //
    //  Book 2D histo to hold table
    //
    m2DHisto = new TH2F("table", title.c_str(),
                        nbinsX, Xmin-binWidthX/2., Xmax+binWidthX/2.,
                        nbinsT, tmin-binWidthT/2., tmax+binWidthT/2.);
    m2DHisto->SetDirectory(0);
    
    return nbinsX*nbinsT; // return total number of bins
}
    
int  Table::numberOfEntries() const    
{    
    if (m3DHisto) {    
        int nx = m3DHisto->GetXaxis()->GetNbins();    
        int ny = m3DHisto->GetYaxis()->GetNbins();    
        int nz = m3DHisto->GetZaxis()->GetNbins();    
        return nx*ny*nz;    
    }    
    else if (m2DHisto) {   // UPC case
        int nx = m2DHisto->GetXaxis()->GetNbins();
        int ny = m2DHisto->GetYaxis()->GetNbins();
        return nx*ny;
    }
    else
        return 0;    
}    
    
void Table::binXYZ(int globalBin, int& binx, int& biny, int& binz) const    
{    
    //    
    //  Find correct bins for each axis.    
    //  Don't use ROOT global bins here,    
    //  they are a mess since they cross    
    //  over underflow and overflow bins.    
    //    
    //  All bins returned count starting     
    //  at 1. The global bin starts at     
    //  0 as usual in C++.     
    //    
    int nx = m3DHisto->GetXaxis()->GetNbins();    
    int ny = m3DHisto->GetYaxis()->GetNbins();    
        
    binx = globalBin%nx;    
    biny = ((globalBin-binx)/nx)%ny;    
    binz = ((globalBin-binx)/nx -biny)/ny;    
    binx++; biny++; binz++;    
}    

//
//  UPC version of binXYZ()
//
void Table::binXY(int globalBin, int& binx, int& biny) const
{
    //
    //  Find correct bins for each axis.
    //  Don't use ROOT global bins here,
    //  they are a mess since they cross
    //  over underflow and overflow bins.
    //
    //  All bins returned count starting
    //  at 1. The global bin starts at
    //  0 as usual in C++.
    //
    int nx = m2DHisto->GetXaxis()->GetNbins();
    
    binx = globalBin%nx;
    biny = (globalBin-binx)/nx;
    binx++; biny++;
}

void Table::binCenter(int i, double& Q2, double& W2, double& t) const    
{    
    int binx, biny, binz;    
    binXYZ(i, binx, biny, binz);    
    // cout << i << '\t' << binx << '\t' << biny << '\t' << binz << endl;    
        
    double x = m3DHisto->GetXaxis()->GetBinCenter(binx);    
    double y = m3DHisto->GetYaxis()->GetBinCenter(biny);    
    double z = m3DHisto->GetZaxis()->GetBinCenter(binz);    
    Q2 = isLogQ2() ?  exp(x) : x;    
    W2 = isLogW2() ?  exp(y) : y;    
    t  = isLogT()  ? -exp(z) : -z;    
    if (t > 0) t = 0; // avoid numeric glitch   
}

//
//   UPC overload
//
void Table::binCenter(int i, double& xpom, double& t) const
{
    int binx, biny;
    binXY(i, binx, biny);
    // cout << i << '\t' << binx << '\t' << biny << endl;
    
    double x = m2DHisto->GetXaxis()->GetBinCenter(binx);
    double y = m2DHisto->GetYaxis()->GetBinCenter(biny);
    xpom = isLogX() ?  exp(x) : x;
    t    = isLogT()  ? -exp(y) : -y;
    if (t > 0) t = 0; // avoid numeric glitch
}

    
void Table::fill(int i, double val, double err)    
{    
    int binx, biny, binz;    
    if (isUPC())
        binXY(i, binx, biny);
    else
        binXYZ(i, binx, biny, binz);

    if (isLogContent()) {   
        if (val == 0) {
            // cout << "Table::fill(): warning, log scale requested but value = 0." << endl;
            val = numeric_limits<float>::min();   
        }   
        val = log(fabs(val));   
    }   
    if (isUPC())
        m2DHisto->SetBinContent(binx, biny, val);
    else
        m3DHisto->SetBinContent(binx, biny, binz, val);
    if(err!=0.) {
        if (isUPC())
            m2DHisto->SetBinError(binx, biny, err);
        else
            m3DHisto->SetBinError(binx, biny, binz, err);
    }
    mFillCounter++;  
      
    //  
    //  Check if backup is due  
    //  
    if (mBackupFrequence)   
        if (mFillCounter%mBackupFrequence == 0) backup(i);  
}    
    
bool Table::fexists(const char* filename) const    
{    
    ifstream ifs(filename);      
    return !ifs.fail();
}    
    
void Table::write(const char* filename)    
{    
    //  
    //  'filename' is optional argument. Null value implies  
    //  that one was already passed via create(). Check  
    //  that this is the case. If one is given, it is used  
    //  no matter if it already was defined or not.  
    //  
    if (filename)  
        mFilename = string(filename);   
    else {  
        if (mFilename.empty()) {  // should be defined but is not  
            cout << "Table::write(): Warning, no filename supplied. Will use 'sartre_table.root'." << endl;    
            mFilename = string("sartre_table.root");  
        }  
    }  
          
    //    
    //  Check if file exist. If so alter     
    //  the filename. Creating tables is     
    //  a time consuming business so we do    
    //  not want to lose precious data if    
    //  something goes wrong here and we do    
    //  want to prevent accidents.    
    //    
    while (fexists(mFilename.c_str())) {    
        mFilename += string(".save");    
    }    
               
    TFile hfile(mFilename.c_str(),"RECREATE");    
    if (isUPC())
        m2DHisto->Write();
    else
        m3DHisto->Write();
    hfile.Close();
           
    cout << "Table::write(): table stored in file '" << mFilename.c_str()<< "'." << endl;    
      
    if (mBackupFrequence && !mLastBackupFilename.empty())  
        remove(mLastBackupFilename.c_str());  
}    
    
bool Table::read(const char* filename)    
{    
    //    
    //  Read histogram into memory    
    //    
    TFile *file = TFile::Open(filename,"READ");    
    if (file == 0) {    
        cout << "Table::read(): failed opening file '" << filename << "'." << endl;     
        return false;    
    }    
    
    //
    //  Clean up
    //
    delete m3DHisto; m3DHisto = 0;
    delete m2DHisto; m2DHisto = 0;
    
    //
    //  Read plain object first and check title which
    //  tells us if it is an UPC table or a std one.
    //  Then cast into the proper object.
    //
    auto ptr = file->Get("table");
    if (ptr == 0) {
        cout << "Table::read(): failed retrieving table from file '" << filename << "'." << endl;
        return false;
    }
    mID = atoll(ptr->GetTitle());  // for isUPC() to work
    
    if (isUPC()) {
        m2DHisto = reinterpret_cast<TH2F*>(ptr);
        m2DHisto->SetDirectory(0);
   }
    else {
        m3DHisto = reinterpret_cast<TH3F*>(ptr);
        m3DHisto->SetDirectory(0);
    }
    if (m2DHisto == 0 && m3DHisto == 0) {
        cout << "Table::read(): failed retrieving table from file '" << filename << "'." << endl;
        return false;
    }

    file->Close();
    
    mFilename = string(filename);
    return true;    
}    
    
int  Table::vectorMesonId() const {return ((mID >> 16) & 0xFFFF);}    
    
unsigned int Table::dipoleModelType() const {return ((mID >> 5) & 0x7);}

unsigned int Table::dipoleModelParameterSet() const {return ((mID >> 42) & 0x7);} 
    
unsigned int Table::A() const {return ((mID >> 8) & 0xFF);}    
    
bool Table::isLongitudinal() const {return !isTransverse();}    
    
bool Table::isTransverse() const {return (mID & (1 << 1));}    
    
GammaPolarization Table::polarization() const {return (isTransverse() ? transverse : longitudinal);}    
    
AmplitudeMoment Table::moment() const {
    if (isLambdaA()) 
        return lambda_real;
    else if (isMeanA())
        return mean_A;
    else if (isVarianceA())
        return variance_A;
    else if (isLambdaSkew())
        return lambda_skew;
    else
        return mean_A2;    
}

bool Table::isMeanA() const {
    // bits 0, 33, 45 and 48 are not set
    return !(mID & 1) && !(mID & (static_cast<uint64_t>(1) << 33)) && !(mID & (static_cast<uint64_t>(1) << 45))
      && !(mID & (static_cast<uint64_t>(1) << 48));
}
    
bool Table::isMeanA2() const {
    // bit 0 is set, bits 33, 45 and 48 are not set
    return (mID & 1) && !(mID & (static_cast<uint64_t>(1) << 33)) && !(mID & (static_cast<uint64_t>(1) << 45))
      && !(mID & (static_cast<uint64_t>(1) << 48));
}
    
bool Table::isLambdaA() const {
    // bit 33 is set, bits 0 and 45 are not set
    return !(mID & 1) && (mID & (static_cast<uint64_t>(1) << 33)) && !(mID & (static_cast<uint64_t>(1) << 45))
      && !(mID & (static_cast<uint64_t>(1) << 48));
}

bool Table::isLambdaSkew() const {
    // bit 48 is set, bits 0, 33, and 45 are not set
    return !(mID & 1) && (mID & (static_cast<uint64_t>(1) << 48)) && !(mID & (static_cast<uint64_t>(1) << 33))
      && !(mID & (static_cast<uint64_t>(1) << 45));
}

bool Table::isVarianceA() const {
    // bit 45 is set, bitss 0, 33, and 48 are not set
    return !(mID & 1) && !(mID & (static_cast<uint64_t>(1) << 33)) && (mID & (static_cast<uint64_t>(1) << 45))
      && !(mID & (static_cast<uint64_t>(1) << 48));
}

bool Table::isLogQ2() const {return (mID & (1 << 4));}
    
bool Table::isLogW2() const {return (mID & (1 << 3));}    
    
bool Table::isLogT() const {return (mID & (1 << 2));}    

bool Table::isLogX() const {return (mID & (static_cast<uint64_t>(1) << 47));}

bool Table::isLogContent() const {return (mID & (static_cast<uint64_t>(1) << 32));}    

unsigned int Table::priority() const {return ((mID >> 34) & 0xFF);}

bool Table::isUPC() const {return (mID & (static_cast<uint64_t>(1) << 46));}

uint64_t Table::id() const {return mID;}
   
double Table::binWidthQ2() const
{
    if (!m3DHisto) return 0;
    return m3DHisto->GetXaxis()->GetBinWidth(1);
}   

double Table::binWidthW2() const
{
    if (!m3DHisto) return 0;
    return m3DHisto->GetYaxis()->GetBinWidth(1);
}

double Table::binWidthT() const
{
    if (isUPC()) {
        return m2DHisto->GetYaxis()->GetBinWidth(1);
    }
    else
        return m3DHisto->GetZaxis()->GetBinWidth(1);
}

double Table::binWidthX() const
{
    if (!m2DHisto) return 0;
    return m2DHisto->GetXaxis()->GetBinWidth(1);
}

double Table::get(double Q2, double W2, double t) const
{    
    if (m3DHisto == 0) return 0;    
        
    //    
    // Transform variables to how they are stored in the table    
    //    
    double x = isLogQ2() ? log(Q2) : Q2;    
    double y = isLogW2() ? log(W2) : W2;    
    t = fabs(t);     
    double z = isLogT() ? log(t) : t;    
            
    //    
    //  Tiny numerical glitches will cause TH3F::Interpolate() to fail    
    //  since it requires that the variables lie between the first    
    //  and last bin center, excluding the centers.    
    //  Here we enforce that this is the case. The downside of this    
    //  is that *all* values of Q2, W2, t will be forced to lie within.    
    //  Hence the user has to make sure that the values are within    
    //  the boundaries (see minQ2(), maxQ2(), minW2() etc.) before     
    //  calling get().    
    //  In this case the corrections below are tiny and are only    
    //  applied to avoid minor glitches that do not affect the results.    
    //      
        
    TAxis *axis = m3DHisto->GetXaxis();    
    x = max(x, axis->GetBinCenter(1)+numeric_limits<float>::epsilon());    
    x = min(x, axis->GetBinCenter(axis->GetNbins())-numeric_limits<float>::epsilon());    
    axis = m3DHisto->GetYaxis();    
    y = max(y, axis->GetBinCenter(1)+numeric_limits<float>::epsilon());    
    y = min(y, axis->GetBinCenter(axis->GetNbins())-numeric_limits<float>::epsilon());    
    axis = m3DHisto->GetZaxis();    
    z = max(z, axis->GetBinCenter(1)+numeric_limits<float>::epsilon());    
    z = min(z, axis->GetBinCenter(axis->GetNbins())-numeric_limits<float>::epsilon());    
        
    // double result = InterpolateGridSpline(x, y, z);   // tmp uncommented until 0's in tables are cleared  
    // double result = m3DHisto->Interpolate(x, y, z);
    double result = modInterpolation(x, y, z);

    if (result == 0 && isLogContent()) {
      cout << "Table::get(): warning, 0 is not a valid table content when working in log scale." << endl;
    }
    if (isLogContent()) result = exp(result);
    return result;    
}    

//
//  UPC version
//
double Table::get(double xpom, double t) const
{
    if (m2DHisto == 0) return 0;

    //
    // Transform variables to how they are stored in the table
    //
    double x = isLogX() ? log(xpom) : xpom;
    t = fabs(t);
    double y = isLogT() ? log(t) : t;

    //
    //  Tiny numerical glitches will cause TH2F::Interpolate() to fail
    //  since it requires that the variables lie between the first
    //  and last bin center, excluding the centers.
    //  Here we enforce that this is the case. The downside of this
    //  is that *all* values of x, t will be forced to lie within.
    //  Hence the user has to make sure that the values are within
    //  the boundaries (see minX(), maxX(), minT() etc.) before
    //  calling get().
    //  In this case the corrections below are tiny and are only
    //  applied to avoid minor glitches that do not affect the results.
    //
    
    TAxis *axis = m2DHisto->GetXaxis();
    x = max(x, axis->GetBinCenter(1)+numeric_limits<float>::epsilon());
    x = min(x, axis->GetBinCenter(axis->GetNbins())-numeric_limits<float>::epsilon());
    axis = m2DHisto->GetYaxis();
    y = max(y, axis->GetBinCenter(1)+numeric_limits<float>::epsilon());
    y = min(y, axis->GetBinCenter(axis->GetNbins())-numeric_limits<float>::epsilon());

    double result = modInterpolation(x, y);  // xpom, t

    if (result == 0 && isLogContent()) {
        cout << "Table::get(): warning, 0 is not a valid table content when working in log scale." << endl;
    }
    if (isLogContent()) result = exp(result);
    return result;
}

double Table::minQ2() const    
{    
    if (m3DHisto == 0) return 0;    
    TAxis *axis = m3DHisto->GetXaxis();    
    double val = axis->GetBinCenter(1);     
    return isLogQ2() ? exp(val) : val;     
}    
    
double Table::maxQ2() const    
{    
    if (m3DHisto == 0) return 0;    
    TAxis *axis = m3DHisto->GetXaxis();    
    double val = axis->GetBinCenter(axis->GetNbins());     
    return isLogQ2() ? exp(val) : val;     
}    
    
double Table::minW2() const    
{    
    if (m3DHisto == 0) return 0;    
    TAxis *axis = m3DHisto->GetYaxis();    
    double val = axis->GetBinCenter(1);     
    return isLogW2() ? exp(val) : val;     
}    
    
double Table::maxW2() const    
{    
    if (m3DHisto == 0) return 0;    
    TAxis *axis = m3DHisto->GetYaxis();    
    double val = axis->GetBinCenter(axis->GetNbins());     
    return isLogW2() ? exp(val) : val;     
}    
    
double Table::minT() const    
{    
    if (m2DHisto == 0 && m3DHisto == 0) return 0;
    TAxis *axis = 0;
    if (isUPC())
        axis = m2DHisto->GetYaxis();
    else
        axis = m3DHisto->GetZaxis();
    double val = axis->GetBinCenter(axis->GetNbins()); // t always as |t| in table
    return isLogT() ? -exp(val) : -val;     
}    
    
double Table::maxT() const    
{    
    if (m2DHisto == 0 && m3DHisto == 0) return 0;
    TAxis *axis = 0;
    if (isUPC())
        axis = m2DHisto->GetYaxis();
    else
        axis = m3DHisto->GetZaxis();
    double val = axis->GetBinCenter(1); // t always as |t| in table    
    double result = isLogT() ? -exp(val) : -val;   
    if (result > 0) {   
        // cout << "Table::maxT(): warning, t > 0, t (" << result << ") set to 0 now." << endl;   
        result = 0;   
    }   
    return result;    
}    

double Table::minX() const
{
    if (m2DHisto == 0) return 0;
    TAxis *axis = m2DHisto->GetXaxis();
    double val = axis->GetBinCenter(1);
    return isLogX() ? exp(val) : val;
}


double Table::maxX() const
{
    if (m2DHisto == 0) return 0;
    TAxis *axis = m2DHisto->GetXaxis();
    double val = axis->GetBinCenter(axis->GetNbins());
    return isLogX() ? exp(val) : val;
}

string Table::filename() const {return mFilename;}    

bool Table::writeToBinaryFile(const string& filename, bool verbose) const
{
    //
    //  Open binary file
    //
    ofstream binaryFile;
    binaryFile.open (filename, ios::out | ios::binary);
    if (binaryFile.fail()) {
        cout << "Table::writeToBinaryFile: failed to open file '" << filename << "'." << endl;
        return false;
    }
    
    //
    //  Loop over all entries
    //
    double Q2, W2, xpom, t;
    int binx, biny, binz;
    double rawContent;
    for (int i=0; i<numberOfEntries(); i++) {
        if (isUPC()) {
            binXY(i, binx, biny);
            binCenter(i, xpom, t);
            rawContent = m2DHisto->GetBinContent(binx, biny);
            if (isLogContent()) {
                rawContent = exp(rawContent);
            }
            binaryFile.write(reinterpret_cast<const char *>(&i), sizeof(i));
            binaryFile.write(reinterpret_cast<const char *>(&xpom), sizeof(xpom));
            binaryFile.write(reinterpret_cast<const char *>(&t), sizeof(t));
            binaryFile.write(reinterpret_cast<const char *>(&rawContent), sizeof(rawContent));
            
            if (verbose)
                cout << "i=" << i << "\txpom=" << xpom << "\tt=" << t  << "\tvalue=" << rawContent << endl;
        }
        else {
            binXYZ(i, binx, biny, binz);
            binCenter(i, Q2, W2, t);
            rawContent = m3DHisto->GetBinContent(binx, biny, binz);
            if (isLogContent()) {
                rawContent = exp(rawContent);
            }
            binaryFile.write(reinterpret_cast<const char *>(&i), sizeof(i));
            binaryFile.write(reinterpret_cast<const char *>(&Q2), sizeof(Q2));
            binaryFile.write(reinterpret_cast<const char *>(&W2), sizeof(W2));
            binaryFile.write(reinterpret_cast<const char *>(&t), sizeof(t));
            binaryFile.write(reinterpret_cast<const char *>(&rawContent), sizeof(rawContent));
            
            if (verbose)
                cout << "i=" << i << "\tQ2=" << Q2 << "\tW2=" << W2 << "\tt=" << t  << "\tvalue=" << rawContent << endl;
        }
    }
    
    binaryFile.close();
    return true;
}

void Table::list(ostream& os, bool printContent, bool printStatistics, int startBin, int endBin) const
{
    //
    //  List table info as compact as possible
    //
    
    ios::fmtflags fmt = os.flags();  // store current I/O flags     
        
    if ((isUPC() && !m2DHisto) || (!isUPC() && !m3DHisto)) {    
        os << "Table::list(): table is undefined." << endl;    
        return;    
    }    
    
    //
    //  First row: table id, A, and content
    //
    os << setw(11) << "Table ID = " << setw(16) << left << mID;
    os << setw(7) << right << "A = " << setw(4) << left << A();
    os << setw(20) << right << "content = " << (isLogContent() ? "log(" : "");
    if (isMeanA()) 
        os << "<A>";
    else if (isMeanA2())
        os << "<A^2>";
    else if (isLambdaA())
        os << "lambda_<A>";
    else if (isVarianceA())
        os << "variance_A";
    else if (isLambdaSkew())
        os << "lambda_skew";
    else
        os << "unknown";
    os << (isLogContent() ? ")" : "") << endl;
    
    //
    //  Second row: vm and polarization
    //
    os << setw(34) << right << "vmId = " << setw(4) << left << vectorMesonId();
    os << setw(20) << right << "polarization = " <<  (isTransverse() ? 'T' : 'L') << endl;
    
    //
    //  Third row: dipole model and parameter set
    //
    os << setw(34) << right << "model = " << setw(7) << left;
    if (dipoleModelType() == bSat)   
        os << "bSat   ";   
    else if (dipoleModelType() == bNonSat)   
        os << "bNonSat";   
    else    
        os << "bCGC   ";
    os << setw(17) << right << "parameter set = ";
    if (dipoleModelParameterSet() == KMW)
        os << "KMW ";
    else if (dipoleModelParameterSet() == HMPZ)
        os << "HMPZ";
    else if (dipoleModelParameterSet() == CUSTOM)
        os << "CUSTOM";
    else
        os << "?   ";
    os << endl;
    
    //
    //  Third row: priority and UPC
    //
    os << setw(34) << right << "priority = " << setw(4) << left << priority();
    os << setw(20) << right << "UPC = " << (isUPC() ? "yes" : "no") << endl;

    //
    //  Next three rows: Q2, W2, t bins, range and log/linear
    //  or for UPC: x, t
    //
    if (isUPC()) {
        os << setw(35) << right << "xp = [" << minX() << ", " << maxX()
        << "; bins = " << m2DHisto->GetXaxis()->GetNbins() << (isLogX() ? "; log]" : "; lin]") << endl;
        os << setw(35) << right << "t = [" << minT() << ", " << maxT()
        << "; bins = " << m2DHisto->GetYaxis()->GetNbins() << (isLogT() ? "; log]" : "; lin]") << endl;
    }
    else {
        os << setw(35) << right << "Q2 = [" << minQ2() << ", " << maxQ2()
        << "; bins = " << m3DHisto->GetXaxis()->GetNbins() << (isLogQ2() ? "; log]" : "; lin]") << endl;
        os << setw(35) << right << "W2 = [" << minW2() << ", " << maxW2()
        << "; bins = " << m3DHisto->GetYaxis()->GetNbins() << (isLogW2() ? "; log]" : "; lin]") << endl;
        os << setw(35) << right << "t = [" << minT() << ", " << maxT()
        << "; bins = " << m3DHisto->GetZaxis()->GetNbins() << (isLogT() ? "; log]" : "; lin]") << endl;
    }
    
    //
    //  Filename at the end
    //
    os << setw(34) << right << "file = " << mFilename.c_str() << endl;
    os << endl;    
   
    //
    //  Print content (optional)
    //
    double Q2, W2, t, xpom;
    int binx, biny, binz;   
   
    if (printContent) {
        int thePrecision = os.precision();
        double rawContent = 0;
        if(endBin==0)  endBin=numberOfEntries();
        for (int i=startBin; i<endBin; i++) {
            if (isUPC()) {
                binXY(i, binx, biny);
                binCenter(i, xpom, t);
                double value = get(xpom, t);
                rawContent = m2DHisto->GetBinContent(binx, biny);
                double error = m2DHisto->GetBinError(binx, biny);
                os << "bin = "     << setw(8)  << left << i;
                os << "xp = "      << setw(13) << left << fixed << setprecision(5) << xpom;
                os << "t = "       << setw(16) << left << scientific << t;
                os << "value = "   << setw(15) << left << scientific << value;
                os << "(binx = "   << setw(5)  << left << binx;
                os << "biny = "    << setw(5)  << left << biny;
                os << "content = " << setw(14) << left << scientific << rawContent;
                os << "error = " << left << error << ')';
            }
            else {
                binXYZ(i, binx, biny, binz);
                binCenter(i, Q2, W2, t);
                double value = get(Q2, W2, t);
                rawContent = m3DHisto->GetBinContent(binx, biny, binz);
                double error = m3DHisto->GetBinError(binx, biny, binz);
                os << "bin = "     << setw(8)  << left << i;
                os << "Q2 = "      << setw(13) << left << fixed << setprecision(5) << Q2;
                os << "W2 = "      << setw(12) << left << fixed << W2;
                os << "t = "       << setw(16) << left << scientific << t;
                os << "value = "   << setw(15) << left << scientific << value;
                os << "(binx = "   << setw(5)  << left << binx;
                os << "biny = "    << setw(5)  << left << biny;
                os << "binz = "    << setw(5)  << left << binz;
                os << "content = " << setw(14) << left << scientific << rawContent;
                os << "error = " << left << error << ')';
            }
            if ( (isLogContent() && rawContent <= log(numeric_limits<float>::min()*2)) || (rawContent == 0)) cout << " Z";
            cout << endl;
        }
        os << endl;
        os.precision(thePrecision);
    }

    //
    //  Print statistics (optional)
    //
    if (printStatistics) {
        int nEmpty = 0;   
        int nInvalid = 0;   
        int nNegative = 0;   
        double sum = 0;   
        double maxContent = 0;   
        double minContent = numeric_limits<float>::max();   
        int minBin, maxBin;   
        minBin = maxBin = 0;  
        double c = 0;
        for (int i=0; i<numberOfEntries(); i++) {
            if (isUPC()) {
                binXY(i, binx, biny);
                c = m2DHisto->GetBinContent(binx, biny);
            }
            else {
                binXYZ(i, binx, biny, binz);
                c = m3DHisto->GetBinContent(binx, biny, binz);
            }

            if (c == 0) nEmpty++;
            
            if ( !isLogContent() && c < 0) nNegative++;
            if (std::isnan(c) || std::isinf(c) || !std::isfinite(c)) nInvalid++;   
            if (isLogContent()) c = exp(c);   
            if (c > maxContent) {maxContent = c; maxBin = i;}   
            if (c >= 0 && c < minContent) {minContent = c; minBin = i;}   
            if (c > 0) sum += c;   
        }   
           
        os << setw(34) << right << "total number of cells = " << numberOfEntries() << endl;
        os << setw(34) << right << "cells with negative content = " << nNegative << endl;
        os << setw(34) << right << "cells with no (0) content = " << nEmpty << endl;
        os << setw(34) << right << "cells with invalid content = " << nInvalid << endl;
        if (isUPC()) {
            binCenter(minBin, xpom, t);
            os << setw(34) << right << "minimum content = " << minContent << " at xp = " << xpom << ", t = " << t << endl;
            binCenter(maxBin, xpom, t);
            os << setw(34) << right << "maximum content = " << maxContent << " at xp = " << xpom << ", t = " << t << endl;
        }
        else {
            binCenter(minBin, Q2, W2, t);
            os << setw(34) << right << "minimum content = " << minContent << " at Q2 = " << Q2 << ", W2 = " << W2 << ", t = " << t << endl;
            binCenter(maxBin, Q2, W2, t);
            os << setw(34) << right << "maximum content = " << maxContent << " at Q2 = " << Q2 << ", W2 = " << W2 << ", t = " << t << endl;
        }
        os << setw(34) << right << "sum = " << sum << endl;
        os << endl;   
    }   
           
    os.flags(fmt);  // restore I/O flags     
}       

double Table::InterpolateGridSpline(double x, double y, double z) const // Q2, W2, t
{       
    //   
    //  The algorithm was taken from Cristian Constantin Lalescu.   
    //  See http://arxiv.org/abs/0905.3564 for details.   
    //   
    //  Grid points:  4, 6, or 8   
    //  Spline order: 3,5 for order 4   
    //                3,5,7,9 for order 6   
    //                3,5,7,9,11,13 for order 8   
    //   
       
    //   
    //  Find cell that refer to the position   
    //   
    int nx = m3DHisto->GetXaxis()->FindBin(x);    
    int ny = m3DHisto->GetYaxis()->FindBin(y);    
    int nz = m3DHisto->GetZaxis()->FindBin(z);    
       
    //   
    //  Define the # of grid points depending on how far we are away   
    //  from the edge. If at the edge we fall back to linear interpolation   
    //  as provided by TH3.   
    //  There must be a smarter way doing this than the code below   
    //   
    int gridPoints = 0;   
  
    if (nx-3 >= 1 && nx+4 <= m3DHisto->GetXaxis()->GetNbins() &&   
        ny-3 >= 1 && ny+4 <= m3DHisto->GetYaxis()->GetNbins() &&   
        nz-3 >= 1 && nz+4 <= m3DHisto->GetZaxis()->GetNbins()) {   
        gridPoints = 8;   
    }   
    else if (nx-2 >= 1 && nx+3 <= m3DHisto->GetXaxis()->GetNbins() &&   
             ny-2 >= 1 && ny+3 <= m3DHisto->GetYaxis()->GetNbins() &&   
             nz-2 >= 1 && nz+3 <= m3DHisto->GetZaxis()->GetNbins()) {   
        gridPoints = 6;   
    }   
    else if (nx-1 >= 1 && nx+2 <= m3DHisto->GetXaxis()->GetNbins() &&   
             ny-1 >= 1 && ny+2 <= m3DHisto->GetYaxis()->GetNbins() &&   
             nz-1 >= 1 && nz+2 <= m3DHisto->GetZaxis()->GetNbins()) {   
        gridPoints = 4;   
    }   
    else {   
        return m3DHisto->Interpolate(x,y,z);   
    }   
      
    //   
    // Find bin centers   
    //   
    double xc = m3DHisto->GetXaxis()->GetBinCenter(nx);   
    double yc = m3DHisto->GetYaxis()->GetBinCenter(ny);   
    double zc = m3DHisto->GetZaxis()->GetBinCenter(nz);   
           
    //   
    // Define the scale for coordinate transformations    
    // grid_spline expects x,y,z in bin units   
    //   
    double xscale = m3DHisto->GetXaxis()->GetBinWidth(1);   
    double yscale = m3DHisto->GetYaxis()->GetBinWidth(1);   
    double zscale = m3DHisto->GetZaxis()->GetBinWidth(1);   
   
    //   
    // Prepare grid spline alogorithm   
    //   
    double result, splineOrder;   
    if (gridPoints == 8) {   
        local_scal_3D<8> lf;   
        for (int i=-3;i<=4;i++) {   
            for(int j=-3;j<=4;j++) {                                     
                lf(-3,i,j) = m3DHisto->GetBinContent(nx-3, ny+i, nz+j);  
                lf(-2,i,j) = m3DHisto->GetBinContent(nx-2, ny+i, nz+j);  
                lf(-1,i,j) = m3DHisto->GetBinContent(nx-1, ny+i, nz+j);  
                lf( 0,i,j) = m3DHisto->GetBinContent(nx  , ny+i, nz+j);  
                lf( 1,i,j) = m3DHisto->GetBinContent(nx+1, ny+i, nz+j);  
                lf( 2,i,j) = m3DHisto->GetBinContent(nx+2, ny+i, nz+j);  
                lf( 3,i,j) = m3DHisto->GetBinContent(nx+3, ny+i, nz+j);  
                lf( 4,i,j) = m3DHisto->GetBinContent(nx+4, ny+i, nz+j);  
            }   
        }   
        splineOrder = 7;   
        result = grid_spline(splineOrder,lf, (x-xc)/xscale,(y-yc)/yscale,(z-zc)/zscale);   
    }   
    else if (gridPoints == 6) {   
        local_scal_3D<6> lf;   
        for (int i=-2;i<=3;i++) {   
            for(int j=-2;j<=3;j++) {                                     
                lf(-2,i,j) = m3DHisto->GetBinContent(nx-2, ny+i, nz+j);  
                lf(-1,i,j) = m3DHisto->GetBinContent(nx-1, ny+i, nz+j);  
                lf( 0,i,j) = m3DHisto->GetBinContent(nx  , ny+i, nz+j);  
                lf( 1,i,j) = m3DHisto->GetBinContent(nx+1, ny+i, nz+j);  
                lf( 2,i,j) = m3DHisto->GetBinContent(nx+2, ny+i, nz+j);  
                lf( 3,i,j) = m3DHisto->GetBinContent(nx+3, ny+i, nz+j);  
            }   
        }           
        splineOrder = 5;   
        result = grid_spline(splineOrder,lf, (x-xc)/xscale,(y-yc)/yscale,(z-zc)/zscale);   
    }           
    else if (gridPoints == 4) {   
        local_scal_3D<4> lf;   
        for (int i=-1;i<=2;i++) {   
            for(int j=-1;j<=2;j++) {                                     
                lf(-1,i,j) = m3DHisto->GetBinContent(nx-1, ny+i, nz+j);  
                lf( 0,i,j) = m3DHisto->GetBinContent(nx  , ny+i, nz+j);  
                lf( 1,i,j) = m3DHisto->GetBinContent(nx+1, ny+i, nz+j);  
                lf( 2,i,j) = m3DHisto->GetBinContent(nx+2, ny+i, nz+j);  
            }   
        }           
        splineOrder = 3;    
        result = grid_spline(splineOrder,lf, (x-xc)/xscale,(y-yc)/yscale,(z-zc)/zscale);   
    }     
    else {  
        cout << "Table::InterpolateGridSpline(): Error, illegal number of grid points." << endl;  
        result = 0;  
    }  
          
    return result;   
}   

double Table::modInterpolation(double x, double y, double z) const{  // Q2, W2, t
    
    //
    //  After testing many points: The best results is had for content in log and rest in lin.
    //  However, if any of the bins that participate in the interpolation are zero
    //  the interpolation has to be linear.
    //
    
    //  Make sure the point is within the histogram:
    if(m3DHisto->GetXaxis()->GetBinCenter(1) > x ||
       m3DHisto->GetXaxis()->GetBinCenter(m3DHisto->GetXaxis()->GetNbins()) < x ||
       m3DHisto->GetYaxis()->GetBinCenter(1) > y ||
       m3DHisto->GetYaxis()->GetBinCenter(m3DHisto->GetYaxis()->GetNbins()) < y ||
       m3DHisto->GetZaxis()->GetBinCenter(1) > z ||
       m3DHisto->GetZaxis()->GetBinCenter(m3DHisto->GetZaxis()->GetNbins()) < z ){
        cout<<"Table::myInterpolation Error: point lies outside the limits of the table!"<<endl;
        return 0;
    }
    
    // Find the bins surrounding the point, lower bins and upper bins
    int lbx = m3DHisto->GetXaxis()->FindBin(x);
    if( x < m3DHisto->GetXaxis()->GetBinCenter(lbx) ) lbx--;
    int ubx=lbx+1;
    
    int lby = m3DHisto->GetYaxis()->FindBin(y);
    if( y < m3DHisto->GetYaxis()->GetBinCenter(lby) ) lby--;
    int uby=lby+1;
    
    int lbz = m3DHisto->GetZaxis()->FindBin(z);
    if( z < m3DHisto->GetZaxis()->GetBinCenter(lbz) ) lbz--;
    int ubz=lbz+1;
    
    //The corresponding bin-centers:
    double ux=m3DHisto->GetXaxis()->GetBinCenter(ubx);
    double lx=m3DHisto->GetXaxis()->GetBinCenter(lbx);
    double uy=m3DHisto->GetYaxis()->GetBinCenter(uby);
    double ly=m3DHisto->GetYaxis()->GetBinCenter(lby);
    double uz=m3DHisto->GetZaxis()->GetBinCenter(ubz);
    double lz=m3DHisto->GetZaxis()->GetBinCenter(lbz);
    
    ux = isLogQ2() ? exp(ux) : ux;
    lx = isLogQ2() ? exp(lx) : lx;
    uy = isLogW2() ? exp(uy) : uy;
    ly = isLogW2() ? exp(ly) : ly;
    uz = isLogT()  ? exp(uz) : uz;
    lz = isLogT()  ? exp(lz) : lz;
    
    x = isLogQ2() ? exp(x) : x;
    y = isLogW2() ? exp(y) : y;
    z = isLogT()  ? exp(z) : z;
    
    //Find the distance to the point:
    double xd = (x-lx)/(ux-lx);
    double yd = (y-ly)/(uy-ly);
    double zd = (z-lz)/(uz-lz);
    
    // Make a vector containing all the values of the 8 surrounding bins
    double v[8] = { m3DHisto->GetBinContent(lbx, lby, lbz),  m3DHisto->GetBinContent(lbx, lby, ubz),
        m3DHisto->GetBinContent(lbx, uby, lbz),  m3DHisto->GetBinContent(lbx, uby, ubz),
        m3DHisto->GetBinContent(ubx, lby, lbz),  m3DHisto->GetBinContent(ubx, lby, ubz),
        m3DHisto->GetBinContent(ubx, uby, lbz),  m3DHisto->GetBinContent(ubx, uby, ubz) };
    
    bool logC=true;
    for(int i=0; i<8; i++){
        if(isLogContent() and v[i] <= log(numeric_limits<float>::min()*2) )
            logC=false;
        if(!isLogContent() and v[i] == 0)
            logC=false;
    }
    
    // Make interpolation in log(content) except for when at least one of the points are zero.
    if(isLogContent() and !logC){
        for(int i=0; i<8; i++){
            if (v[i] <= log(numeric_limits<float>::min()*2))
                v[i] = 0;
            else
                v[i]=exp(v[i]);
        }
    }
    if(!isLogContent() and logC){
        for(int i=0; i<8; i++)
            v[i]=log(v[i]);
    }
    
    //  First make four 1D interpolations in the z-direction
    double i1 = v[0] * (1 - zd) + v[1] * zd;
    double i2 = v[2] * (1 - zd) + v[3] * zd;
    double j1 = v[4] * (1 - zd) + v[5] * zd;
    double j2 = v[6] * (1 - zd) + v[7] * zd;
    
    // Secondly, make two 1D interpolations in the y-direction using the values obtained.
    double w1 = i1 * (1 - yd) + i2 * yd;
    double w2 = j1 * (1 - yd) + j2 * yd;
    
    // Finaly, make 1D interpolation in the x-direction
    double result= w1 * (1-xd) + w2 * xd;
    
    //Reverse exp/log compensation to get real result back:
    if(isLogContent() and !logC){
        if(result == 0)
            result=log(numeric_limits<float>::min());
        else
            result=log(result);
    }
    if(!isLogContent() and logC)
        result=exp(result);
    
    return result;
}

//
//   UPC version
//
double Table::modInterpolation(double x, double y) const{  // x, t
    
    //  Make sure the point is within the histogram:
    if(m2DHisto->GetXaxis()->GetBinCenter(1) > x ||
       m2DHisto->GetXaxis()->GetBinCenter(m2DHisto->GetXaxis()->GetNbins()) < x ||
       m2DHisto->GetYaxis()->GetBinCenter(1) > y ||
       m2DHisto->GetYaxis()->GetBinCenter(m2DHisto->GetYaxis()->GetNbins()) < y) {
        cout<<"Table::myInterpolation Error: point lies outside the limits of the table!"<<endl;
        return 0;
    }
    
    //   Bins that contain the requested values
    int bin_x = m2DHisto->GetXaxis()->FindBin(x);
    int bin_y = m2DHisto->GetYaxis()->FindBin(y);

    //  Find quadrant of the bin we are in
    int quadrant = 0;
    double dx = m2DHisto->GetXaxis()->GetBinUpEdge(bin_x)-x;
    double dy = m2DHisto->GetYaxis()->GetBinUpEdge(bin_y)-y;
    if (dx<=m2DHisto->GetXaxis()->GetBinWidth(bin_x)/2 && dy<=m2DHisto->GetYaxis()->GetBinWidth(bin_y)/2)
        quadrant = 1; // upper right
    if (dx>m2DHisto->GetXaxis()->GetBinWidth(bin_x)/2 && dy<=m2DHisto->GetYaxis()->GetBinWidth(bin_y)/2)
        quadrant = 2; // upper left
    if (dx>m2DHisto->GetXaxis()->GetBinWidth(bin_x)/2 && dy>m2DHisto->GetYaxis()->GetBinWidth(bin_y)/2)
        quadrant = 3; // lower left
    if (dx<=m2DHisto->GetXaxis()->GetBinWidth(bin_x)/2 && dy>m2DHisto->GetYaxis()->GetBinWidth(bin_y)/2)
        quadrant = 4; // lower right

    double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
    switch(quadrant) {
        case 1:
            x1 = m2DHisto->GetXaxis()->GetBinCenter(bin_x);
            y1 = m2DHisto->GetYaxis()->GetBinCenter(bin_y);
            x2 = m2DHisto->GetXaxis()->GetBinCenter(bin_x+1);
            y2 = m2DHisto->GetYaxis()->GetBinCenter(bin_y+1);
            break;
        case 2:
            x1 = m2DHisto->GetXaxis()->GetBinCenter(bin_x-1);
            y1 = m2DHisto->GetYaxis()->GetBinCenter(bin_y);
            x2 = m2DHisto->GetXaxis()->GetBinCenter(bin_x);
            y2 = m2DHisto->GetYaxis()->GetBinCenter(bin_y+1);
            break;
        case 3:
            x1 = m2DHisto->GetXaxis()->GetBinCenter(bin_x-1);
            y1 = m2DHisto->GetYaxis()->GetBinCenter(bin_y-1);
            x2 = m2DHisto->GetXaxis()->GetBinCenter(bin_x);
            y2 = m2DHisto->GetYaxis()->GetBinCenter(bin_y);
            break;
        case 4:
            x1 = m2DHisto->GetXaxis()->GetBinCenter(bin_x);
            y1 = m2DHisto->GetYaxis()->GetBinCenter(bin_y-1);
            x2 = m2DHisto->GetXaxis()->GetBinCenter(bin_x+1);
            y2 = m2DHisto->GetYaxis()->GetBinCenter(bin_y);
            break;
    }

    // Now find the bins we interpolate with
    int bin_x1 = m2DHisto->GetXaxis()->FindBin(x1);
    if (bin_x1<1) bin_x1=1;
    int bin_x2 = m2DHisto->GetXaxis()->FindBin(x2);
    if (bin_x2>m2DHisto->GetXaxis()->GetNbins()) bin_x2=m2DHisto->GetXaxis()->GetNbins();
    int bin_y1 = m2DHisto->GetYaxis()->FindBin(y1);
    if (bin_y1<1) bin_y1=1;
    int bin_y2 = m2DHisto->GetYaxis()->FindBin(y2);
    if (bin_y2>m2DHisto->GetYaxis()->GetNbins()) bin_y2=m2DHisto->GetYaxis()->GetNbins();
 
    //  Get content
    int bin_q22 = m2DHisto->GetBin(bin_x2,bin_y2);
    int bin_q12 = m2DHisto->GetBin(bin_x1,bin_y2);
    int bin_q11 = m2DHisto->GetBin(bin_x1,bin_y1);
    int bin_q21 = m2DHisto->GetBin(bin_x2,bin_y1);
    double q11 = m2DHisto->GetBinContent(bin_q11);
    double q12 = m2DHisto->GetBinContent(bin_q12);
    double q21 = m2DHisto->GetBinContent(bin_q21);
    double q22 = m2DHisto->GetBinContent(bin_q22);
    
    //
    //   As explained in the 3D version we interpolate on linear x, y
    //   but on log content, except when 0's are involved.
    //
    bool logC = true;
    if ((isLogContent() && q11 <= log(numeric_limits<float>::min()*2)) ||
        (!isLogContent() && q11 <= 0)) logC = false;
    if ((isLogContent() && q12 <= log(numeric_limits<float>::min()*2)) ||
        (!isLogContent() && q12 <= 0)) logC = false;
    if ((isLogContent() && q21 <= log(numeric_limits<float>::min()*2)) ||
        (!isLogContent() && q21 <= 0)) logC = false;
    if ((isLogContent() && q22 <= log(numeric_limits<float>::min()*2)) ||
        (!isLogContent() && q22 <= 0)) logC = false;

    if(isLogContent() && !logC) {   // 0's present, use linear content
        if (q11 <= log(numeric_limits<float>::min()*2))
            q11 = 0;
        else
            q11=exp(q11);
        if (q12 <= log(numeric_limits<float>::min()*2))
            q12 = 0;
        else
            q12=exp(q12);
        if (q21 <= log(numeric_limits<float>::min()*2))
            q21 = 0;
        else
            q21=exp(q21);
        if (q22 <= log(numeric_limits<float>::min()*2))
            q22 = 0;
        else
            q22=exp(q22);
    }

    if (!isLogContent() && logC) {
        q11 = log(q11);
        q12 = log(q12);
        q21 = log(q21);
        q22 = log(q22);
    }

    x = isLogX() ? exp(x) : x;
    x1 = isLogX() ? exp(x1) : x1;
    x2 = isLogX() ? exp(x2) : x2;
    y = isLogT() ? exp(y) : y;
    y1 = isLogT() ? exp(y1) : y1;
    y2 = isLogT() ? exp(y2) : y2;

    //  Interpolation
    double d = 1.0*(x2-x1)*(y2-y1);
    double result = 1.0*q11/d*(x2-x)*(y2-y)+1.0*q21/d*(x-x1)*(y2-y)+1.0*q12/d*(x2-x)*(y-y1)+1.0*q22/d*(x-x1)*(y-y1);
  
    //  Reverse exp/log compensation to get real result back:
    if (isLogContent() && !logC){
        if (result == 0)
            result=log(numeric_limits<float>::min());
        else
            result=log(result);
    }
    if (!isLogContent() && logC)
        result=exp(result);

    return result;
}

void Table::setAutobackup(const char* prefix, int freq)   
{  
    mBackupPrefix = string(prefix);  
    mBackupFrequence = freq;  
}  
  
void Table::backup(int backupBin)  
{  
    ostringstream backupFilenameStream;  
    backupFilenameStream << mBackupPrefix.c_str() << "_backup." << static_cast<int>(getpid())   
                         << '.' << mID << '.' << backupBin <<  ".root";  
  
    string filename = backupFilenameStream.str();  
    TFile file(filename.c_str(),"RECREATE");   
    m3DHisto->Write();   
    file.Close();   
      
    if (filename != mLastBackupFilename)  
        remove(mLastBackupFilename.c_str());  
    mLastBackupFilename = filename;  
    time_t now = time(0);  
    cout << "Table::backup(): autobackup performed, file = '" << mLastBackupFilename.c_str() << "', time = " << ctime(&now);  
}  
  
int Table::globalBin(int binx, int biny, int binz) const  
{  
    int nbinx = m3DHisto->GetXaxis()->GetNbins();  
    int nbiny = m3DHisto->GetYaxis()->GetNbins();  
      
    return (binz-1)*nbinx*nbiny+(biny-1)*nbinx+(binx-1);  
}  

//
//   UPC version
//
int Table::globalBin(int binx, int biny) const
{
    int nbinx = m3DHisto->GetXaxis()->GetNbins();
    return (biny-1)*nbinx+(binx-1);
}
