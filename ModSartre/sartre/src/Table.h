//==============================================================================
//  Table.h
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
//  To read/load an existing table:         
//         
//  Table tbl;         
//  tbl.read("filename");  // read table into memory         
//  double content = tbl.get(Q2, W2, t); // return content for Q2, W2, t         
//         
//  To create a table:         
//          
//  Table tbl;         
//  int n = tbl.create(nbinQ2, Q2min, Q2max,         
//                     nBinW2, W2min, W2max,         
//                     nBinT,  tmin,  tmax,         
//                     logQ2, logW2, logT, logC   // all bools         
//                     AmplitudeMoment m, GammaPolarization p,       
//                     massA, vmPDG, model);         
//         
//  for (int i=0; i<n; i++) { // filling table         
//      double Q2, W2, t, content;         
//      tbl.binCenter(i, Q2, W2, t);         
//      // calculate cross-section for Q2, W2, t...         
//      tbl.fill(i, content);         
//  }         
//  tbl.write("filename");         
//         
//  Use Table::list() to query the definition and content of a table.        
//  
//  
//  Backup mechanism:  
//  
//  void setAutobackup(const char* prefix, int freq);  
//  
//  Will store the current table after each freq'th to fill().  
//  Filename is  <prefix>_backup.<pid>.<id>.<bin>.root  
//      where  
//      <prefix> is a unique name passed to setAutobackup()  
//      <pid> is the process ID (see getpid())  
//      <bin> is the number of the last bin filled  
//      <id>  is the table ID  
//  
//  If a new backup file is created the old/previous one gets deleted.  
//  The last valid backup file is deleted after a successful call to write().  
//==============================================================================       
#ifndef Table_h         
#define Table_h         
#include <iostream>         
#include <string>         
#include "Enumerations.h"         
#include <stdint.h>         

using namespace std;         
         
class TH2F;
class TH3F;

class Table {         
public:         
    Table();         
    Table(const Table&);         
    ~Table();         
  
    Table& operator=(const Table&);  
      
    //         
    //  Reading and using a table         
    //         
    bool read(const char*);         
    double get(double Q2, double W2, double t) const;         
    double get(double xpom, double t) const;  // UPC version

    //         
    //  Creating, filling, and storing a table         
    //         
    unsigned int create(int, double, double,         
                        int, double, double,         
                        int, double, double,         
                        bool logQ2, bool logW2, bool logt,         
                        bool logContent,       
                        AmplitudeMoment, GammaPolarization,         
                        unsigned int A, int vm,         
                        DipoleModelType model,
                        DipoleModelParameterSet pset,
                        const char* filename = 0,
                        unsigned char priority = 0);
    
    unsigned int create(int, double, double,        // UPC version
                        int, double, double,
                        bool logx, bool logt,
                        bool logContent,
                        AmplitudeMoment, GammaPolarization,
                        unsigned int A, int vm,
                        DipoleModelType model,
                        DipoleModelParameterSet pset,
                        const char* filename = 0,
                        unsigned char priority = 0);
    
             
    void binCenter(int, double& Q2, double& W2, double& t) const;
    void binCenter(int, double& xpom, double& t) const;  // UPC version

    
    void fill(int, double, double err = 0);  
    void write(const char* filename = 0);         
             
    //  
    //  Backup mechanism  
    //  
    void setAutobackup(const char* prefix, int freq);  
      
    //         
    //  Query functions         
    //         
    unsigned int A() const;      // mass number         
    bool isTransverse() const;         
    bool isLongitudinal() const;         
    GammaPolarization polarization() const;         
             
    bool isMeanA2() const;         
    bool isMeanA() const;         
    bool isLambdaA() const;         
    bool isLambdaSkew() const;         
    bool isVarianceA() const;
    AmplitudeMoment moment() const;
          
    int  vectorMesonId() const;         
    unsigned int  dipoleModelType() const;         
    unsigned int  dipoleModelParameterSet() const;
    int  numberOfEntries() const;
             
    double minQ2() const;           
    double maxQ2() const;         
    double minW2() const;         
    double maxW2() const;         
    double minT() const;         
    double maxT() const;         
    double minX() const;
    double maxX() const;

    double binWidthQ2() const;
    double binWidthW2() const;
    double binWidthT() const;
    double binWidthX() const;

    unsigned int priority() const;
    bool isUPC() const;
       
    uint64_t id() const;         
       
    void list(ostream& = cout, bool = false, bool = false, int = 0, int = 0) const;
    string filename() const;   
      
    int globalBin(int binx, int biny, int binz) const;
    int globalBin(int binx, int biny) const;  // UPC version
    

    void binXYZ(int globalBin, int& binx, int& biny, int& binz) const;         
    void binXY(int globalBin, int& binx, int& biny) const;    // UPC version

    bool writeToBinaryFile(const string&, bool verbose=false) const;
    
private:          
    bool isLogQ2() const;         
    bool isLogW2() const;         
    bool isLogT() const;         
    bool isLogX() const;
    bool isLogContent() const;
    bool fexists(const char*) const;         
           
    double InterpolateGridSpline(double x, double y, double z) const;    
    double modInterpolation(double, double, double) const;    
    double modInterpolation(double, double) const;  // UPC version

    void backup(int);  
             
private:         
    uint64_t       mID;         
    TH3F           *m3DHisto;
    TH2F           *m2DHisto;  // UPC table
    string         mFilename;
    int            mFillCounter;  
    int            mBackupFrequence;  
    string         mBackupPrefix;  
    string         mLastBackupFilename;  
};         
         
#endif   
