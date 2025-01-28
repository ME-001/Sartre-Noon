//==============================================================================
//  TableCollection.h
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
#ifndef TableCollection_h         
#define TableCollection_h         
#include "Enumerations.h"         
#include <vector>         
#include <iostream>         
         
using namespace std;         
         
class Table;         
         
class TableCollection {         
public:         
    TableCollection();         
    TableCollection(const TableCollection&);         
    TableCollection(int A, DipoleModelType typ, DipoleModelParameterSet set, int vmID);
    ~TableCollection();         
   
    TableCollection& operator=(const TableCollection&);         
  
    bool init(int A, DipoleModelType typ, DipoleModelParameterSet set, int vmID);
    
    bool tableExists(GammaPolarization pol, AmplitudeMoment mom) const;         
    bool tableExists(AmplitudeMoment mom) const;   // UPC version

    bool available(double Q2, double W2, double t, GammaPolarization p, AmplitudeMoment m) const;
    bool available(double xpom, double t, AmplitudeMoment m) const;  // UPC version

    double get(double Q2, double W2, double t, GammaPolarization p, AmplitudeMoment m) const;         
    double get(double Q2, double W2, double t, GammaPolarization p, AmplitudeMoment m, Table *&) const;         
    double get(double xpom, double t, AmplitudeMoment m) const;     // UPC version
    double get(double xpom, double t, AmplitudeMoment m, Table *&) const;     // UPC version

    void list(ostream& = cout, bool = false) const;         
         
    double minQ2() const;           
    double maxQ2() const;         
    double minW2() const;         
    double maxW2() const;         
    double minW() const;         
    double maxW() const;         
    double minT() const;         
    double maxT() const;         
    double minX() const;
    double maxX() const;

private:         
    double minimumValue(unsigned int) const;   
    double maximumValue(unsigned int) const;   
  
private:         
    vector<Table*> mTables;         
};         
#endif         
