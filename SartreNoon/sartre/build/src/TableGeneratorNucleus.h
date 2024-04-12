//==============================================================================
//  TableGeneratorNucleus.h
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
#ifndef TableGeneratorNucleus_h         
#define TableGeneratorNucleus_h         
#include "Nucleon.h"         
#include "Nucleus.h"         
#include <vector>         
         
using namespace std;         
         
class TH1D;         
         
class TableGeneratorNucleus : public Nucleus {         
public:         
    TableGeneratorNucleus();         
    TableGeneratorNucleus(unsigned int A);         
    TableGeneratorNucleus(const TableGeneratorNucleus&);         
    ~TableGeneratorNucleus();
         
    TableGeneratorNucleus& operator=(const TableGeneratorNucleus&);         
  
    bool  generate();         
    const TH1D* getRHisto() const;
         
public:             
    vector<Nucleon> configuration;         
             
private:         
    TH1D* mRadialDistributionHistogram;         
};         
#endif         
