//==============================================================================
//  ExclusiveFinalStateGenerator.h
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
//  $Date: 2019-04-01 14:39:33 -0400 (Mon, 01 Apr 2019) $
//  $Author: ullrich $
//==============================================================================
#ifndef ExclusiveFinalStateGenerator_h         
#define ExclusiveFinalStateGenerator_h         
         
#include "FinalStateGenerator.h"         
         
class ExclusiveFinalStateGenerator : public FinalStateGenerator {         
public:         
    ExclusiveFinalStateGenerator();         
    ~ExclusiveFinalStateGenerator();         
             
    bool generate(int id, double t, double y, double Q2,          
                  bool isIncoherent, int A, Event *event);         

    //UPC version:
    bool generate(int id, double t, double xpom,          
                  bool isIncoherent, int A, Event *event);
    double xpomMin(double massVM, double t, TLorentzVector hBeam, TLorentzVector eBeam);
};         
#endif         
