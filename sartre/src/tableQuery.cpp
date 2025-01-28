//==============================================================================
//  tableQuery.cpp
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
//  Allows to query the content of table(s) interactively.  
//  User is prompted for t, Q2, W values, program returns interpolated values.  
//  Useful to check interpolation mechanism.  
//  
//  Usage:  
//  tableQuery file(s)  
//    
//==============================================================================   
#include "Table.h"   
#include <iostream>   
#include <vector>   
#include <unistd.h>   
#include <sstream>  
#include <string>  
#include <ctype.h>  
#define PR(x) cout << #x << " = " << (x) << endl;    

using namespace std;   

template<class T>  
inline void Prompt(const char *text, T& var)  
{  
    string line;
    char   c;
    
    cout << text << " [" << var << "]: ";
    while ((c = cin.get()) && c != '\n') line += c;
    if (line.length() > 0) {
        istringstream ist(line);
        ist >> var;
    }
}  

inline void Prompt(const char *text, bool& var)  
{  
    string line;
    char   c;
    string svar = var ? "true" : "false";
    
    cout << text << " [" << svar.c_str() << "]: ";
    while ((c = cin.get()) && c != '\n') line += c;
    if (line.length() > 0) {
        if (line == "true")
            var = true;
        else if (line == "t")
            var = true;
        else if (line == "yes")
            var = true;
        else if (line == "y")
            var = true;
        else if (line == "on")
            var = true;
        else if (line == "1")
            var = true;
        else
            var = false;
    }
}  

void usage(const char* prog)   
{   
    cout << "Usage: " << prog << " file(s) ..." << endl;
}   

int main(int argc, char **argv)   
{   
    //
    //  Handle command line arguments
    //
    if (argc == 1) {
        usage(argv[0]);
        return 2;
    }
    
    //
    //  Store tables in vector
    //
    vector<Table*> tables;
    vector<Table*> upcTables;
    for (int i=1; i<argc; i++) {
        Table *oneTable = new Table;
        oneTable->read(argv[i]);
        if (oneTable->isUPC())
            upcTables.push_back(oneTable);
        else
            tables.push_back(oneTable);
    }
    
    if (tables.size() == 0 && upcTables.size() == 0) {
        cout << "tableQuery: no tables loaded." << endl;
    }
    
    if (tables.size() && upcTables.size()) {
        cout << "Some of the tables are UPC tables. Will start" << endl;
        cout << "to loop over non-UPC tables then UPC tables." << endl;
    }
    //
    //  Loop until user stops it
    //
    bool loop = true;
    double Q2 = 10;
    double W = 50;
    double t = -0.1;
    double xpom = 0.01;
    double W2;
    string bound;
    
    if (tables.size()) {
        while (loop) {
            Prompt("t", t);
            Prompt("Q2", Q2);
            Prompt("W", W);
            W2 = W*W;
            
            for (unsigned int i=0; i<tables.size(); i++) {
                double c = tables[i]->get(Q2, W2 , t);
                if (t >= tables[i]->minT() && t <= tables[i]->maxT() &&
                    Q2 >= tables[i]->minQ2() && Q2 <= tables[i]->maxQ2()&&
                    W2 >= tables[i]->minW2() && W2 <= tables[i]->maxW2())
                    bound = "";
                else
                    bound = " (outside boundary)";
                
                cout << tables[i]->filename().c_str() << " --> " << c << bound.c_str() << endl;
            }
            Prompt("continue", loop);
        }
    }
    
    loop = true;
    
    if (upcTables.size()) {
        while (loop) {
            Prompt("t", t);
            Prompt("xp", xpom);
            for (unsigned int i=0; i<upcTables.size(); i++) {
                double c = upcTables[i]->get(xpom , t);
                if (t >= upcTables[i]->minT() && t <= upcTables[i]->maxT() &&
                    xpom >= upcTables[i]->minX() && xpom <= upcTables[i]->maxX())
                    bound = "";
                else
                    bound = " (outside boundary)";
                
                cout << upcTables[i]->filename().c_str() << " --> " << c << bound.c_str() << endl;
            }
            Prompt("continue", loop);
        }
    }
    
    return 0;
}   
