//==============================================================================
//  tableInspector.cpp
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
//  Utility program to check the integrity of lookup tables used in Sartre.   
//  Usage:  tableInspector [-x bin -y bin -z bin | -g bin | [-a] [-s] [-b]] file(s) ...  
//          -b  print table id in binary representation first
//          -a  print position and values of each cell   
//          -s  print statstics of table content (min, max, empty bines etc.)   
//          -g bin  returns the ROOT/table bins in x,y,z for a given global bin  
//          -x bin -y bin -z bin   
//                  returns the global bin for a given set of ROOT/table bins  
//==============================================================================   
#include "Table.h"   
#include <iostream>   
#include <unistd.h>   
#include <cstdlib>  
#include <string>  
#include <unistd.h>
   
using namespace std;   
   
template<class T> string binary(T word)
{
    int nbits = 8*sizeof(T);
    string bitpattern(nbits+1, ' ');
    unsigned long long mask = 1;
    mask <<= (nbits-1);
    for (int i=0; i<nbits; i++) {
        bitpattern[i] = (((word&mask) == mask) ? '1' : '0');
        mask >>= 1;
    }
    bitpattern[nbits] = '\0';
    return bitpattern;
}

template<class T> void printBinary(T word)
{
    string bits = binary(word);
    cout << bits.c_str() << endl;
}

void usage(const char* prog)   
{   
    cout << "Usage: " << prog << " [ -x bin -y bin -z bin | -g bin | [-a] [-s] [-b]] file(s) ..." << endl;   
    cout << "                    " << "-b  print table id in binary representation first" << endl;
    cout << "                    " << "-a  print position and values of each cell" << endl;
    cout << "                    " << "-s  print statstics of table content (min, max, empty bines etc.)" << endl;
    cout << "                    " << "-g bin  returns the ROOT/table bins in x,y,z for a given global bin" << endl;
    cout << "                    " << "-x bin -y bin -z bin returns the global bin for a given set of ROOT/table bins" << endl;
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
       
    bool listall = false;   
    bool liststat = false;   
    bool calcXYZ = false;   
    bool calcGlobal = false;   
    bool binaryId = false;   
    int globalBin = 0;  
    int binx = 1;  
    int biny = 1;  
    int binz = 1;  
    int ch;   
    while ((ch = getopt(argc, argv, "absg:x:y:z:")) != -1) {   
        switch (ch) {   
            case 'a':   
                listall = true;   
                break;   
            case 'b':   
                binaryId = true;   
                break;   
            case 's':   
                liststat = true;   
                break;   
            case 'g':  
                globalBin = atol(optarg);   
                calcXYZ = true;  
                break;   
            case 'x':  
                binx = atol(optarg);   
                calcGlobal = true;  
                break;   
            case 'y':  
                biny = atol(optarg);   
                calcGlobal = true;  
                break;   
            case 'z':  
                binz = atol(optarg);   
                calcGlobal = true;  
                break;   
            case '?':   
            default:   
                usage(argv[0]);   
                return 2;   
                break;   
        }   
    }   
    if (optind == argc || ( calcXYZ && calcGlobal ) || ( (calcXYZ || calcGlobal) && (listall || liststat) ) ) {   
        usage(argv[0]);   
        return 2;   
    }   
       
    //   
    //  List table info (optionally list all entries)   
    //   
    for (int index = optind; index < argc; index++) {   
        Table tbl;   
        if (tbl.read(argv[index])) { 
            if (binaryId) {
                uint64_t theId = tbl.id();
                cout << "Table ID = ";
                printBinary(theId);
            }
            if (calcXYZ) {
                cout << argv[index] << ": ";
                if (tbl.isUPC()) {
                    tbl.binXY(globalBin, binx, biny);
                    cout << globalBin << " --> binx=" << binx << " (xp), biny=" << biny << " (t)" << endl;
                }
                else {
                    tbl.binXYZ(globalBin, binx, biny, binz);
                    cout << globalBin << " --> binx=" << binx << " (Q2), biny=" << biny << " (W2), binz=" << binz << " (t)" << endl;
                }
            }
            else if (calcGlobal) {  
                cout << argv[index] << ": ";
                if (tbl.isUPC()) {
                    globalBin = tbl.globalBin(binx, biny);
                    cout << "binx=" << binx << " (xp), biny=" << biny << " (t)  --> " << globalBin  << endl;
                }
                else {
                    globalBin = tbl.globalBin(binx, biny, binz);
                    cout << "binx=" << binx << " (Q2), biny=" << biny << " (W2), binz=" << binz << " (t)  --> " << globalBin  << endl;
                }
            }
            else {  
                tbl.list(cout, listall, liststat);   
            }  
        }          
    }   
    return 0;   
}   
