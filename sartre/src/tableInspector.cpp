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
//  $Date: 2019-09-04 12:48:48 -0400 (Wed, 04 Sep 2019) $
//  $Author: ullrich $
//==============================================================================
//   
//  Utility program to check the integrity of lookup tables used in Sartre.   
//  Usage:  tableInspector [-x bin -y bin -z bin | -g bin | [-a] [-s] [-b] | -e txt] file(s) ...
//          -b  print table id in binary representation first
//          -a  print position and values of each cell   
//          -s  print statstics of table content (min, max, empty bines etc.)   
//          -g bin  returns the ROOT/table bins in x,y,z for a given global bin  
//          -x bin -y bin -z bin   
//                  returns the global bin for a given set of ROOT/table bins  
//          -e txt  returns specific table info elements only. txt can be one of:
//                  A, pol, cont, vmid, vmstr, model, pset, type.
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
    cout << "Usage: " << prog << " [ -x bin -y bin -z bin | -g bin | [-a] [-s] [-b] | -e txt] file(s) ..." << endl;
    cout << "                    " << "-b  print table id in binary representation first" << endl;
    cout << "                    " << "-a  print position and values of each cell" << endl;
    cout << "                    " << "-s  print statstics of table content (min, max, empty bines etc.)" << endl;
    cout << "                    " << "-g bin  returns the ROOT/table bins in x,y,z for a given global bin" << endl;
    cout << "                    " << "-x bin -y bin -z bin returns the global bin for a given set of ROOT/table bins" << endl;
    cout << "                    " << "-e txt  returns specific table info element only. txt can be one of:" << endl;
    cout << "                    " << "        A, pol, cont, vmid, vmstr, model, pset, type. " << endl;
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
    bool showElement = false;
    int globalBin = 0;
    int binx = 1;  
    int biny = 1;  
    int binz = 1;  
    int ch;
    string element;
    while ((ch = getopt(argc, argv, "absg:x:y:z:e:")) != -1) {
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
            case 'e':
                showElement = true;
                element = string(optarg);
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
    //  If user choose to extract only a certain element out of
    //  all available info.
    //
    if (showElement) {
        for (int index = optind; index < argc; index++) {
            Table tbl;
            if (tbl.read(argv[index])) {
                if (element == "A") {
                    cout << tbl.A() << endl;
                }
                else if (element == "pol") {
                    if (tbl.isTransverse()) cout << "T"  << endl;
                    else if (tbl.isLongitudinal()) cout << "L"  << endl;
                    else cout << "?"  << endl;
                }
                else if (element == "cont") {
                    if (tbl.isMeanA2()) cout << "mean_A2" << endl;
                    else if (tbl.isMeanA()) cout << "mean_A" << endl;
                    else if (tbl.isLambdaA()) cout << "lambda_real" << endl;
                    else if (tbl.isLambdaSkew()) cout << "lambda_skew" << endl;
                    else if (tbl.isVarianceA()) cout << "variance_A" << endl;
                    else cout << "?"  << endl;
                }
                else if (element == "vmid") {
                    cout << tbl.vectorMesonId() << endl;
                }
                else if (element == "vmstr") {
                    int id = tbl.vectorMesonId();
                    if (id == 22) cout << "dvcs" << endl;
                    else if (id == 113) cout << "rho" << endl;
                    else if (id == 333) cout << "phi" << endl;
                    else if (id == 443) cout << "jpsi" << endl;
                    else if (id == 553) cout << "ups" << endl;
                    else cout << "?"  << endl;
               }
                else if (element == "model") {
                    if (tbl.dipoleModelType() == bSat) cout << "bSat"  << endl;
                    else if (tbl.dipoleModelType() == bNonSat) cout << "bNonSat"  << endl;
                    else if (tbl.dipoleModelType() == bCGC) cout << "bCGC"  << endl;
                    else cout << "?"  << endl;
                }
                else if (element == "pset") {
                    if (tbl.dipoleModelParameterSet() == KMW) cout << "KMW"  << endl;
                    else if (tbl.dipoleModelParameterSet() == HMPZ) cout << "HMPZ"  << endl;
                    else if (tbl.dipoleModelParameterSet() == STU) cout << "STU"  << endl;
                    else if (tbl.dipoleModelParameterSet() == CUSTOM) cout << "CUSTOM"  << endl;
                    else cout << "?"  << endl;
                }
                else if (element == "type") {
                    if (tbl.isUPC()) cout << "UPC"  << endl;
                    else cout << "EIC"  << endl;
                }
                else {
                    cout << "Invalid request\n" << endl;
                    usage(argv[0]);
                    return 2;
                }
            }
        }
        return 0;
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
