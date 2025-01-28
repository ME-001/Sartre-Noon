//==============================================================================
//  tablePriority.cpp
//
//  Copyright (C) 2021 Tobias Toll and Thomas Ullrich
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
//  $Date: 2021-08-24 16:17:14 -0400 (Tue, 24 Aug 2021) $
//  $Author: ullrich $
//==============================================================================
//   
//  Utility program to display or set the table priority.
//  Usage:  tableInspector [-s priority] file(s) ...
//          -s priority    set the priority to given value
//==============================================================================
#include "Table.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH3F.h"
#include <iostream>
#include <sstream>
#include <unistd.h>   
#include <cstdlib>  
#include <string>  
#include <unistd.h>
   
using namespace std;   
  
#define PR(x) cout << #x << " = " << (x) << endl;

void usage(const char* prog)   
{   
    cout << "Usage: " << prog << " [-s priority] [-d] file(s) ..." << endl;
    cout << "                    " << "-s priority    set the priority to given value" << endl;
    cout << "                    " << "-f             display filename after priority" << endl;
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
       
    bool setPriority = false;
    bool displayFilename = false;
    unsigned int newPriority = 0;
    int  ch;
    
    while ((ch = getopt(argc, argv, "fs:")) != -1) {
        switch (ch) {   
            case 'f':
                displayFilename = true;
                break;
            case 's':
                setPriority = true;
                newPriority = atol(optarg);
                break;
            case '?':
            default:   
                usage(argv[0]);   
                return 2;   
                break;   
        }   
    }   
    if (optind == argc) {
        usage(argv[0]);   
        return 2;   
    }   
      
    if (setPriority && newPriority > 0xFF) {
        cout << "Error, priority cannot be larger than 255." << endl;
        return 1;
    }

    //
    //  Set priority mode
    //
    TFile *file;
    if (setPriority) {
        for (int index = optind; index < argc; index++) {
            //
            //  Open file in read only mode and get table
            //
            file = TFile::Open(argv[index],"READ");
            if (!file) {
                cout << "Error, failed opening file '" << argv[index] << "'." << endl;
                return 1;
            }
            auto ptr = file->Get("table");
            if (ptr == 0) {
                cout << "Error, failed retrieving table from file '" << argv[index] << "'." << endl;
                return 1;
            }
            
            //
            //  Rewrite ID (histo title) using new priority
            //
            uint64_t mID = atoll(ptr->GetTitle());
            uint64_t one = 1;
            int oldPriority = ((mID >> 34) & 0xFF);
            for (int k=34; k<41; k++) mID &= ~(one << k);
            mID |= (static_cast<uint64_t>(newPriority) << 34);
            ostringstream titlestream;
            titlestream << mID;
            string title = titlestream.str();
            
            //
            //  Type of table different for UPC
            //
            bool isUPCTable = (mID & (static_cast<uint64_t>(1) << 46));
            if (isUPCTable) {
                auto hist = reinterpret_cast<TH2F*>(ptr);
                hist->SetDirectory(0);
                hist->SetTitle(title.c_str());
            }
            else {
                auto hist = reinterpret_cast<TH3F*>(ptr);
                hist->SetDirectory(0);
                hist->SetTitle(title.c_str());
            }

            file->Close();
            
            //
            //  Open same file in recreate/new mode and write
            //  updated histos into them.
            //  We need to write a new file since adding it
            //  to the same one (update mode) doubles the size
            //  of the file otherwise.
            //
            file = TFile::Open(argv[index],"RECREATE");
            if (isUPCTable) {
                auto hist = reinterpret_cast<TH2F*>(ptr);
                hist->Write();
            }
            else {
                auto hist = reinterpret_cast<TH3F*>(ptr);
                hist->Write();
            }
            file->Close();

            //
            //  Print out
            //
            cout << oldPriority << " -> " << newPriority;
            if (displayFilename) cout  << '\t' << argv[index];
            cout << endl;
        }
    }
    
    //
    //  List mode only
    //
    else {
        for (int index = optind; index < argc; index++) {
            Table tbl;
            if (tbl.read(argv[index])) {
                cout << tbl.priority();
                if (displayFilename) cout  << '\t' << tbl.filename();
                cout << endl;
            }
        }
    }
    
    return 0;
}
