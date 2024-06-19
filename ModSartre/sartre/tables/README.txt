//==============================================================================
//  tables/README.txt
//
//  Copyright (C) 2019 Tobias Toll and Thomas Ullrich
//
//  This file is part of Sartre.
//
//  Author: Thomas Ullrich
//  Last update: Sep 6, 2019, TU
//==============================================================================

For event generation Sartre needs a set of tables that are stored
in a specific directory structure under sartre/tables/.
Tables are stored as binary ROOT files. Each file contains one table.
Sartre itself does not use the names but queries the content
of the files. Only requirement is that they have the extension .root.
What matters for Sartre is the directory structure under which the
files are stored. The naming convention described here is mostly for
bookkeeping and easy searching.

To-date every table ever created is stored under sartre/tables/. Some
are not complete. If you are not sure what set you need and if the
available one fits your needs, ask the authors. For now, the most complete
set is the one done with the KMW parameter set. The HMPZ data sets for
protons are complete, none for any A>1 exist. The Pb (208) tables for
EIC should not be used. Most UPC tables are for proton and Pb and are
created with either the STU or the HMPZ parameter set.

The following naming convention for Sartre tables is used:
<TYPE>-A<A>-<MOD>-<VM>-<PSET>-<CONT>-<POL>.root
<TYPE>-A<A>-<SAT>-<VM>-<PSET>-<CONT>-<POL>-C=<OPT>.root

where:
<TYPE>     EIC or UPC
           EIC = Electron-Ion Collisions
           UPC = UltraPeripheral Collisions
<A>        is the mass number, 1 for protons
<MOD>      dipole model, either bSat, bNonSat, or bCGC
<VM>       vector meson name (dvcs, rho, phi, jpsi, ups)
<PSET>     parameter set: one of KMW, HMPZ, STU
<CONT>     content of the table, one of
           mean_A2
           mean_A
           lambda_real
           lambda_skew
           variance_A
<POL>      polarization: L or T
<OPT>      option string, options are separated by dashes "-"

Examples:
EIC-A197-bSat-jpsi-KMW-lambda_real-L.root
UPC-A1-bSat-jpsi-STU-variance_A-T.root
EIC-A197-bNonSat-dvcs-KMW-mean_A2-T-C=region2.root
UPC-A208-bNonSat-jpsi-STU-variance_A-T-C=region2-despiked.root
