# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/bhakta/Sartre-Noon/sartreb/cuba"
  "/home/bhakta/Sartre-Noon/sartreb/cuba/src/cuba-build"
  "/home/bhakta/Sartre-Noon/sartreb/cuba"
  "/home/bhakta/Sartre-Noon/sartreb/cuba/tmp"
  "/home/bhakta/Sartre-Noon/sartreb/cuba/src/cuba-stamp"
  "/home/bhakta/Sartre-Noon/sartreb/cuba/src"
  "/home/bhakta/Sartre-Noon/sartreb/cuba/src/cuba-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/bhakta/Sartre-Noon/sartreb/cuba/src/cuba-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/bhakta/Sartre-Noon/sartreb/cuba/src/cuba-stamp${cfgdir}") # cfgdir has leading slash
endif()
