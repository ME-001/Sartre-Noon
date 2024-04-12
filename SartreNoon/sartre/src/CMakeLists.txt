#===============================================================================
#  CMakeLists.txt (src) 
# 
#  Copyright (C) 2010-2019 Tobias Toll and Thomas Ullrich
# 
#  This file is part of Sartre.
# 
#  This program is free software: you can redistribute it and/or modify 
#  it under the terms of the GNU General Public License as published by 
#  the Free Software Foundation.   
#  This program is distributed in the hope that it will be useful, 
#  but without any warranty; without even the implied warranty of 
#  merchantability or fitness for a particular purpose. See the 
#  GNU General Public License for more details. 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#   
#  Author: Thomas Ullrich 
#  Last update: 
#  $Date: 2019-04-01 14:39:33 -0400 (Mon, 01 Apr 2019) $
#  $Author: ullrich $
#===============================================================================
include(ExternalProject)
cmake_minimum_required (VERSION 3.1)

#
#  Compiler flags for release and debug version
#
set(CMAKE_C_FLAGS "-W")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O")
message("COMPILER = ${CMAKE_CXX_COMPILER_ID}")
set(CMAKE_CXX_FLAGS "-W -Wall -Wextra -pedantic -Wno-long-long")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-potentially-evaluated-expression")
endif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


#
#  Find external required packages
#  (see also FindGSL.cmake and FindROOT.cmke in cmake/modules)
#  Herer we need only the root library to create the table
#  tools.
#

# GSL
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

# ROOT
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})
set(LIBS ${LIBS} ${ROOT_LIBRARIES})

#BOOST
if (MULTITHREADED)
   set(Boost_USE_STATIC_LIBS ON)
   set(Boost_USE_MULTITHREADED ON)
   find_package(Boost 1.39 COMPONENTS thread REQUIRED)
   if(Boost_FOUND)
       include_directories(${Boost_INCLUDE_DIR})
   endif(Boost_FOUND)
endif (MULTITHREADED)

#
#  Include files from sartre package
#
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/gemini)
include_directories(${PROJECT_SOURCE_DIR}/cuba)

#
#  Cuba is built using the autoconf shipped with it
#
ExternalProject_Add(
    cuba
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/cuba ${PROJECT_BINARY_DIR}/cuba
    SOURCE_DIR ${PROJECT_BINARY_DIR}/cuba
    PREFIX ${PROJECT_BINARY_DIR}/cuba
    CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/cuba/configure --prefix=${PROJECT_BINARY_DIR}/cuba/build
    BUILD_COMMAND make lib
    BUILD_IN_SOURCE 1
)

#
#  Defines source files for sartre library
#
set(SARTRE_SRC "AlphaStrong.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Amplitudes.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "BreakupProduct.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "CrossSection.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "DglapEvolution.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "DipoleModel.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "DipoleModelParameters.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Event.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "ExclusiveFinalStateGenerator.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "FinalStateGenerator.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "FrangibleNucleus.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "EventGeneratorSettings.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Integrals.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Kinematics.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "ModeFinderFunctor.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Nucleon.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Nucleus.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "PhotonFlux.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Sartre.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Settings.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Table.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "TableCollection.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "TableGeneratorNucleus.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "TableGeneratorSettings.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "WaveOverlap.cpp")
add_library(sartre ${SARTRE_SRC})

#
#  Table tools (all stand alone programs) 
#
add_executable(tableInspector tableInspector.cpp)
add_executable(tableMerger tableMerger.cpp)
add_executable(tableQuery tableQuery.cpp)
add_executable(tableDumper tableDumper.cpp)
add_executable(tableVarianceMaker tableVarianceMaker.cpp)
target_link_libraries(tableInspector sartre ${LIBS})
target_link_libraries(tableMerger sartre ${LIBS})
target_link_libraries(tableQuery sartre ${LIBS})
target_link_libraries(tableDumper sartre ${LIBS})
target_link_libraries(tableVarianceMaker sartre ${LIBS})
add_dependencies(tableInspector sartre)
add_dependencies(tableMerger sartre)
add_dependencies(tableQuery sartre)
add_dependencies(tableDumper sartre)
add_dependencies(tableVarianceMaker sartre)

#
#  Install library and include files (make install) within
#  the distribution tree. Top level CMakeLists.txt will install
#  Sartre in final destination.
#
install(TARGETS sartre DESTINATION sartre/lib)
install(FILES "${PROJECT_BINARY_DIR}/cuba/build/lib/libcuba.a" DESTINATION sartre/lib)

FILE(GLOB AllIncludeFiles *.h)
install(FILES ${AllIncludeFiles} DESTINATION sartre/include)
install(FILES "${PROJECT_BINARY_DIR}/cuba/build/include/cuba.h" DESTINATION sartre/include)

install(TARGETS tableInspector DESTINATION sartre/bin)
install(TARGETS tableMerger DESTINATION sartre/bin)
install(TARGETS tableQuery DESTINATION sartre/bin)
install(TARGETS tableDumper DESTINATION sartre/bin)
install(TARGETS tableVarianceMaker DESTINATION sartre/bin)
