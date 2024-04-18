# Install script for directory: /home/bhakta/sartre/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/bhakta/sartre/sartreinstall")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/lib" TYPE STATIC_LIBRARY FILES "/home/bhakta/sartre/sartreinstall/src/libsartre.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/lib" TYPE FILE FILES "/home/bhakta/sartre/sartreinstall/cuba/build/lib/libcuba.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/include" TYPE FILE FILES
    "/home/bhakta/sartre/src/._Event.h"
    "/home/bhakta/sartre/src/._Integrals.h"
    "/home/bhakta/sartre/src/._PhotonFlux.h"
    "/home/bhakta/sartre/src/._Table.h"
    "/home/bhakta/sartre/src/._TableCollection.h"
    "/home/bhakta/sartre/src/._Version.h"
    "/home/bhakta/sartre/src/AlphaStrong.h"
    "/home/bhakta/sartre/src/Amplitudes.h"
    "/home/bhakta/sartre/src/BreakupProduct.h"
    "/home/bhakta/sartre/src/Constants.h"
    "/home/bhakta/sartre/src/CrossSection.h"
    "/home/bhakta/sartre/src/DglapEvolution.h"
    "/home/bhakta/sartre/src/DipoleModel.h"
    "/home/bhakta/sartre/src/DipoleModelParameters.h"
    "/home/bhakta/sartre/src/Enumerations.h"
    "/home/bhakta/sartre/src/Event.h"
    "/home/bhakta/sartre/src/EventGeneratorSettings.h"
    "/home/bhakta/sartre/src/ExclusiveFinalStateGenerator.h"
    "/home/bhakta/sartre/src/FinalStateGenerator.h"
    "/home/bhakta/sartre/src/FrangibleNucleus.h"
    "/home/bhakta/sartre/src/GridSpline.h"
    "/home/bhakta/sartre/src/GridSplineBetaPolynomials.h"
    "/home/bhakta/sartre/src/GridSplineInterpolateSums.h"
    "/home/bhakta/sartre/src/Integrals.h"
    "/home/bhakta/sartre/src/IntegrandWrappers.h"
    "/home/bhakta/sartre/src/Kinematics.h"
    "/home/bhakta/sartre/src/ModeFinderFunctor.h"
    "/home/bhakta/sartre/src/NeutronGenerator.h"
    "/home/bhakta/sartre/src/Nucleon.h"
    "/home/bhakta/sartre/src/Nucleus.h"
    "/home/bhakta/sartre/src/PhotonFlux.h"
    "/home/bhakta/sartre/src/Sartre.h"
    "/home/bhakta/sartre/src/Settings.h"
    "/home/bhakta/sartre/src/Table.h"
    "/home/bhakta/sartre/src/TableCollection.h"
    "/home/bhakta/sartre/src/TableGeneratorNucleus.h"
    "/home/bhakta/sartre/src/TableGeneratorSettings.h"
    "/home/bhakta/sartre/src/Version.h"
    "/home/bhakta/sartre/src/WaveOverlap.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/include" TYPE FILE FILES "/home/bhakta/sartre/sartreinstall/cuba/build/include/cuba.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableInspector" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableInspector")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableInspector"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/bin" TYPE EXECUTABLE FILES "/home/bhakta/sartre/sartreinstall/src/tableInspector")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableInspector" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableInspector")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableInspector")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableMerger" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableMerger")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableMerger"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/bin" TYPE EXECUTABLE FILES "/home/bhakta/sartre/sartreinstall/src/tableMerger")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableMerger" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableMerger")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableMerger")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableQuery" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableQuery")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableQuery"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/bin" TYPE EXECUTABLE FILES "/home/bhakta/sartre/sartreinstall/src/tableQuery")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableQuery" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableQuery")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableQuery")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableDumper" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableDumper")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableDumper"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/bin" TYPE EXECUTABLE FILES "/home/bhakta/sartre/sartreinstall/src/tableDumper")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableDumper" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableDumper")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableDumper")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableVarianceMaker" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableVarianceMaker")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableVarianceMaker"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/bin" TYPE EXECUTABLE FILES "/home/bhakta/sartre/sartreinstall/src/tableVarianceMaker")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableVarianceMaker" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableVarianceMaker")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableVarianceMaker")
    endif()
  endif()
endif()

