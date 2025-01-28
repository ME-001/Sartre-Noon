# Install script for directory: /home/bhakta/Sartre-Noon/sartre/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/bhakta/Sartre-Noon/sartrei")
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

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/lib" TYPE STATIC_LIBRARY FILES "/home/bhakta/Sartre-Noon/sartreb/src/libsartre.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/bhakta/Sartre-Noon/sartreb/src/CMakeFiles/sartre.dir/install-cxx-module-bmi-release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/lib" TYPE FILE FILES "/home/bhakta/Sartre-Noon/sartreb/cuba/build/lib/libcuba.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/include" TYPE FILE FILES
    "/home/bhakta/Sartre-Noon/sartre/src/._AlphaStrong.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._Amplitudes.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._BreakupProduct.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._Constants.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._CrossSection.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._DglapEvolution.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._DipoleModel.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._DipoleModelParameters.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._EicSmearFormatWriter.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._Enumerations.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._Event.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._EventGeneratorSettings.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._ExclusiveFinalStateGenerator.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._FinalStateGenerator.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._FrangibleNucleus.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._GridSpline.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._Integrals.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._Kinematics.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._ModeFinderFunctor.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._Nucleon.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._Settings.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._Table.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._TableCollection.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._TableGeneratorNucleus.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._TableGeneratorSettings.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._TwoBodyVectorMesonDecay.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._VectorMesonDecayMass.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._Version.h"
    "/home/bhakta/Sartre-Noon/sartre/src/._WaveOverlap.h"
    "/home/bhakta/Sartre-Noon/sartre/src/AlphaStrong.h"
    "/home/bhakta/Sartre-Noon/sartre/src/Amplitudes.h"
    "/home/bhakta/Sartre-Noon/sartre/src/BreakupProduct.h"
    "/home/bhakta/Sartre-Noon/sartre/src/Constants.h"
    "/home/bhakta/Sartre-Noon/sartre/src/CrossSection.h"
    "/home/bhakta/Sartre-Noon/sartre/src/DglapEvolution.h"
    "/home/bhakta/Sartre-Noon/sartre/src/DipoleModel.h"
    "/home/bhakta/Sartre-Noon/sartre/src/DipoleModelParameters.h"
    "/home/bhakta/Sartre-Noon/sartre/src/EicSmearFormatWriter.h"
    "/home/bhakta/Sartre-Noon/sartre/src/Enumerations.h"
    "/home/bhakta/Sartre-Noon/sartre/src/Event.h"
    "/home/bhakta/Sartre-Noon/sartre/src/EventGeneratorSettings.h"
    "/home/bhakta/Sartre-Noon/sartre/src/ExclusiveFinalStateGenerator.h"
    "/home/bhakta/Sartre-Noon/sartre/src/FinalStateGenerator.h"
    "/home/bhakta/Sartre-Noon/sartre/src/FrangibleNucleus.h"
    "/home/bhakta/Sartre-Noon/sartre/src/GridSpline.h"
    "/home/bhakta/Sartre-Noon/sartre/src/GridSplineBetaPolynomials.h"
    "/home/bhakta/Sartre-Noon/sartre/src/GridSplineInterpolateSums.h"
    "/home/bhakta/Sartre-Noon/sartre/src/Integrals.h"
    "/home/bhakta/Sartre-Noon/sartre/src/IntegrandWrappers.h"
    "/home/bhakta/Sartre-Noon/sartre/src/Kinematics.h"
    "/home/bhakta/Sartre-Noon/sartre/src/ModeFinderFunctor.h"
    "/home/bhakta/Sartre-Noon/sartre/src/Nucleon.h"
    "/home/bhakta/Sartre-Noon/sartre/src/Nucleus.h"
    "/home/bhakta/Sartre-Noon/sartre/src/PhotonFlux.h"
    "/home/bhakta/Sartre-Noon/sartre/src/Sartre.h"
    "/home/bhakta/Sartre-Noon/sartre/src/Settings.h"
    "/home/bhakta/Sartre-Noon/sartre/src/Table.h"
    "/home/bhakta/Sartre-Noon/sartre/src/TableCollection.h"
    "/home/bhakta/Sartre-Noon/sartre/src/TableGeneratorNucleus.h"
    "/home/bhakta/Sartre-Noon/sartre/src/TableGeneratorSettings.h"
    "/home/bhakta/Sartre-Noon/sartre/src/TwoBodyVectorMesonDecay.h"
    "/home/bhakta/Sartre-Noon/sartre/src/VectorMesonDecayMass.h"
    "/home/bhakta/Sartre-Noon/sartre/src/Version.h"
    "/home/bhakta/Sartre-Noon/sartre/src/WaveOverlap.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/include" TYPE FILE FILES "/home/bhakta/Sartre-Noon/sartreb/cuba/build/include/cuba.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableInspector" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableInspector")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableInspector"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/bin" TYPE EXECUTABLE FILES "/home/bhakta/Sartre-Noon/sartreb/src/tableInspector")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableInspector" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableInspector")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableInspector")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/bhakta/Sartre-Noon/sartreb/src/CMakeFiles/tableInspector.dir/install-cxx-module-bmi-release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableMerger" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableMerger")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableMerger"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/bin" TYPE EXECUTABLE FILES "/home/bhakta/Sartre-Noon/sartreb/src/tableMerger")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableMerger" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableMerger")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableMerger")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/bhakta/Sartre-Noon/sartreb/src/CMakeFiles/tableMerger.dir/install-cxx-module-bmi-release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableQuery" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableQuery")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableQuery"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/bin" TYPE EXECUTABLE FILES "/home/bhakta/Sartre-Noon/sartreb/src/tableQuery")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableQuery" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableQuery")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableQuery")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/bhakta/Sartre-Noon/sartreb/src/CMakeFiles/tableQuery.dir/install-cxx-module-bmi-release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableDumper" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableDumper")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableDumper"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/bin" TYPE EXECUTABLE FILES "/home/bhakta/Sartre-Noon/sartreb/src/tableDumper")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableDumper" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableDumper")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableDumper")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/bhakta/Sartre-Noon/sartreb/src/CMakeFiles/tableDumper.dir/install-cxx-module-bmi-release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableVarianceMaker" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableVarianceMaker")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableVarianceMaker"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/bin" TYPE EXECUTABLE FILES "/home/bhakta/Sartre-Noon/sartreb/src/tableVarianceMaker")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableVarianceMaker" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableVarianceMaker")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tableVarianceMaker")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/bhakta/Sartre-Noon/sartreb/src/CMakeFiles/tableVarianceMaker.dir/install-cxx-module-bmi-release.cmake" OPTIONAL)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tablePriority" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tablePriority")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tablePriority"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/sartre/bin" TYPE EXECUTABLE FILES "/home/bhakta/Sartre-Noon/sartreb/src/tablePriority")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tablePriority" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tablePriority")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/sartre/bin/tablePriority")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  include("/home/bhakta/Sartre-Noon/sartreb/src/CMakeFiles/tablePriority.dir/install-cxx-module-bmi-release.cmake" OPTIONAL)
endif()

