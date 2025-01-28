#=========================================================================
# 
# This cmake is to include the neutron generator part
# It contains th NeutronGenerator.cxx and header files
# 
#
#
# Author: bhakta
# Date  : 13th jan 2025
#=========================================================================

cmake_minimum_required (VERSION 3.1)

#
#  Make sure required external libraries are available
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules")

#
# release or debug and compiler flags
#
set(CMAKE_C_FLAGS_DEBUG "-g -W -Wall -Wextra")
set(CMAKE_C_FLAGS_RELEASE "-O -W -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -W -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O -W -Wall -Wextra")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#
# ROOT
# 
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})

# plan on changing:  not including it and instead use another program to convert the endf txt to root file
# =============================================
#add_custom_command(
#  OUTPUT ${PROJECT_BINARY_DIR}/converted.root
#  COMMAND root ${CMAKE_SOURCE_DIR}/convert_txt_to_root.py ${CMAKE_SOURCE_DIR}/data/input.txt ${PROJECT_BINARY_DIR}/converted.root
#  DEPENDS ${CMAKE_SOURCE_DIR}/data/input.txt
#  COMMENT "Converting .txt to .root"
#)
#==============================================


set(NOON_SRC "NeutronGenerator.cxx")
add_library(noon ${NOON_SRC})


install(TARGETS noon DESTINATION sartre/lib)


FILE(GLOB AllIncludeFiles *.h)
install(FILES ${AllIncludeFiles} DESTINATION sartre/include)

