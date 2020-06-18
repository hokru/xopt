cmake_minimum_required(VERSION 3.9)



#--- OPTIONS
option(OPENMP "ON - enables adding of OpenMP targets, OFF - disables adding of OpenMP targets" ON)
option(INTEL "Turn ON/OFF intel compiler" OFF)
option(GNU "Turn ON/OFF gcc compiler" OFF)
option(PGI "Turn ON/OFF pgi compiler" OFF)
option(CLANG "Turn ON/OFF clang compiler" OFF)
option(Enable_D3 "Use dftd3 lib" OFF)
option(Enable_GCP "Use gcplib" OFF)
option(Enable_LBFGS "Use bfgs lib" OFF)
option(Enable_EXTRAS "Use all 3 above libs" OFF)

if(INTEL)
  set(CMAKE_Fortran_COMPILER ifort)
  list(APPEND D_FLAGS INTEL)
elseif(GNU)
  set(CMAKE_Fortran_COMPILER gfortran)
  list(APPEND D_FLAGS GNU)
elseif(PGI)
  set(CMAKE_Fortran_COMPILER pgfortran)
  list(APPEND D_FLAGS PGI)
#elseif(CLANG)
#   set(CMAKE_Fortran_COMPILER clang)
#   list(APPEND D_FLAGS CLANG)
endif()

if(Enable_D3)
  list(APPEND D_FLAGS DFTD3)
endif()
if(Enable_GCP)
  list(APPEND D_FLAGS GCP)
endif()
if(Enable_LBFGS)
  list(APPEND D_FLAGS GCP)
endif()
if(Enable_EXTRAS)
  list(APPEND D_FLAGS GCP DFTD3 BFGS)
endif()
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# set project data
project(xopt Fortran)
set(xopt_URL "https://github.com/hokru/xopt")
set(xopt_LICENSE "GNU Lesser General Public License, version 3 (LGPL-3.0)")
set(xopt_DESCRIPTION "Open-Source External Optimizer for Quantum Chemistry")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(custom_color_messages)
set(CMAKE_Fortran_FORMAT FREE)


find_package(OpenMP)
if(OpenMP_FC_FOUND)
  target_link_libraries(MyTarget PUBLIC OpenMP::OpenMP_FC)
endif()

# BLAS AND LAPACK DETECTION
set(BLAS "MKL" CACHE STRING "Default Lapack/Blas library")
set_property(CACHE BLAS PROPERTY STRINGS "MKL;OpenBLAS;Generic")
message(STATUS "Searching for Lapack/Blas backend: " ${BLAS})
if(BLAS STREQUAL "MKL")
  find_package(MKL REQUIRED)
  list(APPEND LAPACK_LIBRARIES ${MKL_LIBRARY})
elseif(BLAS STREQUAL "OpenBLAS")
  find_package(OpenBLAS REQUIRED)
  include_directories(LAPACK_INCLUDE_DIR ${OpenBLAS_INCLUDE_DIR})
  list(APPEND LAPACK_LIBRARIES ${OpenBLAS_LIB})
elseif(BLAS STREQUAL "Generic")
  find_package(LAPACK REQUIRED)
endif()


# find_package(BLAS REQUIRED)
#find_package(LAPACK REQUIRED)

message(STATUS "${Cyan}Using LAPACK${ColourReset}: ${LAPACK_LIBRARIES};...")

#--- Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/distclean.cmake
)


# Write version.dat file
execute_process(COMMAND git describe --abbrev=4 --dirty --always --tags
  RESULT_VARIABLE git_return
  OUTPUT_VARIABLE commit)
if(git_return)
  set(GIT_VERSION "unknown")
else()
  string(REGEX REPLACE "\n$" "" GIT_VERSION ${commit})
endif()
execute_process(COMMAND git rev-parse --abbrev-ref HEAD
  RESULT_VARIABLE git_return
  OUTPUT_VARIABLE branch)
if(git_return)
  set(GIT_VERSION "unknown")
else()
  string(REGEX REPLACE "\n$" "" GIT_BRANCH ${branch})
endif()
string(TIMESTAMP DATE "%H:%M:%S %Y-%m-%d ")
set(USER $ENV{USER})
set(HOST $ENV{HOSTNAME})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/version.dat.in ${CMAKE_CURRENT_SOURCE_DIR}/src/version.dat)



# compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-fbounds-check -Og -g  -Wunused-dummy-argument -std=f2008 -fall-intrinsics -fbacktrace")
  # set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-check bounds -O0 -g traceback -std08 -debug -fp-model=consistent -fp-speculation=save ")
  # set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif()


# Have the .mod files placed in a special directory
set(mod-dir ${CMAKE_BINARY_DIR}/mod)
file(MAKE_DIRECTORY ${mod-dir})
set(CMAKE_Fortran_MODULE_DIRECTORY ${mod-dir})


add_executable(xopt "")
add_subdirectory(src)
# target_link_libraries(xopt
#   tgt::lapack)
target_link_libraries(xopt
  ${LAPACK_LIBRARIES})

message(STATUS "${Cyan}pre-processor flags${ColourReset}:  ${D_FLAGS}")
# message(pre-processor flags:  ${D_FLAGS})
# set pre-processor flags  
target_compile_definitions(xopt
    PRIVATE
        ${D_FLAGS}
)  

#--- Install 
# if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
#     set(CMAKE_INSTALL_PREFIX "/usr/local/xopt" CACHE PATH "Install path" FORCE)
# endif()
# message(STATUS "xopt install: ${CMAKE_INSTALL_PREFIX}")

# set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage)
# message(STATUS "xopt staging: ${STAGED_INSTALL_PREFIX}")

install(DIRECTORY ${STAGED_INSTALL_PREFIX}/
        DESTINATION ${CMAKE_INSTALL_PREFIX}
        USE_SOURCE_PERMISSIONS)
