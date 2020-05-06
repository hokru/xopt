cmake_minimum_required(VERSION 3.9)


#--- OPTIONS
option(OPENMP "ON - enables adding of OpenMP targets, OFF - disables adding of OpenMP targets" ON)
option(INTEL "Turn ON/OFF intel compiler" OFF)
option(GNU "Turn ON/OFF gcc compiler" OFF)
option(PGI "Turn ON/OFF gcc compiler" OFF)

if(INTEL)
  set(CMAKE_Fortran_COMPILER ifort)
elseif(GNU)
  set(CMAKE_Fortran_COMPILER gfortran)
elseif(PGI)
  set(CMAKE_Fortran_COMPILER pgfortran)
endif()


project(xopt Fortran)
set(xopt_URL "https://github.com/hokru/xopt")
set(xopt_LICENSE "GNU Lesser General Public License, version 3 (LGPL-3.0)")
set(xopt_DESCRIPTION "Open-Source External Optimizer for Quantum Chemistry")






list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(custom_color_messages)
set(CMAKE_Fortran_FORMAT FREE)


set(TargetOpenMP_FIND_COMPONENTS "Fortran")
find_package(TargetLAPACK CONFIG REQUIRED)
get_property(_ill TARGET tgt::lapk PROPERTY INTERFACE_LINK_LIBRARIES)
list(GET _ill 0 _ill0)
get_property(_cd TARGET tgt::lapack PROPERTY INTERFACE_COMPILE_DEFINITIONS)
if(${_cd} MATCHES "USING_LAPACK_MKL")
    set(_isMKL " MKL")
endif()
get_property(_illb TARGET tgt::blas PROPERTY INTERFACE_LINK_LIBRARIES)
list(APPEND _addons ${_ill} ${_illb})
message(STATUS "${Cyan}Using LAPACK${_isMKL}${ColourReset}: ${_ill0};...")

# find_package(LAPACK REQUIRED)
# find_package(BLAS REQUIRED)


#--- Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/distclean.cmake
)



# Populate version.dat
# set(version ${PROJECT_VERSION})
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
# string(TIMESTAMP DATEs "%d.%m.%Y %S%H")
string(TIMESTAMP DATE "%H:%M:%S %Y-%m-%d ")
# string(TIMESTAMP DATE)
set(USER $ENV{USER})
set(HOST $ENV{HOSTNAME})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/version.dat.in ${CMAKE_CURRENT_SOURCE_DIR}/src/version.dat)



if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # set(dialect " -ffree-line-length-none ")
  set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # set(dialect "-xHost ")
  set(bounds "-check bounds")
endif()



# Have the .mod files placed in a special directory
set(mod-dir ${CMAKE_BINARY_DIR}/mod)
file(MAKE_DIRECTORY ${mod-dir})
set(CMAKE_Fortran_MODULE_DIRECTORY ${mod-dir})

message("BLAS_LIBRARIES ${BLAS_LIBRARIES}")
message("LAPACK_LIBRARIES ${LAPACK_LIBRARIES}")

FIND_PACKAGE(LAPACK REQUIRED)



add_executable(xopt "")
add_subdirectory(src)
target_link_libraries(xopt
  tgt::lapack)


#--- Install 
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "/usr/local/xopt" CACHE PATH "Install path" FORCE)
endif()
message(STATUS "xopt install: ${CMAKE_INSTALL_PREFIX}")

set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage)
message(STATUS "xopt staging: ${STAGED_INSTALL_PREFIX}")

install(DIRECTORY ${STAGED_INSTALL_PREFIX}/
        DESTINATION ${CMAKE_INSTALL_PREFIX}
        USE_SOURCE_PERMISSIONS)