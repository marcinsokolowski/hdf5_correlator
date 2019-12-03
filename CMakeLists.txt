cmake_minimum_required(VERSION 2.6)

project(hdf5_correlator)

find_library(CFITSIO_LIB cfitsio)
find_path(FITSIO_INCLUDE_DIR NAMES fitsio.h)
find_library(LIBNOVA_LIB nova)
find_library(LIBMSFITS_LIB msfitslib)
# find_library(LIBHDF5 hdf5_cpp)
find_package(HDF5 COMPONENTS C CXX REQUIRED)

find_library(FFTW3_LIB fftw3 REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-noexcept-type -DNDEBUG -march=native -std=c++11 -D_UNIX")

include_directories(${FITSIO_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})


# The following stuff will set the "rpath" correctly, so that
# LD_LIBRARY_PATH doesn't have to be set.

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

add_executable(hdf5_correlator main.cpp hdf5_commons.cpp eda2_tpm_coefficients.cpp calsolutions.cpp)
target_link_libraries(hdf5_correlator ${CFITSIO_LIB} ${LIBNOVA_LIB} ${LIBMSFITS_LIB} ${LIBHDF5} ${HDF5_LIBRARIES} ${LIBNOVA_LIB} ${FFTW3_LIB} ${CFITSIO_LIB})

install (TARGETS hdf5_correlator DESTINATION bin)
