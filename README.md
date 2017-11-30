# Introduction

These are the requisite libraries for optimizing the attenutation t\* and depth estimation parameter for primary waveforms observed at teleseismic distances.  This library was written for CTBTO's Expert Technical Analysis group by Instrumental Software Technology's Inc.

# Installation Requirements

The following list is the requisite software libraries for building tdsearch.

  - A C 2011 compliant compiler
  - A Fortran90 compiler
  - [CMake](https://cmake.org/) for generation of Makefiles.
  - [Message Passing Interface](https://www.open-mpi.org/).
  - [Intel Performance Primitives](https://software.intel.com/en-us/intel-ipp).
  - [LAPACK(E)](http://www.netlib.org/lapack/) and [(C)BLAS](http://www.netlib.org/blas/).  If using Intel then the [MKL](https://software.intel.com/en-us/mkl) should be used instead of the libraries available through a package manager.  If using IBM Power architecture then the [ESSL](https://www-03.ibm.com/systems/power/software/essl/) library should be used.
  - The high-performance file IO library [HDF5](https://support.hdfgroup.org/HDF5/).
  - The geodetic computation [GeographicLib](https://geographiclib.sourceforge.io/).
  - ISTI's [libiscl](https://github.com/bakerb845/libiscl).
  - The [ttimes](https://github.com/bakerb845/libttimes) for theoretical travel-times.
  - ISTI's [SAC IO](https://github.com/bakerb845/sacio).
  - ISTI's signals processing library [ISPL](https://github.com/bakerb845/ispl/).
  - The library version of Programs in Seismology's waveform modeling routines [libcps](https://github.com/bakerb845/libcps).
  - The initialization file parsing library [iniparser](https://github.com/ndevilla/iniparser).
  - The C moment tensor manipulation library [compearth](https://github.com/bakerb845/compearth).
  - The Parallel Moment Tensor estimation library [parmt](https://github.com/bakerb845/parmt).

# Configuring CMake

An example CMake compilation script is given by

    #!/bin/sh
    cmake ./ -DCMAKE_BUILD_TYPE=DEBUG \
    -DCMAKE_INSTALL_PREFIX=./ \
    -DCMAKE_C_FLAGS="-g3 -O2 -fopenmp" \
    -DMPI_C_INCLUDE_PATH=/path/to/mpi/include \
    -DMPI_C_LIBRARIES=/path/to/mpi/lib/libmpi.so \
    -DTDSEARCH_USE_INTEL=TRUE \
    -DCOMPEARTH_LIBRARY=/path/to/compearth/lib/libcompearth_shared.so \
    -DCOMPEARTH_INCLUDE_DIR=/path/to/compearth/include \
    -DPARMT_LIBRARY="/path/to/parmt/lib//libparmt_shared.so;/path/to/parmt//lib/libparmtUtils_shared.so;/path/to/parmt/lib/libprepmt_shared.so" \
    -DPARMT_INCLUDE_DIR=/path/to/parmt/include \
    -DMKL_LIBRARY="/path/to/mkl/libmkl_intel_lp64.a;/path/to/mkl/lib/libmkl_sequential.a;/path/to/mkl/lib/libmkl_core.a" \
    -DIPP_LIBRARY="/path/to/ipp/lib/libipps.a;/path/to/ipp/lib/libippvm.a;/path/to/ipp/lib/libippcore.a" \
    -DISCL_LIBRARY=/path/to/libiscl/lib/libiscl_shared.so \
    -DISCL_INCLUDE_DIR=/path/to/lib/iscl/include \
    -DISPL_LIBRARY=/path/to/ispl/lib/libispl_shared.so \
    -DISPL_INCLUDE_DIR=/path/to/ispl/include \
    -DCPS_LIBRARY=/path/to/cps/lib/libcps_shared.so \
    -DCPS_INCLUDE_DIR=/path/to/cps/include \
    -DTTIMES_LIBRARY=/path/to/libttimes/lib/libttimes_shared.so \
    -DTTIMES_INCLUDE_DIR=/path/to/libttimes/include \
    -DSACIO_INCLUDE_DIR=/path/to/sacio/include \
    -DSACIO_LIBRARY=/path/to/sacio/lib/libsacio_shared.so \
    -DH5_C_INCLUDE_DIR=/path/to/hdf5/include \
    -DH5_LIBRARY=/path/to/hdf5/lib/libhdf5.so \
    -DINIPARSER_INCLUDE_DIR=/path/to/iniparser/src \
    -DINIPARSER_LIBRARY=/path/to/iniparser/libiniparser.a

and can be run in the root source directory by

    ./config.sh

If the Doxygen is available, then the Doxygen documentation can be obtained by typing

    doxygen Doxyfile

in the root source directory.

