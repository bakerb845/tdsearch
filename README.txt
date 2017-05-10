These are the requisite libraries for optimizing the attenutation t* and
depth estimation parameter for primary waveforms observed at teleseismic
distances.  This library was written for CTBTO's Expert Technical Analysis
group by Instrumental Software Technology's Inc.

Installation Requirements:

 (1)  A C 2011 compliant compiler
 (2)  A Fortran90 compiler
 (3)  CMake
        https://cmake.org/
 (4)  Message Passing Interface:
        https://www.open-mpi.org/
 (5)  HDF5 
        https://support.hdfgroup.org/HDF5/
 (6)  GeographicLib
        https://geographiclib.sourceforge.io/
 (7)  ISTI's libiscl
        https://github.com/bakerb845/libiscl
 (8)  The ttimes library for theoretical travel time compuations 
        https://github.com/bakerb845/libttimes
 (9)  ISTI's sacio library
        https://github.com/bakerb845/sacio
 (10) ISTI's signals processing library
        https://github.com/bakerb845/ispl/
 (11) The library version of Programs in Seismology's waveform modeling routines
        https://github.com/bakerb845/libcps
 (12) The iniparser library
        https://github.com/ndevilla/iniparser
 (13) The parmt library (currently contains compearth)
        https://github.com/bakerb845/parmt

Recommended libraries:
 (1) Intel MKL and Performance Primitives:
       https://software.intel.com/en-us/performance-libraries

As an example: I initialize cmake with the following script:

#!/bin/sh
cmake ./ -DCMAKE_BUILD_TYPE=DEBUG \
-DCMAKE_INSTALL_PREFIX=./ \
-DCMAKE_C_COMPILER=icc \
-DCMAKE_CXX_COMPILER=icpc \
-DCMAKE_C_FLAGS="-g3 -O2 -qopenmp -Wall -Wextra -Wcomment -Wcheck" \
-DCMAKE_CXX_FLAGS="-g3 -O2 -Wall -qopenmp" \
-DMPI_C=/opt/intel/impi/2017.1.132/bin64/mpiicc \
-DMPI_C_INCLUDE_PATH=/opt/intel/impi/2017.1.132/include64 \
-DMPI_C_LIBRARIES=/opt/intel/impi/2017.1.132/lib64/libmpi.so \
-DTDSEARCH_USE_INTEL=TRUE \
-DCOMPEARTH_LIBRARY=/home/bakerb25/C/parmt/lib/libcompearth_shared.so \
-DCOMPEARTH_INCLUDE_DIR=/home/bakerb25/C/parmt/compearth \
-DPARMT_LIBRARY="/home/bakerb25/C/parmt/lib/libparmt_shared.so;/home/bakerb25/C/parmt/lib/libparmtUtils_shared.so" \
-DPARMT_INCLUDE_DIR=/home/bakerb25/C/parmt/include \
-DMKL_LIBRARY="/opt/intel/mkl/lib/intel64_lin/libmkl_intel_lp64.a;/opt/intel/mkl/lib/intel64_lin/libmkl_sequential.a;/opt/intel/mkl/lib/intel64_lin/libmkl_core.a" \
-DIPP_LIBRARY="/opt/intel/ipp/lib/intel64_lin/libipps.a;/opt/intel/ipp/lib/intel64_lin/libippvm.a;/opt/intel/ipp/lib/intel64_lin/libippcore.a" \
-DIPP_INCLUDE_DIR=/opt/intel/ipp/include \
-DMKL_INCLUDE_DIR=/opt/intel/mkl/include \
-DISCL_LIBRARY=/home/bakerb25/C/libiscl/libiscl_shared.so \
-DISCL_INCLUDE_DIR=/home/bakerb25/C/libiscl/include \
-DISPL_LIBRARY=/home/bakerb25/C/ispl/libispl_shared.so \
-DISPL_INCLUDE_DIR=/home/bakerb25/C/ispl/include \
-DCPS_LIBRARY=/home/bakerb25/C/libcps/libcps_shared.so \
-DCPS_INCLUDE_DIR=/home/bakerb25/C/libcps/include \
-DTTIMES_LIBRARY=/home/bakerb25/C/libttimes/lib/libttimes_shared.so \
-DTTIMES_INCLUDE_DIR=/home/bakerb25/C/libttimes/include \
-DSACIO_INCLUDE_DIR=/home/bakerb25/C/sacio/include \
-DSACIO_LIBRARY=/home/bakerb25/C/sacio/libsacio_shared.so \
-DH5_C_INCLUDE_DIR=/home/bakerb25/C/hdf5-1.10.0_intel/include \
-DH5_LIBRARY=/home/bakerb25/C/hdf5-1.10.0_intel/lib/libhdf5.so \
-DINIPARSER_INCLUDE_DIR=/home/bakerb25/C/iniparser/src \
-DINIPARSER_LIBRARY=/home/bakerb25/C/iniparser/libiniparser.a

--------------------------------------------------------------------------------

In-code documentation can be generated with doxygen by typing:

doxygen Doxyfile

in the root source directory.

