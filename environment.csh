#!/bin/tcsh

########################################################################
# PREFIX and non-system programs/libraries
########################################################################

### prefix area
setenv PREFIX /group/clas/builds

### non-system builds of programs and libraries
setenv GCC /apps/gcc/4.8.0
#setenv ROOT /apps/root/5.34.05
setenv CERN /apps/cernlib/x86_64_rhel6_4.7.2/2005
setenv PYTHON /apps/python
setenv SCONS /apps/scons
setenv BOOST /group/clas/boost/boost-1.53.0

########################################################################
# PATH
########################################################################

setenv PATH .:${PREFIX}/bin
setenv PATH ${PATH}:${PREFIX}/scripts

setenv PATH ${PATH}:${GCC}/bin
#setenv PATH ${PATH}:${ROOT}/root/bin
setenv PATH ${PATH}:${PYTHON}/bin
setenv PATH ${PATH}:${SCONS}/bin

### standard system paths
setenv PATH ${PATH}:/site/bin:/apps/bin
setenv PATH ${PATH}:/usr/bin:/bin:/usr/sbin:/sbin

setenv PATH ${PATH}:./bin:./build/bin

########################################################################
# LD_LIBRARY_PATH
########################################################################

### run-time library loading path
setenv LD_LIBRARY_PATH .:${PREFIX}/lib

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${GCC}/lib64
#setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${ROOT}/root/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PYTHON}/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${SCONS}/lib

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${BOOST}/lib

########################################################################
# PYTHONPATH
########################################################################

### python modules search path

setenv PYTHONPATH ${PREFIX}/lib/python

########################################################################
# sources for build directories
########################################################################

setenv MYSQLINC /usr/include/mysql
setenv MYSQLLIB /usr/lib64/mysql

setenv BOOSTINC ${BOOST}
setenv BOOSTLIB ${BOOST}/lib

setenv CERNLIB ${CERN}/lib

setenv CLAS6INC ${PREFIX}/include
setenv CLAS6LIB ${PREFIX}/lib

########################################################################
# misc
########################################################################

setenv CLAS_PARMS /group/clas/parms
setenv JLAB_ROOT /site/12gev_phys
setenv JLAB_VERSION 1.3
source $JLAB_ROOT/$JLAB_VERSION/ce/jlab.csh

rehash

