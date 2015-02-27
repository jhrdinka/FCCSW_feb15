#!/bin/sh -u

if [[ "x$GAUDI" == "x" ]]; then
echo "Need to set the GAUDI environment variable to the path of the Gaudi software directory (contains GaudiKernel/)."
return 1
else
echo "Gaudi   :    $GAUDI"
fi
if [[ "x$FCCBASE" == "x" ]]; then
echo "Need to set the FCCBASE environment variable to root path of the software (contains both Gaudi and the FCC software)."
return 1
else
echo "FCC root:    $FCCBASE"
fi
#

export COMOPT=dbg #opt

# set up CMake:
export PATH=/afs/cern.ch/sw/lcg/contrib/CMake/2.8.12.2/Linux-i386/bin:$PATH
#export CMAKEFLAGS='-DCMAKE_USE_CCACHE=ON'
export CMAKE_PREFIX_PATH=$GAUDI/cmake:$FCCBASE:/afs/cern.ch/exp/fcc/sw/0.2
#export CMAKE_PREFIX_PATH=$GAUDI/cmake:$FCCBASE:/afs/cern.ch/sw/lcg/releases:/afs/cern.ch/exp/fcc/sw/0.2/DD4hep
export CMTCONFIG=x86_64-slc6-gcc48-$COMOPT

# set up the compilers
export PATH=/afs/cern.ch/lhcb/software/releases/LBSCRIPTS/LBSCRIPTS_v8r0/InstallArea/scripts:$PATH
export LCG_hostos=x86_64-slc6
#export LCG_external_area=/afs/cern.ch/sw/lcg/external
#export PATH=/afs/cern.ch/sw/lcg/contrib/ninja/1.4.0/x86_64-slc6:$PATH

#export PATH=/afs/cern.ch/user/j/jhrdinka/FCC/GAUDI/GAUDI_v25r2/InstallArea/x86_64-slc6-gcc48-$COMOPT/scripts:$PATH

# rootmap file
export LD_LIBRARY_PATH=/afs/cern.ch/user/j/jhrdinka/FCC/FCCSW/build.x86_64-slc6-gcc48-$COMOPT/lib:$LD_LIBRARY_PATH



