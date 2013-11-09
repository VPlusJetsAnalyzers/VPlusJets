#!/bin/sh

SRCDIR=$1
shift
cd $SRCDIR
#export SCRAM_ARCH=slc5_amd64_gcc462
export SCRAM_ARCH=slc5_amd64_gcc472
eval `scram runtime -sh`
#if [ -e pyroot_logon.py ]
#then
#    python pyroot_logon.py
#fi
cd -
echo "$*"
`$*`
