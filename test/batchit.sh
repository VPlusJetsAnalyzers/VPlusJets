#!/bin/sh

SRCDIR=$1
shift
cd $SRCDIR
source /uscmst1/prod/sw/cms/shrc cvmfs
eval `scram runtime -sh`
cd -
uname -n
echo "$*"
eval "$*"
