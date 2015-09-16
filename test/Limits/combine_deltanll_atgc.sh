#!/bin/sh

echo $#

if [ $# -lt 2 ]
then
    echo "Usage: $0 wkspace.root par1 [par2 par3]"
    echo "       where parx={lZ|dkg|dg1} and par3 is the parameter to fix to 0"
    exit
fi

WKSPACEROOT=$1

#DEBUGLEV=99
DEBUGLEV=1
#DEBUGLEV=0

RANGE1=$2
PAR1=${RANGE1%=*}
VALS1=${RANGE1#*=}
VALS1=${VALS1/,/ }

if [ $# -ge 4 ]
then
# ranges in form parX=min,max
    RANGE2=$3
    RANGE3=$4

    PAR2=${RANGE2%=*}
    PAR3=${RANGE3%=*}
    VALS2=${RANGE2#*=}
    VALS3=${RANGE3#*=}
    VALS2=${VALS2/,/ }
fi

suffix=${WKSPACEROOT##wkspace_}
suffix=`basename $suffix .root`

#COMMON_ARGS="-M MultiDimFit --redefineSignalPOIs $PAR1,$PAR2 --setPhysicsModelParameters ${RANGE3} --setPhysicsModelParameterRanges ${RANGE1}:${RANGE2} --freezeNuisances ${PAR3} --floatOtherPOIs=0 -v${DEBUGLEV}"

COMMON_ARGS="-M MultiDimFit --floatOtherPOIs=0 -v${DEBUGLEV}"

SCAN_ARGS=" --algo=grid --points=10201 --minimizerStrategy=2"
# -v99"

COMMON_ARGS="$COMMON_ARGS $SCAN_ARGS"

#ONED_ARGS="--algo=singles --robustFit 1 --do95 1"
TWOD_ARGS="-P $PAR1 -P $PAR2 --setPhysicsModelParameterRanges ${RANGE1}:${RANGE2}"

# "-t -1" requests asimov dataset
EXPD_ARGS="-t -1 --expectSignal 1"

########################################
# expected 1D limit:

suffix1d="${suffix}_1Dexp_${PAR1}"
echo "combine $WKSPACEROOT $COMMON_ARGS $EXPD_ARGS -P $PAR1 -n _${suffix1d}"        2>&1 | tee    limit_${suffix1d}.log
combine       $WKSPACEROOT $COMMON_ARGS $EXPD_ARGS -P $PAR1 -n _${suffix1d}         2>&1 | tee -a limit_${suffix1d}.log
echo "python build1DInterval.py ${VALS1} higgsCombine_${suffix1d}*.root ${PAR1} -b" 2>&1 | tee -a limit_${suffix1d}.log
python build1DInterval.py ${VALS1} higgsCombine_${suffix1d}*.root ${PAR1} -b        2>&1 | tee -a limit_${suffix1d}.log

if [ "$PAR2" != "" ]
then
    suffix1d="${suffix}_1Dexp_${PAR2}"
    echo "combine $WKSPACEROOT $COMMON_ARGS $EXPD_ARGS -P $PAR2 -n _${suffix1d}"        2>&1 | tee    limit_${suffix1d}.log
    combine       $WKSPACEROOT $COMMON_ARGS $EXPD_ARGS -P $PAR2 -n _${suffix1d}         2>&1 | tee -a limit_${suffix1d}.log
    echo "python build1DInterval.py ${VALS2} higgsCombine_${suffix1d}*.root ${PAR2} -b" 2>&1 | tee -a limit_${suffix1d}.log
    python build1DInterval.py ${VALS2} higgsCombine_${suffix1d}*.root ${PAR2} -b        2>&1 | tee -a limit_${suffix1d}.log
fi

########################################
# observed 1D limit:

suffix1d="${suffix}_1Dobs_${PAR1}"
echo "combine $WKSPACEROOT $COMMON_ARGS -P $PAR1 -n _${suffix1d}"                   2>&1 | tee    limit_${suffix1d}.log
combine       $WKSPACEROOT $COMMON_ARGS -P $PAR1 -n _${suffix1d}                    2>&1 | tee -a limit_${suffix1d}.log
echo "python build1DInterval.py ${VALS1} higgsCombine_${suffix1d}*.root ${PAR1} -b" 2>&1 | tee -a limit_${suffix1d}.log
python build1DInterval.py ${VALS1} higgsCombine_${suffix1d}*.root ${PAR1} -b        2>&1 | tee -a limit_${suffix1d}.log

if [ "$PAR2" != "" ]
then
    suffix1d="${suffix}_1Dobs_${PAR2}"
    echo "combine $WKSPACEROOT $COMMON_ARGS -P $PAR2 -n _${suffix1d}"                   2>&1 | tee    limit_${suffix1d}.log
    combine       $WKSPACEROOT $COMMON_ARGS -P $PAR2 -n _${suffix1d}                    2>&1 | tee -a limit_${suffix1d}.log
    echo "python build1DInterval.py ${VALS2} higgsCombine_${suffix1d}*.root ${PAR2} -b" 2>&1 | tee -a limit_${suffix1d}.log
    python build1DInterval.py ${VALS2} higgsCombine_${suffix1d}*.root ${PAR2} -b        2>&1 | tee -a limit_${suffix1d}.log
fi

if [ $# -ge 4 ]
then

########################################
# observed 2D limit:
    echo "combine $WKSPACEROOT $COMMON_ARGS $TWOD_ARGS -n _${suffix}_2Dobs" 2>&1 | tee    limit_${suffix}_2Dobs.log
    combine       $WKSPACEROOT $COMMON_ARGS $TWOD_ARGS -n _${suffix}_2Dobs  2>&1 | tee -a limit_${suffix}_2Dobs.log

# # expected 2D limit:
    echo "combine $WKSPACEROOT $COMMON_ARGS $TWOD_ARGS $EXPD_ARGS -n _${suffix}_2Dexp" 2>&1 | tee    limit_${suffix}_2Dexp.log
    combine       $WKSPACEROOT $COMMON_ARGS $TWOD_ARGS $EXPD_ARGS -n _${suffix}_2Dexp  2>&1 | tee -a limit_${suffix}_2Dexp.log

fi

if [ $DEBUGLEV -gt 1 ]
then
    gzip limit_${suffix}*.log
fi
