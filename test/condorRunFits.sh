#!/bin/bash
JAKESRCDIR=$(pwd)
executable=${JAKESRCDIR}/batchit.sh

printheader() {
cat >$submit <<EOF
Executable = $executable
Universe = vanilla
Requirements = Memory > 250 && FileSystemDomain=="fnal.gov" && Disk > 500000 && Arch=="X86_64"
Notification = ERROR
Should_Transfer_Files = YES
transfer_input_files = $JAKESRCDIR/HWWFitting.tar.gz,$JAKESRCDIR/runScanBatch.sh
WhenToTransferOutput = ON_EXIT
EOF
}

printargs() {
    echo "Error = fitmu_HWW${masspt}_${cutval}.stderr" >>$submit
    echo "Output = fitmu_HWW${masspt}_${cutval}.stdout" >>$submit
    echo "Arguments = $JAKESRCDIR ./runScanBatch.sh ./runHiggsNew.py --doShapes --doFit --expLimit 20 --obsLimit --xrootd $masspt ${modeArg}${mvaArg}${unbinned}${fixZeros}" >>$submit
    echo "Queue" >>$submit
    echo "Error = fitel_HWW${masspt}_${cutval}.stderr" >>$submit
    echo "Output = fitel_HWW${masspt}_${cutval}.stdout" >>$submit
    echo "Arguments = $JAKESRCDIR ./runScanBatch.sh ./runHiggsNew.py --doShapes --doFit --electrons --expLimit 20 --obsLimit --xrootd $masspt ${modeArg}${mvaArg}${unbinned}${fixZeros}" >>$submit
    echo "Queue" >>$submit
}

mode=""
mvaVal=0.0
unbinned=""
binnedLabel=""

while getopts "M:F:U" opt ; do
    case $opt in
	M)
	    mvaVal=$OPTARG
	    ;;
	F)
	    mode=$OPTARG
	    ;;
	U)
	    unbinned=" --unbinned "
	    binnedLabel="unbinned"
	    ;;
	?)
	    echo "invalid option"
	    exit
	    ;;
    esac
done

shift $((OPTIND -1))

if [ "$#" -gt 0 ] ; then
    masses=$@
else
    masses="170 180 190 200 250 300 350 400 450 500 550 600"
fi
echo "masses: $masses"

python pyroot_logon.py

echo -n "tarring ... "
tar czf HWWFitting.tar.gz -T filesForHWWFit.txt -X excludeForHWWFit.txt
echo "done"

modeArg="--mode HWW1D2FitsConfig${mode}"
mvaArg=""
DateString=$(date "+%d%b%Y")
if [ $(echo "$mvaVal > 0.1" | bc) -eq 1 ] ; then
    mvaArg=" --mva ${mvaVal}"
fi
cutval=HWW${mode}Fits${DateString}${binnedLabel}_mva${mvaVal}
cutdir=${cutval}
submit=submit_${cutval}.txt
mkdir -p $cutdir
cd $cutdir
printheader
for masspt in $masses
do
  fixZeros=""
  # if [ $masspt -gt 401 ] ; then
  #   if [ "$mode" == "Poly" ] ; then
  #     fixZeros=" --fixToZero"
  #   fi
  # fi
  printargs
done
condor_submit $submit
cd -
