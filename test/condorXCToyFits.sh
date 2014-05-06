#!/bin/bash
# BASEDIR=/uscms/home/pdudero/data/HWW/HWWmvaScan131024
# JAKESRCDIR=/uscms_data/d2/andersj/Wjj/2013/CMSSW_6_1_2/src/ElectroWeakAnalysis/VPlusJets/test/
JAKESRCDIR=$(pwd)
executable=${JAKESRCDIR}/batchit.sh

echo $JAKESRCDIR

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
    echo "Error = fitToys_M${masspt}_Sig${cutval}_${start}_${end}.stderr" >>$submit
    echo "Output = fitToys_M${masspt}_Sig${cutval}_${start}_${end}.stdout" >>$submit
    echo "Arguments = $JAKESRCDIR ./runScanBatch.sh ./runToyBatch.py --mH ${masspt} -s ${start} -e ${end} --injectS ${cutval} --genConfig HWW1D2FitsConfig${Alternative}${highlow} --xc 1D2FitsParameters/WpJHWW${masspt}Parameters_mWW_${alternative}.txt --mode HWW1D2FitsConfig${BaseLine}${highlow} --WpJ 1D2FitsParameters/WpJHWW${masspt}Parameters_mWW_${baseline}.txt ${unbinned}" >>$submit
    echo "Queue" >>$submit
}

cutval=1.0
Alternative="Alt"
BaseLine="Nominal"
highlow=""
unbinned=""
binnedLabel="binned"

while getopts "US:G:F:H:" opt ; do
    case $opt in
	S)
	    cutval=$OPTARG
	    ;;
	G)
	    Alternative=$OPTARG
	    ;;
	F)
	    BaseLine=$OPTARG
	    ;;
	H)
	    highlow=$OPTARG
	    ;;
	U)
	    unbinned="--unbinned "
	    binnedLabel="unbinned"
	    ;;
	?)
	    echo "invalid option"
	    exit
	    ;;
    esac
done

shift $((OPTIND -1))
# echo "$#"
# echo "$*"
if [ "$#" -gt 0 ] ; then
    masses=$@
else
    masses="170 180 190 200 250 300 350 400 450 500 550 600"
fi
echo "masses: $masses"

#rm HWWFitting.tar.gz	   
python pyroot_logon.py
echo -n "tarring ... "
tar czf HWWFitting.tar.gz -T filesForHWWFit.txt -X excludeForHWWFit.txt
echo "done"

alternative="$(echo $Alternative | tr '[:upper:]' '[:lower:]')"
baseline="$(echo $BaseLine | tr '[:upper:]' '[:lower:]')"
DateString=$(date "+%d%b%Y")
cutdir=toyFit${DateString}_XC_${cutval}Sig_${binnedLabel}_${Alternative}${BaseLine}${highlow}
if [ "$Alternative" = "Nominal" ]; then
    Alternative=""
fi
if [ "$BaseLine" = "Nominal" ]; then
    BaseLine=""
fi
echo "injecting ${cutval}x signal"
echo "generating: $alternative fitting: $baseline highlow: $highlow"
submit=submit_${cutdir}.txt
mkdir -p $cutdir
cd $cutdir
printheader
# for masspt in 170 180 190 200 250 300 350 400 450 500 550 600
for masspt in $masses
do
    echo "mass: $masspt"
    for (( start=0; start<=450; start+=50 ))
    do
	let "end = $start+50"
        # echo "mH: $masspt start: $start end: $end"
	printargs
    done
done
condor_submit $submit
cd -
