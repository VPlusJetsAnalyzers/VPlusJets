#!/bin/bash
# BASEDIR=/uscms/home/pdudero/data/HWW/HWWmvaScan131024
JAKESRCDIR=/uscms_data/d2/andersj/Wjj/2013/CMSSW_6_1_2/src/ElectroWeakAnalysis/VPlusJets/test/
MYSRCDIR=/uscms/home/pdudero/work/CMSSW_6_1_2/src/ElectroWeakAnalysis/VPlusJets/test/Limits/
executable=${MYSRCDIR}/batchit.sh

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
    echo "Arguments = $JAKESRCDIR ./runScanBatch.sh ./runToyBatch.py --mH ${masspt} -s ${start} -e ${end} --injectS ${cutval} --genConfig HWW1D2FitsConfig" >>$submit
    echo "Queue" >>$submit
}

cutval=0.0
cutdir=toyFit${cutval}Sig_allGen
submit=submit_${cutdir}.txt
mkdir -p $cutdir
cd $cutdir
printheader
# for masspt in 170 180 190 200 250 300 350 400 450 500 550 600
for masspt in 180 200 300 400 600
do
  for (( start=0; start<=450; start+=50 ))
  do
    let "end = $start+50"
    # echo "mH: $masspt start: $start end: $end"
    printargs
  done
done
condor_submit $submit
cd -


# do
#   outputfile="explimvscutval-M=${masspt}_${channel}.tab"
#   echo "#Cutval	explim" >$outputfile
#   for cutvalue in `ls -d $BASEDIR/cutvalue=0.* | egrep -o "cutvalue=[0-9.]+" | egrep -o "[0-9.]+"`
#     do
#     LOGFILE=$BASEDIR/cutvalue=${cutvalue}/limit*${channel}*${masspt}*log
#     explim=""
#     if stat -t $LOGFILE >/dev/null 2>&1
#     then
# 	explim=`grep -h "Expected 50" $LOGFILE | egrep -o "r < [0-9.]+" | egrep -o "[0-9.]+"`
#     fi
#     if [ "$explim" = "" ]
#     then
# 	explim=-1
#     fi
#     echo "$cutvalue	$explim" >>$outputfile
#   done
# done
