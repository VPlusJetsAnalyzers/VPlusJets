#!/bin/bash
BASEDIR=/uscms/home/pdudero/data/HWW/HWWmvaScan131024
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
transfer_input_files = $JAKESRCDIR/HWWFitting.tar.gz,$MYSRCDIR/runScanBatch.sh
WhenToTransferOutput = ON_EXIT
EOF
}

printargs() {
    echo "Error = fitmu_M=${masspt}_MVA${cutval}.stderr" >>$submit
    echo "Output = fitmu_M=${masspt}_MVA${cutval}.stdout" >>$submit
    echo "Arguments = $JAKESRCDIR ./runScanBatch.sh ./runHiggsNew.py --doShapes --doFit --expLimit 20 --mva $cutval $masspt" >>$submit
    echo "Queue" >>$submit
    echo "Error = fitel_M=${masspt}_MVA${cutval}.stderr" >>$submit
    echo "Output = fitel_M=${masspt}_MVA${cutval}.stdout" >>$submit
    echo "Arguments = $JAKESRCDIR ./runScanBatch.sh ./runHiggsNew.py --doShapes --doFit --electrons --expLimit 20 --mva $cutval $masspt" >>$submit
    echo "Queue" >>$submit
}

for cutval in 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
  cutdir=MVA${cutval}
  submit=submit_cut=${cutval}.txt
  mkdir -p $cutdir
  cd $cutdir
  printheader
  for masspt in 170 180 190 200 250 300 350 400 450 500 550 600
  do
    printargs
  done
  condor_submit $submit
  cd -
done


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
