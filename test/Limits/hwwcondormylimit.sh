#!/bin/bash
#BASEDIR=/uscms/home/pdudero/data/HWW/HWWmvaOptFit130219
SRCDIR=/uscms/home/pdudero/work/CMSSW_6_1_2/src/ElectroWeakAnalysis/VPlusJets/test/Limits/

#executable=${SRCDIR}/hwwrunCombine.sh
executable=${SRCDIR}/batchit.sh
cmbscript1=${SRCDIR}/combine_asymp.sh

if [ $# -lt 1 ]
then
    echo "Usage: $0 <fileglob> ( e.g., datacard*.txt)"
    exit
fi

if [ "$LOCALRT" == "" ]
then
    echo "ERROR: cms environment LOCALRT not set"
    exit
fi

submit=submit.txt

printheader() {
cat >$submit <<EOF
Executable = $executable
Universe = vanilla
Requirements = Memory > 250 && FileSystemDomain=="fnal.gov" && Disk > 500000 && Arch=="X86_64"
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Notification = ERROR
EOF
}

#transfer_input_files = $SRCDIR/ggHtable8tev.txt,$SRCDIR/twikiBRtable.txt,$SRCDIR/vbfHtable8tev.txt,$SRCDIR/jetbinerrtable.txt

queuecards() {
    echo "transfer_input_files = $xferlist" >>$submit
    echo "Error = hwwlimit_${cfgtag}.log" >>$submit
    echo "Output = hwwlimit_${cfgtag}.log" >>$submit
    echo "Arguments = $arglist" >>$submit
    echo "Queue" >>$submit
}

printheader

for file in $*
do
  ROOTFILES=`grep -o "hwwinputdir.*root" $file | sort -u |xargs|sed 's# #,#g'`
  ADDTLFILES=",$cmbscript1,$ROOTFILES"
  xferlist="$file${ADDTLFILES}"
  scriptbase=`basename $cmbscript1`
  filebase=`basename $file`
  arglist="$LOCALRT ln -s . hwwinputdir; ./$scriptbase $filebase"
  base=`basename $file ${FILESUFFIX}`
  cfgtag=${base##datacard_}
  queuecards
done

# for rootf in $BASEDIR/*.root
#   do
#   mH=`echo $rootf | egrep -o "HWW[0-9]{3}" | egrep -o "[0-9]{3}"`
#   printargs
# done

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
