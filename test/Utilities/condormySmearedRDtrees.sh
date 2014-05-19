#!/bin/bash
#
#
MYSRCDIR=${LOCALRT}/src/ElectroWeakAnalysis/VPlusJets/
executable=${MYSRCDIR}/test/Utilities/batchit.sh
MYRUN=mkSmearedJetPtTree
INDIR=/eos/uscms/store/user/lnujj/VBF_Higgs_v1/VBF_Higgs_5Aug_v1/

submit=submit.txt

printheader() {
cat >$submit <<EOF
Executable = $executable
Universe = vanilla
Requirements = Memory > 250 && FileSystemDomain=="fnal.gov" && Disk > 500000 && Arch=="X86_64"
Notification = ERROR
Should_Transfer_Files = YES
transfer_input_files =${MYSRCDIR}/test/Utilities/${MYRUN}.C
WhenToTransferOutput = ON_EXIT
EOF
}

printargs() {
cat >>$submit<<EOF
Error = ${MYRUN}_${BASE}.stderr
Output = ${MYRUN}_${BASE}.stdout
Arguments = ${MYSRCDIR} root -n -b -q -l ${MYRUN}.C(\"${INDIR}\",\".\",\"${FILE}\",\"GroomedJet_CA8_pt\",\"GroomedJet_CA8_eta\")
Queue
EOF
}

printheader

########################################
## make file list of your choice here:
########################################
for FILE in `ls $INDIR | grep "WZ_minPt150_amcnlo_add" | grep -v txt | grep -o "RD_.*root"`
do
  echo $FILE
  BASE=`basename  $FILE .root`
  printargs
done
