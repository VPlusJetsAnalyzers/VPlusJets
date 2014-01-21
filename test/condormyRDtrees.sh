#!/bin/bash
#
# This script generates a text file that can be read by a
# "condor_submit" command, submitting jobs to the condor queue for
# processing. The jobs are dedicated to 'reducing' lnujj ntuples using
# the kana*.C scripts. The beauty of this setup is that it does not
# require you to have a working area yourself; you can copy someone
# else's working area entirely (assuming you have read permissions) and
# run these scripts in batch to create a set of reduced trees that end
# up in whatever destination you choose, or mix and match; take the
# code from one area and the tables from another.Caveat emptor; this
# script cannot compensate for dumb choices.
# 
# Instructions for use:
#
# 1. Create/Locate a ElectroWeakAnalysis/VplusJets/test working area. The regular
#    MyRunElec.C and MyRunMuon.C scripts should have been used to compile the
#    complete set of code necessary to run the associated kana*.C scripts.
#
# 2. Edit the environment variables below to tell the script where to find things.
#    The "transfer_input_files" field is a comma-separated list that tells condor
#    what files to transfer to the worker node before starting processing.
#
# 3. Edit the list of tree codes in the loop at bottom for the trees you want
#    processed. The codes are interpreted by kana*.C.
#
# 4. Run this script, generate 'submit.txt'
#
# 5. 'cp submit.txt' to a directory that has the disk/quota space to receive all
#    the output once it's done. I suggest 3DayLifetime scratch space. Condor will
#    deliver your output to this directory. You can copy from there to eos or dcache
#    once you've verified the output.
# 
# 6. 'cd' to the directory with the submit.txt file.
#
# 7. Let 'er rip: "condor_submit submit.txt"
#
MYSRCDIR=${LOCALRT}/src/ElectroWeakAnalysis/VPlusJets/
IFTABLEDIR=${MYSRCDIR}/test/InterferenceTable2012
EFFTABLEDIR=${MYSRCDIR}/test/EffTable2012
CLASSDIR=${MYSRCDIR}/test/ClassifierOut
executable=${MYSRCDIR}/test/Limits/batchit.sh

if [ $# -ne 1 ]
then
    echo "Usage: $0 [el|mu]"
    exit
fi

FLAV=$1
#FLAV=mu
#FLAV=el

submit=submit${FLAV}.txt
#
# Construct long list of TMVA libraries to load
#
CLASSLIST=""
for nJ in nJ2 nJ3 VBF; do
  for MASS in 170 180 190 200 250 300 350 400 450 500 550 600; do
    CLASSLIST=$CLASSLIST,${CLASSDIR}/TMVAClassification_${MASS}_${nJ}_${FLAV}_Likelihood.class_C.so
  done
done

CLASSLIST=$CLASSLIST,${CLASSDIR}/TMVAClassification_126_VBF_${FLAV}_Likelihood.class_C.so

for nJ in nJ2 nJ3; do
  for MASS in 400 450 500 550 600 700 800 900 1000; do
    for TYPE in interferencedown interferencenominal interferenceup; do
      CLASSLIST=$CLASSLIST,${CLASSDIR}/TMVAClassification_${MASS}_${nJ}_${FLAV}_${TYPE}_Likelihood.class_C.so
    done
  done
done
for nJ in nJ2 nJ3; do
  for QG in noqg withqg; do
    CLASSLIST=$CLASSLIST,${CLASSDIR}/TMVAClassification_${QG}_${nJ}_${FLAV}_BDT.class_C.so
  done
done

TABLES=""
for MASS in 400 450 500 550 600 700 800 900 1000; do
    TABLES=$IFTABLEDIR/ratio${MASS}.txt
done

if [ "$FLAV" == "mu" ]
then
    KANA=kanamuon
    MYRUN=MyRunMuonBatch
    TABLES=${TABLES},$EFFTABLEDIR/scaleFactor-Run2012ABCD-RecoToIso.txt,$EFFTABLEDIR/efficiency-Run2012ABCD-IsoToIsoMuHLT.txt
elif [ "$FLAV" == "el" ]
then
    KANA=kanaelec
    MYRUN=MyRunElecBatch
    TABLES=${TABLES},$EFFTABLEDIR/scaleFactor-Run2012ABCD-GsfElectronToId.txt,$EFFTABLEDIR/scaleFactor-Run2012ABCD-SCToElectron.txt,$EFFTABLEDIR/efficiency-Run2012ABCD-WP80ToHLTEle.txt,$EFFTABLEDIR/FullyEfficient.txt,$EFFTABLEDIR/FullyEfficient_Jet2NoJet1.txt,$EFFTABLEDIR/FullyEfficient_MHT.txt
else
    echo "Unknown lepton flavor $FLAV"
    echo "Usage: $0 [el|mu]"
    exit
fi

TABLES=${TABLES},$MYSRCDIR/test/Data190456-208686_PileupHistogram.root,$MYSRCDIR/test/sampledb2012.txt

COMMONSRC=${MYSRCDIR}/test/${MYRUN}.C,${MYSRCDIR}/test/${KANA}_C.so,${MYSRCDIR}/test/Resolution_cc.so,${MYSRCDIR}/src/METzCalculator_cc.so,${MYSRCDIR}/src/QGLikelihoodCalculator_C.so,${MYSRCDIR}/test/EffTableReader_cc.so,${MYSRCDIR}/test/EffTableLoader_cc.so

printheader() {
cat >$submit <<EOF
Executable = $executable
Universe = vanilla
Requirements = Memory > 250 && FileSystemDomain=="fnal.gov" && Disk > 500000 && Arch=="X86_64"
Notification = ERROR
Should_Transfer_Files = YES
transfer_input_files =${COMMONSRC},${CLASSLIST},${TABLES}
WhenToTransferOutput = ON_EXIT
EOF
}

printargs() {
cat >>$submit<<EOF
Error = ${KANA}_2012${treecode}.stderr
Output = ${KANA}_2012${treecode}.stdout
Arguments = ${MYSRCDIR} root -n -b -q -l ${MYRUN}.C(2012${treecode},0,0)
Queue
EOF
}

printheader
#for treecode in 0000  
#for treecode in 3126 3127 10241 1010
#for treecode in 0000 1002 1003 1004 1005 1006 1007 1010 10201 10202 1021 1022 1023 1024 10241 1030 1031 1032 1040 3126 3127
for treecode in 0001 0002
do
  printargs
done
