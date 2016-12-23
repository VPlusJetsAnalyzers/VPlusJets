#! /bin/sh
echo "Starting" 
pwd
cd /uscms/home/ilyao/EWKVPLUSJETS/CMSSW_5_3_2_patch4/src/ElectroWeakAnalysis/VPlusJets/test
##source /uscmst1/prod/grid/gLite_SL5.csh
source /uscmst1/prod/sw/cms/shrc prod
pwd
##eval `scram runtime -csh`
eval `scram runtime -sh`
##cmsenv
cd -
which root
tar -xvf AllToyFitDiboson8TeV.tar.gz
##source root/bin/thisroot.csh
##source root/bin/thisroot.sh
##which root

#### Define the process and the number of toy datasets to fit
Process="${BASH_ARGV[1]}"
NToys="${BASH_ARGV[0]}"
echo "Process=${Process}"
echo "NToys=${NToys}"
StartToy="$(( ${Process} * $NToys ))"
EndToy="$(( ${NToys} + ${StartToy} - 1))"
echo "StartToy=NToys*Process-1=$StartToy"
echo "EndToy=NToys+StartToy-1=$EndToy"
#### Run the Toy Fits
echo "muons, Boosted:"
initRand="$(( ${Process} * 2 )+33807)"
echo "initRand=$initRand"
CstrNoBtag="runToyMCFits8TeV_Diboson.cc+(${StartToy},${EndToy},false,\"FitBoostedValInput_mu.log\",${initRand})"
echo "CstrNoBtag=$CstrNoBtag"
eval "root -l -b -q '$CstrNoBtag'"
##echo "3jets:"
##Cstr3j="runToyMCFits.cc+(3,${StartToy},${EndToy})"
##echo "Cstr3j=$Cstr3j"
##eval "root -l -b -q '$Cstr3j '"
##echo "Finished"
