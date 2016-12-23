#! /bin/sh
echo "Starting" 
pwd
##site_arch_choice=slc4_ia32_gcc345
##export BUILD_ARCH=`echo $site_arch_choice | awk -F_ '{print $1"_"$2}'`
##export SCRAM_ARCH=$site_arch_choice
export BUILD_ARCH=slc5_amd64_gcc462
export SCRAM_ARCH=slc5_amd64_gcc462
##setenv X509_USER_PROXY /uscms/home/aperloff/x509up_u43841
tar -xvf FullRelease5_3_15.tar.gz
##source /cvmfs/cms.cern.ch/cmsset_default.sh
##scram p CMSSW_5_3_2_patch4
mv *.* CMSSW_5_3_15_patch1/src/ElectroWeakAnalysis/VPlusJets/test/
cd CMSSW_5_3_15_patch1/src/ElectroWeakAnalysis/VPlusJets/test/
scramv1 b ProjectRename
scram b
pwd
##cd /uscms/home/ilyao/EWKVPLUSJETS/CMSSW_5_3_2_patch4/src/ElectroWeakAnalysis/VPlusJets/test
##source /uscmst1/prod/grid/gLite_SL5.csh
##source /uscmst1/prod/sw/cms/shrc prod
##pwd
##eval `scram runtime -csh`
eval `scram runtime -sh`
##cmsenv
##cd -
which root
##tar -xvf AllToyFitDiboson8TeV_GeneralFiles_NoLib.tar.gz
##source root/bin/thisroot.csh
##source root/bin/thisroot.sh
##which root
##ls

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
echo "muons, Btag Resolved. Fit the default fit result with the default configuration"
initRand="$(( 2*${Process} + 3066302 ))"
echo "initRand=${initRand}"
CstrNoBtag="runToyMCFits.cc+(${StartToy},${EndToy},3,${initRand},\"Fit_DibosonBtag_mu_log.txt\",\"python -i runDiboson8TevFit.py -j 2 -m Diboson8TeVConfig --runPdfGenToySim --genParamFile DibosonBtagFitOutPars_mu.txt --btag --noNull --nosig\",\"topDibosonBtagParameters_mu.txt dibosonDibosonBtagParameters_mu.txt WpJDibosonBtagParameters_mu.txt WHbbDibosonBtagParameters_mu.txt \",\"./\",\"BtagFitMuon_ShapeBased_VariableYield_GenWpJDataDefaultFitDefault_\",true,true)"
echo "CstrNoBtag=$CstrNoBtag"
eval "root -l -b -q '$CstrNoBtag'"
mv *.log ../../../../../
pwd