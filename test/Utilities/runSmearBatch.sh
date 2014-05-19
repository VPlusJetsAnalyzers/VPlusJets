#!/bin/bash
executable=$LOCALRT/src/ElectroWeakAnalysis/VPlusJets/test/Limits/batchit.sh
INDIR=/eos/uscms/store/user/lnujj/VBF_Higgs_v1/VBF_Higgs_5Aug_v1/

submit=submit.txt

if [ $# -eq 0 ]
then
    echo "Usage: $0 <eosdir>"
    echo "       set eosdir to subdir of /eos/uscms/store/user/${USER}/..."
    exit
fi

printheader()
{
cat >$submit <<EOF
Executable = $executable
Universe = vanilla
Requirements = Memory > 250 && FileSystemDomain=="fnal.gov" && Disk > 500000 && Arch=="X86_64"
Notification = ERROR
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
transfer_input_files = ${PWD}/runandcp2eos.sh,${PWD}/mkSmearedJetPtTree.C
EOF
}

printargs() {
cat >>$submit <<EOF
Error = ${PWD}/log/${3}.stderr
Output = ${PWD}/log/${3}.stdout
Arguments = ${LOCALRT} ./runandcp2eos.sh SmearedTrees ${2}/${3}.root root -n -l -b -q mkSmearedJetPtTree.C(\"$1\",\"${2}\",\"${3}.root\",\"${4}\",\"${5}\")
Queue
EOF
}

printheader

for  i in   RD_el_WW_CMSSW532.root RD_el_WWToAnything_CMSSW532.root RD_el_WZ_CMSSW532.root RD_el_WpJ_PT180_Madgraph_CMSSW532.root RD_el_TTJets_poheg_CMSSW532.root RD_el_Mtt700to1000_CMSSW532.root RD_el_Mtt1000toinf_CMSSW532.root RD_el_STopS_T_CMSSW532.root RD_el_STopT_T_CMSSW532.root RD_el_STopTW_T_CMSSW532.root RD_el_STopS_Tbar_CMSSW532.root RD_el_STopT_Tbar_CMSSW532.root RD_el_STopTW_Tbar_CMSSW532.root  RD_mu_WW_CMSSW532.root RD_mu_WWToAnything_CMSSW532.root RD_mu_WZ_CMSSW532.root RD_mu_WpJ_PT180_Madgraph_CMSSW532.root RD_mu_TTJets_poheg_CMSSW532.root RD_mu_Mtt700to1000_CMSSW532.root RD_mu_Mtt1000toinf_CMSSW532.root RD_mu_STopS_T_CMSSW532.root RD_mu_STopT_T_CMSSW532.root RD_mu_STopTW_T_CMSSW532.root RD_mu_STopS_Tbar_CMSSW532.root RD_mu_STopT_Tbar_CMSSW532.root RD_mu_STopTW_Tbar_CMSSW532.root
do
  echo Processing $i ...
  dname=`dirname $i`
  fname=`basename $i .root`
  printargs $INDIR . $fname GroomedJet_CA8_pt GroomedJet_CA8_eta
done
