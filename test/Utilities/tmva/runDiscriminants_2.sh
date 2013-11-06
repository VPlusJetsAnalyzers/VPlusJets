#!/bin/bash

source cmslpc_standalone_setup.sh

elements=1
startMass=( 250)
#elements=11
#startMass=( 170 180 190 200 250 300 350 400 450 500 600)

outdir="WWlvjj_TMVAoutput_final_w"

if [ -e $outdir ]; then
echo $outdir already exists
else
mkdir $outdir
fi

echo ${elements}
a=0
while ((a < ${elements}))
do

## ------- active elements --------
mStart=${startMass[a]}
dummy='""'
chanmu='"mu"'
chanel='"el"'
echo $mStart ", " $chanmu ", " $chanel
# mu, 2jets
root -l -b -q "WWTMVAClassification.C(${dummy},${mStart},2,${chanmu})"   
# mu, 3jets
#root -l -b -q "WWTMVAClassification.C(${dummy},${mStart},3,${chanmu})"   
## el, 2jets
root -l -b -q "WWTMVAClassification.C(${dummy},${mStart},2,${chanel})"   
## el, 3jets
#root -l -b -q "WWTMVAClassification.C(${dummy},${mStart},3,${chanel})"   
#
### make all the plots
fname1='"TMVA_'${mStart}'_nJ2_mu.root"'
#fname2='"TMVA_'${mStart}'_nJ3_mu.root"'
fname3='"TMVA_'${mStart}'_nJ2_el.root"'
#fname4='"TMVA_'${mStart}'_nJ3_el.root"'
### plot variables
root -l -b -q "variables.C(${fname1})"
#root -l -b -q "variables.C(${fname2})"
root -l -b -q "variables.C(${fname3})"
#root -l -b -q "variables.C(${fname4})"
### plot correlations
root -l -b -q "correlations.C(${fname1})"
#root -l -b -q "correlations.C(${fname2})"
root -l -b -q "correlations.C(${fname3})"
#root -l -b -q "correlations.C(${fname4})"
### plot efficiencies
root -l -b -q "efficiencies.C(${fname1})"
#root -l -b -q "efficiencies.C(${fname2})"
root -l -b -q "efficiencies.C(${fname3})"
#root -l -b -q "efficiencies.C(${fname4})"
### plot mvas
root -l -b -q "mvas.C(${fname1})"
#root -l -b -q "mvas.C(${fname2})"
root -l -b -q "mvas.C(${fname3})"
#root -l -b -q "mvas.C(${fname4})"


#root -l -b -q "mvaeffs.C(${fname1})"
#root -l -b -q "mvas.C(${fname2})"
#root -l -b -q "mvaeffs.C(${fname3})"
#root -l -b -q "mvas.C(${fname4})"


#
## ------- active elements --------

if [ -e $outdir/m${mStart} ]; then
echo $outdir/m${mStart} already exists
else
mkdir $outdir/m${mStart}
fi

mv TMVA_${mStart}_* $outdir/m${mStart}/.
#mv plots/TMVA_${mStart}_* $outdir/m${mStart}/.
mv weights/TMVAClassification_${mStart}_* $outdir/m${mStart}/.

let a=$a+1
done

