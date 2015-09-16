#!/bin/bash
#
# start with el_boosted.root, mu_boosted.root, and ATGC_shape_coefficients.root in
#  $CMSSW_BASE/src/ElectroWeakAnalysis/VPlusJets/test/TGC
#

if [ $# -lt 1 ]
then
    echo "Usage: $0 configtag"
    exit
fi

#TUCHAN=WVsemilep2chan
MUCHAN=muboosted
ELCHAN=elboosted
TUCHAN=elmu
# MUCHAN=WV_mu
# ELCHAN=WV_el
CFGTAG=$1
HZ2DM=CombinedEWKAnalysis.CommonTools.HagiwaraAndZeppenfeldTwoDimensionalModel
#DC=datacard_8TeV_ATGC
DC=WVsemileptonic

# buildWVworkspace.py found in 
# https://github.com/lgray/CombinedEWKAnalysis/blob/master/CommonTools/test/buildWVworkspace.py
#
buildws()
{
    plane=$1
    python buildWVworkspace.py muon $plane ${CFGTAG} -b
    python buildWVworkspace.py electron $plane ${CFGTAG} -b
}

# this is now obsolete
#    ./makeDataCardShapes4tgc.exe WVsemileptonic_muboosted_${plane}_ws.root ${plane}_${CFGTAG}
#    ./makeDataCardShapes4tgc.exe WVsemileptonic_elboosted_${plane}_ws.root ${plane}_${CFGTAG}

makecards()
{
    plane=$1
    combineCards.py ${DC}_${ELCHAN}=${DC}_${ELCHAN}_${plane}_${CFGTAG}.txt ${DC}_${MUCHAN}=${DC}_${MUCHAN}_${plane}_${CFGTAG}.txt  >${DC}_${TUCHAN}_${plane}_${CFGTAG}.txt 
}

makefinalws()
{
    plane=$1
    dc="${DC}_${TUCHAN}_${plane}_${CFGTAG}.txt"
    text2workspace.py -m 126 $dc -o wkspace_ATGC_${TUCHAN}_${plane}_${CFGTAG}.root -P ${HZ2DM}:${plane}Model
}

makefinalEFT1Dws()
{
    plane=$1
    range=$2
    dc="${DC}_${TUCHAN}_${plane}_${CFGTAG}.txt"
    text2workspace.py -m 126 $dc -o wkspace_ATGC_${TUCHAN}_${plane}_EFT1D_${CFGTAG}.root -P ${HZ2DM}:${plane}Model --PO ${range}
}

dolimit()
{
    plane=$1
    if [ "$4" == "" ]
    then
	var1=$2
	./combine_deltanll_atgc.sh  wkspace_ATGC_${TUCHAN}_${plane}_EFT1D_${CFGTAG}.root $var1
    else
	var1=$2
	var2=$3
	var3=$4
	./combine_deltanll_atgc.sh  wkspace_ATGC_${TUCHAN}_${plane}_${CFGTAG}.root $var1 $var2 $var3
    fi
}

# I do it in this order because I like to make all the ws' at the same time,
# and all the cards at the same time, etc. It looks nicer in a time-ordered dir listing

buildws dkglZ
buildws dg1lZ
buildws dkgdg1
buildws dkzdg1zToCW
buildws dkzdg1zToCB
buildws dg1zlZToCWWW

makecards dkglZ
makecards dg1lZ
makecards dkgdg1
makecards dkzdg1zToCW
makecards dkzdg1zToCB
makecards dg1zlZToCWWW

makefinalws dkglZ
makefinalws dg1lZ
makefinalws dkgdg1
makefinalEFT1Dws dkzdg1zToCW range_dg1=-7,7
makefinalEFT1Dws dkzdg1zToCB range_dkg=-30,30
makefinalEFT1Dws dg1zlZToCWWW range_lZ=-5,5

dolimit dkglZ lZ=-0.025,0.025 dkg=-0.15,0.15 dg1=0.0
dolimit dg1lZ lZ=-0.025,0.025 dg1=-0.06,0.06 dkg=0.0
#dolimit dkgdg1 dkg=-0.15,0.15 dg1=-0.06,0.06 lZ=0.0
dolimit dkgdg1 dg1=-0.06,0.06 dkg=-0.15,0.15 lZ=0.0
dolimit dkzdg1zToCW dg1=-7,7
dolimit dkzdg1zToCB dkg=-30,30
dolimit dg1zlZToCWWW lZ=-5,5

# # for slices:
# dolimit dkglZ lZ=-0.01,0.01 dkg=-0.05,0.05 dg1=0.0
# dolimit dg1lZ lZ=-0.01,0.01 dg1=-0.01,0.01 dkg=0.0
# #dolimit dkgdg1 dkg=-0.05,0.05 dg1=-0.01,0.01 lZ=0.0
# dolimit dkgdg1 dg1=-0.01,0.01 dkg=-0.05,0.05 lZ=0.0

#root
#.L atgcplotLimit.C+
#atgcplotLimit("higgsCombine_ATGC*2D*par1*par2*.root")
