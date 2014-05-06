#!/bin/bash
#
#DIR=/uscms_data/d2/andersj/Wjj/2012/CMSSW_5_3_3_patch2/src/ElectroWeakAnalysis/VPlusJets/test/HWWFit24July2013
#DIR=/uscms_data/d2/andersj/Wjj/2013/CMSSW_6_2_0/src/ElectroWeakAnalysis/VPlusJets/test/HWWTest_HCP
#DIR=/uscms_data/d2/andersj/Wjj/2013/CMSSW_6_2_0/src/ElectroWeakAnalysis/VPlusJets/test/HWWTweakFit
#DIR=/uscms_data/d2/andersj/Wjj/2013/CMSSW_6_2_0/src/ElectroWeakAnalysis/VPlusJets/test/HWWHCP_NewFit
#DIR=/uscms_data/d2/andersj/Wjj/2013/CMSSW_6_1_2/src/ElectroWeakAnalysis/VPlusJets/test/HWWHCP_scan
#DIR=/uscms_data/d2/andersj/Wjj/2013/CMSSW_6_1_2/src/ElectroWeakAnalysis/VPlusJets/test/HWWFit4Oct

CMSSW_REL=CMSSW_6_1_2
#JAKEDIR=HWWFitAllOneCut
#JAKEDIR=MVANominal30Oct
#JAKEDIR=MVANominal18Nov2013 
#JAKEDIR=FloatFit18Nov2013 
#JAKEDIR=FloatFit300Wider18Nov2013 
#JAKEDIR=FloatFit300Wider20Nov2013 
#JAKEDIR=toyFit1.0Sig_newGen
#JAKEDIR=NomFits19Feb2014
#JAKEDIR=AltFits19Feb2014
#JAKEDIR=Fits4Feb2014
#JAKEDIR=HWWPolyFits20Mar2014_mva0.0 
#JAKEDIR=HWWFits07Apr2014_mva0.0
#JAKEDIR=HWWPolyFits07Apr2014_mva0.0 
#JAKEDIR=HWWPolyFits11Apr2014_mva0.0
JAKEDIR=HighMassFittingFiles

#mvacut=0.001
#mvacut=0.4
#mvacut=0.5
#mvacut=0.6
#mvacut=0.7
#mvacut=0.8
#mvacut=0.9

#JAKEDIR=MVA${mvacut}

#cfgtag=jake${JAKEDIR}

DIR=hwwinputdir

ln -s -f /uscms_data/d2/andersj/Wjj/2013/${CMSSW_REL}/src/ElectroWeakAnalysis/VPlusJets/test/${JAKEDIR} $DIR


for i in 1 2 3 5 7 10
do
  cprime=`echo "scale=1; $i/10" | bc; exit`
  for (( j=0; j<=5; j++ ))
  do
    brnew=`echo "scale=1; $j/10" | bc; exit`
    cfgtag=`printf "Cprime_%3.1f_BRnew_%3.1f" ${cprime} ${brnew}`
#   simulated points:
#    for mass in 170 180 190 200 250 300 350 400 450 500 550 600

#   interpolated points:
    for mass in 175 185 195 205 210 215 220 225 230 235 240 245 255 260 265 270 275 280 285 290 295 300 310 320 330 340 360 370 380 390 400 420 440 460 480 520 540 560 580
    do
      for flavor in electron muon
      do
	rootfile="${DIR}/HWW${mass}lnujj_${flavor}_1D2Fit_output_${cfgtag}.root"
	if [ -f $rootfile ]
	then
	    ./makeDataCardFromRooWorkspace.exe $rootfile $cfgtag 2>&1 | tee makedc_${mass}lnujj_${flavor}_${cfgtag}.log
	else
	    echo "Couldn't find $rootfile, skipping"
	fi
      done
# combine the flavors:
      ./movem 's#M=#M-#g' datacard_8TeV_hww??nu2j_${cfgtag}-M=${mass}.txt
      hwwelnu2j=datacard_8TeV_hwwelnu2j_${cfgtag}-M-${mass}.txt
      hwwmunu2j=datacard_8TeV_hwwmunu2j_${cfgtag}-M-${mass}.txt
      if [ -f $hwwelnujj -a -f $hwwmunujj ]
      then
	  echo "combineCards.py hwwmunu2j=${hwwmunu2j} hwwelnu2j=${hwwelnu2j} >datacard_8TeV_hww2chan_${cfgtag}-M=${mass}.txt"
	  combineCards.py hwwmunu2j=${hwwmunu2j} hwwelnu2j=${hwwelnu2j} >datacard_8TeV_hww2chan_${cfgtag}-M=${mass}.txt
	  ./movem 's#M-#M=#g' datacard_8TeV_hww??nu2j_${cfgtag}-M-${mass}.txt
      fi
    done
  done
done

