#import sys
#sys.path.append('/uscms/home/andersj/pyroot')
#del sys
import math
#import root_logon
#import pyroot_fwlite.py
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-H', '--mH', dest='mH', type='float',
                  help='Higgs Mass Point')
(opts, args) = parser.parse_args()

import pyroot_logon
from ROOT import TFile, TH1F, TCanvas, TLorentzVector, TTree, gSystem, gROOT
from DataFormats.FWLite import Events, Handle
from array import array

import HWWSignalShapes

#files = ['WHZH_GENSIM.root']
#from WHZHFiles import files

#files = ['dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/andersj/WHZH_Wjj_4_2_3_SIM/WHZH_Wjj_4_2_3_SIM/c3d0044b4728087b67b7d100b833fffe/WHZH_GENSIM_2_1_EUL.root']

# files = ['dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lnujj/PatTuples_8TeV_53X-v1/shuai/GluGluToHToWWToLAndTauNuQQ_M-350_8TeV-powheg-pythia6/SQWaT_PAT_53X_ggH350_central/829f288d768dd564418efaaf3a8ab9aa/pat_53x_v03_60_1_cPR.root']

files = [ x.rstrip() for x in open(args[0]) ]

print files[0]

events = Events (files)
handle = Handle('LHEEventProduct')
label = ('source')

def resetVars():

    Higgs_mass[0] = -999.
    NonRunning_wgt[0] = -999.
    Running_wgt[0] = -999.
    Higgs_i = -1

outputFile = TFile('output.root', 'recreate')
outTree = TTree("outTree", "outTree")

Higgs_mass = array('d', [0.])
NonRunning_wgt = array('d', [0.])
Running_wgt = array('d', [0.])

outTree.Branch('Higgs_mass', Higgs_mass, 'Higgs_mass/D')
outTree.Branch('NonRunning_wgt', NonRunning_wgt, 'NonRunning_wgt/D')
outTree.Branch('Running_wgt', Running_wgt, 'Running_wgt/D')

decayCnts = {}
sumWgts = 0.
sumNonRunning = 0.
Higgs_i = -1

for (eventN,event) in enumerate(events):
    resetVars()
    event.getByLabel(label, handle)
    source = handle.product()
    genParts = source.hepeup()

    if (eventN % 100 == 0):
        print "Record:",eventN,
        print "Run:",event.object().id().run(),
        print "event:",event.object().id().event()
    for parti in range(0,genParts.NUP):
        if (abs(genParts.IDUP.at(parti)) == 25):
            Higgs_i = parti
            Higgs_mass[0] = genParts.PUP.at(parti).x[4]
            Running_wgt[0] = HWWSignalShapes.HiggsCPWeight(opts.mH,
                                                           Higgs_mass[0])
            sumWgts += Running_wgt[0]

            NonRunning_wgt[0] = HWWSignalShapes.HiggsCPWeight(opts.mH, 
                                                              Higgs_mass[0], 0)
            sumNonRunning += NonRunning_wgt[0]
            # print eventN, "parti:", parti, 
            # print "id:", genParts.IDUP.at(parti),
            # print "mass:", Higgs_mass[0],
            # print "status:", genParts.ISTUP[parti]

        if (genParts.MOTHUP.at(parti).first >= Higgs_i) and \
                (genParts.MOTHUP.at(parti).second < Higgs_i):
            tmpID = abs(genParts.IDUP.at(parti))
            if not (tmpID in decayCnts):
                decayCnts[tmpID] = 0
            decayCnts[tmpID] += 1
            # print eventN, "parti:", parti, 
            # print "id:", genParts.IDUP.at(parti),
            # print "status:", genParts.ISTUP[parti],
            # print "moms:", genParts.MOTHUP.at(parti).first, genParts.MOTHUP.at(parti).second

    outTree.Fill()
    # if eventN >= 100:
    #     break

print decayCnts
#print "n events:",eventN
print "avg CP weight:",sumWgts/(eventN+1)
print "avg non-running width weight:", sumNonRunning/(eventN+1)

outputFile.Write()
outputFile.Close()
