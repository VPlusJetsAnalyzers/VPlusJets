from optparse import OptionParser
parser = OptionParser()
parser.add_option('-H', '--mH', dest='mH', type='float',
                  help='Higgs Mass Point')
parser.add_option('-u', '--update', dest='update', action='store_true',
                  help='update instead of recreate output file')
parser.add_option('-b', dest='background', action='store_true',
                  help='no windows with plots')
(opts, args) = parser.parse_args()

import pyroot_logon

import ROOT as r
from HWWSignalShapes import HiggsWidth
import numpy as np

fcps = r.TFile("HWW%dNewLHE.root" % opts.mH)
fnocps = r.TFile("HWW%dOldLHE.root" % opts.mH)

cpsTree = fcps.Get("outTree")
nocpsTree = fnocps.Get("outTree")

nocpsTree.Print()

minm = int(min(cpsTree.GetMinimum('Higgs_mass'), 
              nocpsTree.GetMinimum('Higgs_mass'))-0.5)
maxm = int(max(cpsTree.GetMaximum('Higgs_mass'), 
              nocpsTree.GetMaximum('Higgs_mass'))+1.5)
# minm = max(minm, 30)

lowEdge = int(opts.mH - max(2.*HiggsWidth[int(opts.mH)],1.01) + 0.5)
upEdge = int(opts.mH + max(2.*HiggsWidth[int(opts.mH)], 1.01) + 0.5)
print 'mass range:',minm, maxm
print 'signal range:', lowEdge, upEdge


bins = np.arange(minm, lowEdge, (lowEdge - minm)/20.)
bins = np.append(bins, np.arange(lowEdge, upEdge, (upEdge-lowEdge)/50.))
bins = np.append(bins, np.arange(upEdge, maxm, (maxm-upEdge)/20.))
bins = np.append(bins, maxm)

print bins

mww = r.RooRealVar("Higgs_mass", "m_{WW}", bins[0], bins[-1])
mH = r.RooRealVar("mH", "mH", opts.mH, 160., 1200.)
Gamma = r.RooRealVar("Gamma", "#Gamma", HiggsWidth[int(opts.mH)], 0., 200.)
Gamma.setError(Gamma.getVal()*0.1)
mH.setError(mH.getVal()*0.1)
bw = r.RooRelBWRunningWidth("bw", "bw", mww, mH, Gamma)

mH_no = r.RooRealVar('mH_no', "mH", opts.mH, 160., 1200.)
Gamma_no = r.RooRealVar('Gamma_no', "#Gamma", HiggsWidth[int(opts.mH)], 
                        0., 200.)
bw_no = r.RooRelBWRunningWidth("bw_no", "bw_no", mww, mH_no, Gamma_no)

nocps = r.TH1F("nocps_HWW%d" % opts.mH, "nocps", len(bins)-1, bins)
nocps.Sumw2()
cps = r.TH1F("cps_HWW%d" % opts.mH, "cps", len(bins)-1, bins)
cps.Sumw2()
ratio = r.TH1F("ratio_HWW%d" % opts.mH, "ratio", len(bins)-1, bins)
ratio.Sumw2()
# ratioGraph = r.TGraph(len(bins)-1)
# ratioGraph.SetName('pdfRatioGraph')
# ratioGraph.SetLineColor(r.kBlue)
# ratioGraph.SetLineWidth(2)

# wgtVar = "Running_wgt"
# wgtVar = "NonRunning_wgt"
wgtVar = ""
nocpsTree.Draw("Higgs_mass>>+nocps_HWW%d" % opts.mH, wgtVar)
cpsTree.Draw("Higgs_mass>>+cps_HWW%d" % opts.mH)
# r.gPad.Update()
# r.gPad.WaitPrimitive()

dataWeight = r.RooRealVar(wgtVar if len(wgtVar)>0 else 'wgt', 'wgt', 1.0)
# data = r.RooDataHist("data", "data", r.RooArgList(mww), cps)
u_data = r.RooDataSet("u_data", "newCPS", cpsTree, r.RooArgSet(mww))
w_data = r.RooDataSet('w_data', 'oldCPS', nocpsTree, r.RooArgSet(mww,dataWeight), '', wgtVar)
# data.Print()
u_data.Print()
w_data.Print()


datapdfs = [ (u_data, bw), (w_data, bw_no) ]
nocps.Scale(1., 'width')
nocps.Scale(1./nocps.Integral())

cps.Scale(1., 'width')
cps.Scale(1./cps.Integral())

ratio.Divide(cps, nocps)

c1 = r.TCanvas("c1", "no cps")
nocps.Draw()

c2 = r.TCanvas("c2", "cps")
cps.Draw()

c3 = r.TCanvas("c3", "ratio")
ratio.Draw("e")

# fr = bw.fitTo(data, r.RooFit.Save(True),
#               # r.RooFit.Range(200, 908),
#               r.RooFit.PrintLevel(0),
#               )
# fr.SetName("binnedFit_HWW%d" % opts.mH)

# fitplot = mww.frame()
# data.plotOn(fitplot)
# bw.plotOn(fitplot)

# c4 = r.TCanvas("c4", "fit")
# fitplot.Draw()

cans = []
plots = []
frs = []

for (data, pdf) in datapdfs:
    s = 1.0
    fr = None
    mH.setVal(opts.mH)
    Gamma.setVal(HiggsWidth[int(opts.mH)])
    Gamma.setError(Gamma.getVal()*0.1)
    mH.setError(mH.getVal()*0.1)
    mww.setRange('fit', int(mH.getVal()-max(Gamma.getVal()*s,0.51)+ 0.5),
                 int(mH.getVal()+max(Gamma.getVal()*s, 0.51)+0.5))
    fr = pdf.fitTo(data, r.RooFit.Save(True),
                   r.RooFit.Range('fit'),
                   r.RooFit.SumW2Error(True)
                   )

    fr.SetName("unbinnedFit_HWW%d_%s" % (opts.mH,data.GetTitle()))
    frs.append(fr)

    fitplot = mww.frame(lowEdge,
                        upEdge, 50)
    fitplot.SetName("fitplot_HWW%d_%s" % (opts.mH,data.GetTitle()))
    data.plotOn(fitplot)
    pdf.plotOn(fitplot)

    c5 = r.TCanvas('c%s' % data.GetTitle(), data.GetTitle())
    fitplot.Draw()
    cans.append(c5)
    plots.append(fitplot)

nonZero = False
firstNonZero = None
lastNonZero = ratio.GetNbinsX()
for bin in range(1,ratio.GetNbinsX()+1):
    if not nonZero and ratio.GetBinContent(bin) > 0:
        print 'non-zero bin:',bin
        nonZero = True
        if not firstNonZero:
            firstNonZero = bin
    if nonZero and ratio.GetBinContent(bin) == 0:
        print 'zero bin:',bin
        nonZero = False
        lastNonZero = bin-1

print "non-zero bin range:", firstNonZero, lastNonZero
ratio.GetXaxis().SetRange(firstNonZero,lastNonZero)
ratio.Smooth(10, 'R')
c3.Modified()
c3.Update()

sumW = 0.
for event in nocpsTree:
    sumW += ratio.Interpolate(event.Higgs_mass)
print 'sumw',sumW
print 'average cpw:',sumW/nocpsTree.GetEntries()

#     # print 'cps', cps.GetBinContent(bin),'nocps', nocps.GetBinContent(bin),
#     # print 'ratio', ratio.GetBinContent(bin)
#     if cps.GetBinContent(bin) == 0 and nocps.GetBinContent(bin) > 0:
#         ratio.SetBinContent(bin, 0.)
#     mww.setVal(ratio.GetBinCenter(bin))
#     ratioGraph.SetPoint(bin-1, mww.getVal(), bw.getVal(r.RooArgSet(mww))/bw_no.getVal(r.RooArgSet(mww)))

# c3.cd()
# ratio.Draw()
# ratioGraph.Draw('l')
# c3.Modified()
# c3.Update()
    
update = 'recreate'
if opts.update:
    update = 'update'

print "%sing" % update[:-1], "output file"

output = r.TFile("ComplexPoleWeights/CPSHiggs%dShapes.root" % (opts.mH), update)

avgCPS = r.TTree('avgCPS', 'avgCPS')
avgcps = np.array([sumW/nocpsTree.GetEntries()], dtype=np.float32)
print avgcps
avgCPS.Branch('avgcps', avgcps, 'avgcps/F')
avgCPS.Fill()
avgCPS.Write()
nocps.Write()
cps.Write()
ratio.Write()
# ratioGraph.Write()
for fitplot in plots:
    fitplot.Write()
# fr.Write()
for fr in frs:
    fr.Write()
output.ls()
output.Close()
