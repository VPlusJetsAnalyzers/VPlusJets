from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False,
                  help='no X11 windows')
parser.add_option('--bn', dest='bn', default='ShapeParameters',
                  help='basis name for the output files.')
parser.add_option('--basis', dest='basis', type='int',
                  help='basis function')
parser.add_option('--plotName', dest='plotName', 
                  default='fit_mlvjj_plot_with_errors',
                  help='name of RooPlot to pull out of the file')
(opts, args) = parser.parse_args()

import pyroot_logon
import ROOT as r
import pulls

if opts.basis:
    basisModel = opts.basis
else:
    basisModel = args[0]

colorMap = { 
    'nominal'  : r.kBlue,
    'alt' : r.kBlack,
    'poly' : r.kRed+2,
}

curves = []

fbasis = r.TFile.Open('%s_%s.root' % (opts.bn, basisModel))
plot = fbasis.Get(opts.plotName)
plot.Draw()
r.gPad.Update()
plot.Print()
# plot.getCurve('diboson').Print('v')
# basisCurve = plot.getCurve('diboson')
basisCurve = pulls.clipCurve(plot.getCurve('diboson'))
basisError = plot.getCurve('fitErrors')
# topCurve = plot.getCurve('top')
# WpJCurve = plot.getCurve('WpJ')

# basisError.Print()
xmin = plot.GetXaxis().GetXmin()
xmax = plot.GetXaxis().GetXmax()
xtitle = plot.GetXaxis().GetTitle()

varName = 'fit_mlvjj'
wsName = 'w_mWW'
basis_ws = fbasis.Get(wsName)
variable = basis_ws.var(varName)
n_WpJ = basis_ws.var('n_WpJ_mu').getVal()

SignalCenter = basis_ws.var('mean_ggH_fit_mlvjj_core_mu').getVal()
SignalWidth = basis_ws.var('sigma_ggH_fit_mlvjj_core_mu').getVal()
Z = 1.5

SigWindLow = max(SignalCenter - SignalWidth*Z, variable.getMin())
SigWindHigh = min(SignalCenter + SignalWidth*Z, variable.getMax())

print "%.1f sigma signal window: [%.1f, %.1f]" % (Z, SigWindLow, SigWindHigh)
# onlyWpJCurve = pulls.subtractCurves(pulls.clipCurve(WpJCurve),
#                                     pulls.clipCurve(topCurve),
#                                     # xmin-0.8, xmax+0.8
#                                     )

basisNormalCurve = r.RooCurve(basis_ws.pdf('WpJ'), variable, variable.getMin(),
                              variable.getMax(), variable.getBins(), 
                              n_WpJ,
                              r.RooArgSet(variable))
basisNormalCurve = pulls.clipCurve(basisNormalCurve)

basisSigYield = basisNormalCurve.average(SigWindLow,SigWindHigh)
# variable.setRange('SignalWidowRange', SigWindLow, SigWindHigh)
# basisIntegral = basis_ws.pdf('WpJ').createIntegral(r.RooArgSet(variable),
#                                                    r.RooArgSet(variable),
#                                                    "SignalWindowRange")
# overallScale = basisNormalCurve.interpolate(xmin+10)/onlyWpJCurve.interpolate(xmin+10)
overallScale = 1./plot.getFitRangeBinW()

print 'n_WpJ:', n_WpJ, 'overall scale:', overallScale
sigDiffs = []

envelopeCurve = pulls.subtractCurves(basisError,basisCurve,
                                     xmin-0.8,xmax+0.8,
                                     debug=False, scale = overallScale)
basisDataHist = basis_ws.pdf('WpJ').generateBinned(r.RooArgSet(variable),
                                                   n_WpJ+0.5, True)
# basisDataHist.Print()
basisHist = basisDataHist.createHistogram('basisHist', basis_ws.var(varName))
# basisHist.Print()

for model in args:
    f = r.TFile.Open('%s_%s.root' % (opts.bn, model))
    ws_ = f.Get(wsName)
    var_ = ws_.var(varName)
    modelCurve = r.RooCurve(ws_.pdf('WpJ'), var_, variable.getMin(),
                            variable.getMax(), variable.getBins(), 
                            n_WpJ,
                            r.RooArgSet(var_))
    modelDataHist = ws_.pdf('WpJ').generateBinned(r.RooArgSet(ws_.var(varName)),
                                                  n_WpJ+0.5, True)
    # modelDataHist.Print()
    modelHist = modelDataHist.createHistogram("%sHist" % model, 
                                              ws_.var(varName))
    # modelHist.Print()
    model_pval = modelHist.Chi2Test(basisHist, 'WWP')
    modelYield = modelCurve.average(SigWindLow,SigWindHigh)
    curves.append(pulls.subtractCurves(modelCurve, basisNormalCurve, 
                                       xmin-0.8, xmax+0.8
                                       ))
    # modelCurve = pulls.clipCurve(f.Get(opts.plotName).getCurve('diboson'))
    # curves.append(pulls.subtractCurves(modelCurve, basisCurve,
    #                                    xmin-0.8,xmax+0.8, debug=False,
    #                                    scale=overallScale))
    curves[-1].SetLineColor(colorMap[model])
    curves[-1].SetName('p_%s_minus_%s' % (model, basisModel))
    curves[-1].SetTitle('%s minus %s' % (model, basisModel))
    sigDiffs.append((modelYield-basisSigYield)*(SigWindHigh-SigWindLow))
    print "Signal window difference %s - %s: %.2f" % (model, basisModel,
                                                      sigDiffs[-1])
    f.Close()

sigCurve = pulls.scaleCurve(plot.getCurve('signal_HWW'), 0.5*overallScale)
sigCurve.SetLineColor(r.kViolet)
sigCurve.SetTitle("Higgs Signal")
curves.append(sigCurve)
HiggsSigYield = sigCurve.average(SigWindLow,SigWindHigh)*(SigWindHigh-SigWindLow)
sigDiffs.append(HiggsSigYield)
print "Signal window Higgs: %.2f +/- %.2f" % (HiggsSigYield, 
                                              r.TMath.Sqrt(HiggsSigYield))

leg = r.TLegend(0.6, 0.95, 0.95, 0.75, "", "NDC")

envelopeCurve.Draw('af')
envelopeCurve.GetXaxis().SetLimits(xmin, xmax)
envelopeCurve.GetXaxis().SetTitle(xtitle)
envelopeCurve.GetYaxis().SetTitle('Events / GeV')
leg.AddEntry(envelopeCurve, "", "f")

print sigDiffs
for i,curve in enumerate(curves):
    curve.Draw('l')
    leg.AddEntry(curve, "", 'l')
    print curve.GetTitle(),
    print 'bias relative to Higgs yield in a %.1f sigma window: %.3f' % (Z, sigDiffs[i]/HiggsSigYield)

leg.Draw()
r.gPad.Update()
r.gPad.Print("%s_comparisons%s.png" % (opts.bn, basisModel))
