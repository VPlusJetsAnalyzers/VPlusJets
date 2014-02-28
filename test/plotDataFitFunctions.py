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
# plot.Print()
# plot.getCurve('diboson').Print('v')
basisCurve = pulls.clipCurve(plot.getCurve('diboson'))
basisError = plot.getCurve('fitErrors')
# basisError.Print()
xmin = plot.GetXaxis().GetXmin()
xmax = plot.GetXaxis().GetXmax()
xtitle = plot.GetXaxis().GetTitle()
ytitle = plot.GetYaxis().GetTitle()

envelopeCurve = pulls.subtractCurves(basisError,basisCurve,xmin-0.8,xmax+0.8,
                                     debug=False)

for model in args:
    f = r.TFile.Open('%s_%s.root' % (opts.bn, model))
    modelCurve = pulls.clipCurve(f.Get(opts.plotName).getCurve('diboson'))
    curves.append(pulls.subtractCurves(modelCurve, basisCurve,
                                       xmin-0.8,xmax+0.8, debug=False))
    curves[-1].SetLineColor(colorMap[model])
    curves[-1].SetName('p_%s_minus_%s' % (model, basisModel))
    curves[-1].SetTitle('%s minus %s' % (model, basisModel))
    f.Close()

leg = r.TLegend(0.6, 0.95, 0.95, 0.75, "", "NDC")

envelopeCurve.Draw('af')
envelopeCurve.GetXaxis().SetLimits(xmin, xmax)
envelopeCurve.GetXaxis().SetTitle(xtitle)
envelopeCurve.GetYaxis().SetTitle(ytitle)
leg.AddEntry(envelopeCurve, "", "f")

for curve in curves:
    curve.Draw('l')
    leg.AddEntry(curve, "", 'l')

leg.Draw()
r.gPad.Update()
r.gPad.Print("%s_comparisons%s.png" % (opts.bn, basisModel))
