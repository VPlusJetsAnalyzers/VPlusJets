from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False,
                  help='no X11 windows')
parser.add_option('--bn', dest='bn', default='ShapeParameters',
                  help='basis name for the output files.')
parser.add_option('--basis', dest='basis', type='int',
                  help='basis function')
parser.add_option('--plotName', dest='plotName', default='fit_mlvjj_Plot',
                  help='name of RooPlot to pull out of the file')
(opts, args) = parser.parse_args()

import pyroot_logon
import ROOT as r
import pulls

args = [ int(x) for x in args ]

if opts.basis:
    basisModel = opts.basis
else:
    basisModel = args[0]

colorMap = { 
    8  : r.kBlue,
    10 : r.kGreen+2,
    12 : r.kOrange+2,
    14 : r.kBlack,
    23 : r.kRed+2,
    33 : r.kViolet+1,
    35 : r.kYellow+3,
    37 : r.kSpring-1,
}

nameMap = {
    8  : "2",
    10 : "5",
    12 : "6",
    14 : "1",
    23 : "poly",
    33 : "4",
    35 : "4",
    37 : "3",

}             

curves = []

fbasis = r.TFile.Open('%s_%i.root' % (opts.bn, basisModel))
fbasis.Get('fit_mlvjj_Plot').Draw()
r.gPad.Update()
basisCurve = fbasis.Get('fit_mlvjj_Plot').getCurve('fitCurve')
basisError = fbasis.Get('fit_mlvjj_Plot').getCurve('fitErrors')
xmin = fbasis.Get('fit_mlvjj_Plot').GetXaxis().GetXmin()
xmax = fbasis.Get('fit_mlvjj_Plot').GetXaxis().GetXmax()
ytitle = fbasis.Get('fit_mlvjj_Plot').GetYaxis().GetTitle()

envelopeCurve = pulls.subtractCurves(basisError, basisCurve)

for model in args:
    f = r.TFile.Open('%s_%i.root' % (opts.bn, model))
    modelCurve = f.Get('fit_mlvjj_Plot').getCurve('fitCurve')
    curves.append(pulls.subtractCurves(modelCurve, basisCurve))
    curves[-1].SetLineColor(colorMap[model])
    curves[-1].SetName('p_%s_minus_%s' % (nameMap[model], nameMap[basisModel]))
    curves[-1].SetTitle('%s minus %s' % (nameMap[model], nameMap[basisModel]))
    f.Close()

leg = r.TLegend(0.6, 0.95, 0.95, 0.75, "", "NDC")

envelopeCurve.Draw('af')
envelopeCurve.GetXaxis().SetLimits(xmin, xmax)
envelopeCurve.GetXaxis().SetTitle('m_{l#nujj} (GeV)')
envelopeCurve.GetYaxis().SetTitle(ytitle)
leg.AddEntry(envelopeCurve, "", "f")

for curve in curves:
    curve.Draw('l')
    leg.AddEntry(curve, "", 'l')

leg.Draw()
r.gPad.Update()
r.gPad.Print("%s_comparisons%i.png" % (opts.bn, basisModel))
