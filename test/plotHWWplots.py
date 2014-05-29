from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False,
                  help='no X11 windows')
parser.add_option('-p', action='store_false', dest='prelim', default=True,
                  help='turn off the preliminary portion of the CMS label')
parser.add_option('-s', dest='scaler', type=float,
                  help='scale to multiply the signal by')
(opts, args) = parser.parse_args()

import pyroot_logon
import re
from ROOT import TFile, gPad, TCanvas, TF2, SetOwnership

cs = TCanvas('cs', 'stacked')
cp = TCanvas('cp', 'pull')
cp.SetGridy()

scaleFunction = None
for fname in args:
    print fname
    fname_minusPath = fname.split('/')[-1]
    fname_path = '/'.join(fname.split('/')[0:-1])
    fname_parts = fname_minusPath.split('_')
    if len(fname_path) > 0:
        fname_path += '/'
    print 'filename:', fname_minusPath, 'path:', fname_path

    f = TFile(fname)
    match = re.search(r'\d+', fname_parts[0])
    mH = int(match.group(0))
    mjj_plot = f.Get('Mass2j_PFCor_stacked')
    sigCurve = mjj_plot.getCurve('signal_HWW')
    if sigCurve and (opts.scaler != None):
        if scaleFunction == None:
            scaleFunction = TF2("scaleFunction", "y*%.3f" % opts.scaler,
                                mjj_plot.GetXaxis().GetXmin(),
                                mjj_plot.GetXaxis().GetXmax(),
                                mjj_plot.GetMinimum(), mjj_plot.GetMaximum())
        sigCurve.Apply(scaleFunction)
        sigEntry = mjj_plot.findObject('theLegend').GetListOfPrimitives().Last().SetLabel('H(%i)#times%.0f' %(mH,opts.scaler*2))
#         SetOwnership(sigEntry, False)
# #         print sigEntry
# #         sigEntry.Print()
#         sigEntry.SetLabel('H(%i)#times%.0f' %(mH,opts.scaler*2))
        
        
    cs.cd()
    cs.SetLogy(False)
    mjj_plot.Draw()
    pyroot_logon.cmsLabel(cs, lumi = 19.2, prelim = opts.prelim)
    gPad.Update()
    gPad.Print('%s%s_%s_mjj_stacked.pdf' % (fname_path,fname_parts[0],
                                            fname_parts[1]))
    gPad.Print('%s%s_%s_mjj_stacked.png' % (fname_path,fname_parts[0],
                                            fname_parts[1]))

    mjj_pull = f.Get('Mass2j_PFCor_pull')
    cp.cd()
    mjj_pull.Draw('ap')
    pyroot_logon.cmsLabel(cp, lumi = 19.2, prelim = opts.prelim)
    gPad.Update()
    gPad.Print('%s%s_%s_mjj_pull.pdf' % (fname_path,fname_parts[0],
                                         fname_parts[1]))
    gPad.Print('%s%s_%s_mjj_pull.png' % (fname_path,fname_parts[0],
                                         fname_parts[1]))

    mWW_plot = f.Get('fit_mlvjj_plot_stacked')
    sigCurve = mWW_plot.getCurve('signal_HWW')
    if sigCurve and (opts.scaler != None):
        sigCurve.Apply(scaleFunction)
        sigEntry = mWW_plot.findObject('theLegend').GetListOfPrimitives().Last().SetLabel('H(%i)#times%.0f' %(mH,opts.scaler*2))
#         SetOwnership(sigEntry, False)       
# #         print sigEntry
# #         sigEntry.Print()
#         sigEntry.SetLabel('H(%i)#times%.0f' %(mH,opts.scaler*2))
    cs.cd()
    mWW_plot.Draw()
    pyroot_logon.cmsLabel(cs, lumi = 19.2, prelim = opts.prelim)
    gPad.Update()
    gPad.Print('%s%s_%s_mWW_stacked.pdf' % (fname_path,fname_parts[0],
                                            fname_parts[1]))
    gPad.Print('%s%s_%s_mWW_stacked.png' % (fname_path,fname_parts[0],
                                            fname_parts[1]))
    if mH > 400:
        cs.SetLogy(True)
        mWW_plot.SetAxisRange(0.1, 1e5, 'Y')
        gPad.Update()
        gPad.Print('%s%s_%s_mWW_stacked_log.pdf' % \
                       (fname_path,fname_parts[0], fname_parts[1]))
        gPad.Print('%s%s_%s_mWW_stacked_log.png' % \
                       (fname_path,fname_parts[0], fname_parts[1]))


    mWW_pull = f.Get('mWW_pull')
    cp.cd()
    mWW_pull.Draw('ap')
    pyroot_logon.cmsLabel(cp, lumi = 19.2, prelim = opts.prelim)
    gPad.Update()
    gPad.Print('%s%s_%s_mWW_pull.pdf' % (fname_path, fname_parts[0], 
                                         fname_parts[1]))
    gPad.Print('%s%s_%s_mWW_pull.png' % (fname_path, fname_parts[0], 
                                         fname_parts[1]))
