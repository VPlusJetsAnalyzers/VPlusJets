from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False,
                  help='no X11 windows')
parser.add_option('-m', '--mode', default="HWW1D2FitsConfig", 
                  dest='modeConfig',
                  help='which config to select look at HWW2DConfig.py for ' +\
                      'an example.  Use the file name minus the .py extension.'
                  )
parser.add_option('-H', '--mH', dest='mH', default=350, type='int',
                  help='Higgs Mass Point')
parser.add_option('-j', '--Njets', dest='Nj', default=2, type='int',
                  help='Number of jets.')
parser.add_option('--electrons', dest='isElectron', action='store_true',
                  default=False, help='do electrons instead of muons')
parser.add_option('--nosig', dest='includeSignal', action='store_false',
                  default=True, help='include signal shape in fit')
parser.add_option('--toy', dest='toy', action='store_true',
                  help='use pseudo-data instead of data file')
parser.add_option('--toyOut', dest='toyOut', help='filename for toy output')
parser.add_option('--seed', dest='seed', type='int', help='random seed')
parser.add_option('--ws', dest='ws', 
                  help='filename that contains workspace to be used cloned ' +\
                      'for use')
parser.add_option('--debug', dest='debug', action='store_true', default=False,
                  help='turn on extra debugging information')
parser.add_option('--limit', dest='doLimit', type='int',
                  default=0,
                  help='calculate expected limits using profile likelihood, ' +\
                      'number of toys to use.')
parser.add_option('--obsLimit', dest='obsLimit', action='store_true',
                  default=False,
                  help='calculate observed limit too')
parser.add_option('--mva', dest='mvaCut', type='float',
                  help='override cut value for mva')
parser.add_option('--xrootd', dest='xrootd', action='store_true',
                  help='use xrootd file opening.')
parser.add_option('--injectS', type='float', dest='sigInject',
                  help='amount of signal to inject')
parser.add_option('--genConfig',
                  dest='genConfig',
                  help='which config to select look at HWW2DConfig.py for ' +\
                      'an example.  Use the file name minus the .py extension.'
                  )
parser.add_option('--xc', dest='xc',
                  help='use cross-check background to generate')
parser.add_option('--sideband', dest='sb', type='int',
                  help='use sideband dataset and model instead')
parser.add_option('--unbinned', dest='binned', action='store_false',
                  help='unbinned m_lvjj fit instead of binned ML.')
parser.add_option('--qRooFit', dest='qRF', action='store_true',
                  help='quiet RooFit warnings.')

(opts, args) = parser.parse_args()

import pyroot_logon

#import HWW2DConfig
config = __import__(opts.modeConfig)
import RooWjj2DFitter

from ROOT import TCanvas, RooFit, RooLinkedListIter, TMath, RooRandom, TFile, \
    RooDataHist, RooMsgService, TStopwatch, RooAbsPdf, TBox, kBlack, kRed, \
    kBlue, kOrange, kDashed, RooAbsCollection, RooArgSet, RooStats, Double, \
    RooAbsReal, RooAbsData, TH1F, RooArgList
import pulls

timer = TStopwatch()
timer.Start()

#RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9)
#RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9)
if not opts.debug:
    if opts.qRF:
        RooMsgService.instanc().setGlobalKillBelow(RooFit.FATAL)
    else:
        RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)
    # RooMsgService.instance().addStream(RooFit.ERROR,
    #                                    RooFit.Prefix(True), 
    #                                    RooFit.ClassName('RooExpPoly'),
    #                                    RooFit.OutputFile('/dev/null'))
    # RooMsgService.instance().Print('v')

if hasattr(opts, "seed") and (opts.seed >= 0):
    print "random seed:", opts.seed
    RooRandom.randomGenerator().SetSeed(opts.seed)

mvaCutOverride = None
if hasattr(opts, "mvaCut"):
    mvaCutOverride = opts.mvaCut

mjjArgs = []
sideArgs = []
mWWArgs = []
for arg in args:
    if arg[-8:] == '_mjj.txt':
        mjjArgs.append(arg)
    elif arg[-8:] == '_mWW.txt':
        mWWArgs.append(arg)
    else:
        mWWArgs.append(arg)

print mjjArgs

pars = config.theConfig(Nj = opts.Nj, mH = opts.mH, 
                        isElectron = opts.isElectron, initFile = mjjArgs,
                        includeSignal = False, MVACutOverride = mvaCutOverride,
                        xrootd = opts.xrootd)

if opts.toy:
    pars.blind = False

fitter = RooWjj2DFitter.Wjj2DFitter(pars)
fitter.ws.SetName("w_mjj")

totalPdf = fitter.makeFitter()
sigScaler = 2.
mode = 'muon'
if opts.isElectron:
    mode = 'electron'
fr = None

chi2 = 0
ndf = 1

def mjjfitting(fitter):

    if opts.ws:
        fitter.loadWorkspaceFromFile(opts.ws, getFloatPars = False,
                                     wsname = 'w_mjj')

    fitter.readParametersFromFile()
    fitter.expectedFromPars()
    fitter.resetYields()

    startpars = totalPdf.getParameters(fitter.ws.set('obsSet'))
    fitter.ws.defineSet("params", startpars)
    fitter.ws.saveSnapshot("initPars", startpars)

    if opts.toy and not opts.ws:
        #generate toy dataset
        print 'Generated parameters'
        fitter.ws.set('params').Print('v')
        fitter.ws.saveSnapshot("genPars", startpars)

        data = totalPdf.generate(fitter.ws.set('obsSet'), RooFit.Name('data_obs'),
                                 RooFit.Extended())
        if fitter.pars.binData:
            data = RooDataHist('data_obs', 'data_obs', fitter.ws.set('obsSet'),
                               data)
            data.Print('v')
        getattr(fitter.ws, 'import')(data)
    else:
        data = fitter.loadData()

    data.Print()
    startpars.IsA().Destructor(startpars)

    print 'Time elapsed: %.1f sec' % timer.RealTime()
    print 'CPU time used: %.1f sec' % timer.CpuTime()
    print 'starting fitting routine'
    timer.Continue()
    fitter.ws.var('diboson_nrm').setConstant()
    fitter.ws.var('top_nrm').setConstant()
    fitter.ws.var('r_signal').setVal(1.0)
    fitter.ws.var('r_signal').setError(0.04)

    fr = fitter.fit(overrideRangeCmd=True)

    plot1 = fitter.stackedPlot(pars.var[0])

    sigPdf = fitter.ws.pdf('ggH')
    n_sig = fitter.ws.var('n_ggH').getVal() + fitter.ws.var('n_qqH').getVal()
    print "N HWW:", n_sig

    sigPdf.plotOn(plot1, RooFit.Range('plotRange'),
                  RooFit.NormRange('plotRange'),
                  RooFit.Normalization(n_sig*sigScaler, RooAbsReal.NumEvent),
                  RooFit.Name('signal_HWW'),
                  RooFit.LineColor(kBlue+2))
    plot1.getCurve().SetTitle("H(%i)#times%.0f" % (opts.mH, sigScaler))

    leg1 = RooWjj2DFitter.Wjj2DFitter.legend4Plot(plot1)

    c1 = TCanvas('c1', pars.var[0] + ' plot')
    plot1.addObject(leg1)
    plot1.Draw()
    c1.Update()

    ndf = 0

    if fr:
        floatVars = [ fr.floatParsFinal().at(i).GetName() \
                          for i in range(0, fr.floatParsFinal().getSize()) ]
        fitter.ws.defineSet('floatingParams', ','.join(floatVars))
        fitter.ws.saveSnapshot("fitPars", ','.join(floatVars))
        ndf = fr.floatParsFinal().getSize() - \
            fitter.ws.set('constraintSet').getSize()
        # fr.Print('v')

    print '%i free parameters in the fit' % ndf

    firstCurve1 = plot1.getObject(1)
    # firstCurve1.Print()

    (chi2_1, ndf_1) = pulls.computeChi2(plot1.getHist('theData'), firstCurve1)
    pull1 = pulls.createPull(plot1.getHist('theData'), firstCurve1)


    chi2 = chi2_1
    ndf = ndf_1-ndf

    cp1 = TCanvas("cp1", pars.var[0] + ' pull')
    pull1.Draw('ap')
    pull1.SetName(pars.var[0] + "_pull")
    cp1.SetGridy()
    cp1.Update()
    pull1.GetXaxis().SetLimits(pars.varRanges[pars.var[0]][1], 
                               pars.varRanges[pars.var[0]][2])
    pull1.GetXaxis().SetTitle(fitter.ws.var(pars.var[0]).getTitle(True).Data())
    pull1.GetYaxis().SetTitle("pull (#sigma)")

    blinder = plot1.findObject('TBox')
    blinder2 = None
    if blinder:
        blinder2 = TBox(pars.exclude[pars.var[0]][0], pull1.GetYaxis().GetXmin(),
                        pars.exclude[pars.var[0]][1], pull1.GetYaxis().GetXmax())
        blinder2.SetFillColor(blinder.GetFillColor())
        blinder2.SetFillStyle(blinder.GetFillStyle())

    if blinder2:
        #blinder2.Print()
        blinder2.Draw()

    cp1.Update()

    print 'Time elapsed: %.1f sec' % timer.RealTime()
    print 'CPU time used: %.1f sec' % timer.CpuTime()

    timer.Continue()

    return (plot1, pull1, chi2, ndf, fr)

if not opts.toy:
    (plot1, pull1, chi2, ndf, fr) = mjjfitting(fitter)
else:
    fitter.readParametersFromFile()
    fitter.expectedFromPars()
    fitter.resetYields()
    fitter.ws.var('WpJ_nrm').setError(0.006)


import HWW1D2FitsConfig_mWW
print mWWArgs

configArgs = {
    'Nj': opts.Nj,
    'mH': opts.mH,
    'isElectron': opts.isElectron,
    'initFile': mWWArgs,
    'includeSignal': True,
    'mjj_config': opts.modeConfig,
    'MVACutOverride': mvaCutOverride,
    'xrootd': opts.xrootd,
    }
if opts.sb > 0:
    sb = 'high'
    configArgs['sideband'] = sb
elif opts.sb < 0:
    sb = 'low'
    configArgs['sideband'] = sb
if opts.binned != None:
    configArgs['binned'] = opts.binned
# pars_mWW = HWW1D2FitsConfig_mWW.theConfig(**configArgs)
# pars_mWW.yieldConstraints['WpJ'] = fitter.ws.var('WpJ_nrm').getError()

# sb = None

def prepFitter(config = opts.modeConfig, initFiles = mWWArgs):
    theArgs = dict(configArgs)
    theArgs['mjj_config'] = config
    theArgs['initFile'] = initFiles    
    pars_mWW = HWW1D2FitsConfig_mWW.theConfig(**theArgs)
    pars_mWW.yieldConstraints['WpJ'] = fitter.ws.var('WpJ_nrm').getError()

    if opts.binned == False:
        pars_mWW.constrainShapes = []

    if opts.toy or opts.sb:
        pars_mWW.blind = False

    # add systematic errors multipliers to ggH signals
    if (opts.mH >= 400):
        pars_mWW.ggHdoSystMult = True

    fitter_mWW = RooWjj2DFitter.Wjj2DFitter(pars_mWW)
    fitter_mWW.ws.SetName("w_mWW")
    totalPdf_mWW = fitter_mWW.makeFitter()
    fitter_mWW.readParametersFromFile()

    regionName = 'signalRegion'
    if opts.sb:
        regionName = 'sidebandRegion'
        fitter.ws.var(fitter.pars.var[0]).setRange(regionName,
                                                   pars_mWW.regionLow,
                                                   pars_mWW.regionHigh)
    sigInt = fitter.ws.pdf('WpJ').createIntegral(fitter.ws.set('obsSet'),
                                                 fitter.ws.set('obsSet'),
                                                 regionName)
    fullInt = fitter.ws.pdf('WpJ').createIntegral(fitter.ws.set('obsSet'),
                                                  fitter.ws.set('obsSet'))
    n_WpJ_sig = sigInt.getVal()/fullInt.getVal()*fitter.ws.function('f_WpJ_norm').getVal()
    print 'W+jets integral in signal region:',sigInt.getVal()/fullInt.getVal(),
    print 'W+jets yield in the signal region:',n_WpJ_sig

    WpJ_mWW = fitter_mWW.ws.var('n_WpJ')
    WpJ_mWW.setVal(n_WpJ_sig)
    fitter_mWW.expectedFromPars()
    fitter_mWW.resetYields()

    if opts.sigInject:
        fitter_mWW.ws.var('r_signal').setVal(opts.sigInject)
    fitter_mWW.ws.var('r_signal').setError(0.1)
    fitter_mWW.ws.var('r_signal').setRange(-3., 9.)
    fitter_mWW.ws.var('r_signal').setConstant(False)

    return fitter_mWW, pars_mWW

(fitter_mWW,pars_mWW) = prepFitter()
if opts.toy:
    pars_mWW.blind = False
totalPdf_mWW = fitter_mWW.makeFitter()

fitter_gen = fitter_mWW
full_pdf = fitter_mWW.makeConstrainedFitter()
genPdf = full_pdf
if opts.toy and opts.genConfig:
    iFiles = mWWArgs
    if opts.xc:
        iFiles = [ 
            opts.xc if (fn.find('WpJ')>=0) else fn for fn in mWWArgs ]
    (fitter_gen,pars_gen) = prepFitter(opts.genConfig, iFiles)
    genPdf = fitter_gen.makeConstrainedFitter()

mWWCut = '((%s>%.0f)&&(%s<%.0f))' % \
    (pars.var[0], pars.exclude[pars.var[0]][0],
     pars.var[0], pars.exclude[pars.var[0]][1])
print 'signal region cut:',mWWCut

params_mWW = totalPdf_mWW.getParameters(fitter_mWW.ws.set('obsSet'))

predictedPars = params_mWW.snapshot()

params_mWW.Print("v")
fitter_mWW.ws.defineSet("params", params_mWW)

# if opts.ws:
#     fitter_mWW.loadWorkspaceFromFile(opts.ws, getFloatPars = False,
#                                      wsname = 'w_mWW')

fitter_mWW.ws.saveSnapshot("genPars", params_mWW)
if opts.toy and not opts.ws:
    print "Generate global nuisance parameters"
    constraints = genPdf.getAllConstraints(fitter_gen.ws.set('obsSet'),
                                           fitter_gen.ws.set('constrainedSet'))
    # nuisIter = fitter_gen.ws.set('constraintSet').createIterator()
    nuisIter = constraints.createIterator()
    c = nuisIter.Next()
    while c:
        obs = c.getObservables(fitter_gen.ws.set('constrainedSet'))
        c.Print()
        if obs.selectByName('*_nrm').getSize() == 1:
            print "generating for %s" % obs.first().GetName()
            nuisData = c.generate(obs, RooFit.NumEvents(1))
            # nuisData.Print('v')
            fitter_gen.ws.var(nuisData.get().first().GetName()).setVal(nuisData.get().first().getVal())
            nuisData.IsA().Destructor(nuisData)
            obs.Print('v')
        obs.IsA().Destructor(obs)
        c = nuisIter.Next()
    constraints.IsA().Destructor(constraints)
    #generate toy dataset
    print 'Generated parameters'
    params_gen = genPdf.getParameters(fitter_gen.ws.set('obsSet'))
    params_gen.Print('v')
    if not fitter_gen.ws.loadSnapshot('genPars'):
        fitter_gen.ws.saveSnapshot('genPars', params_gen)

    data = genPdf.generate(fitter_gen.ws.set('obsSet'),
                           RooFit.Name('data_obs'),
                           RooFit.Extended())
    if fitter_mWW.pars.binData:
        data = RooDataHist('data_obs', 'data_obs', fitter_mWW.ws.set('obsSet'),
                           data)
    data.Print('v')
    getattr(fitter_mWW.ws, 'import')(data)
else:
    fitter_mWW.loadDataFromWorkspace(fitter.ws, mWWCut)

#compute limits
import limits
from array import array

upperHist = None

if not full_pdf:
    full_pdf = totalPdf_mWW

fitter_mWW.ws.var('r_signal').setVal(0.)

if opts.doLimit:
    # parIter = params_mWW.createIterator()
    # p = parIter.Next()
    # while p:
    #     if p.GetName()[-4:] == '_nrm':
    #         p.setVal(1.0)
    #         p.setRange(-1., 5.)
    #     p = parIter.Next()
    (expectedLimit, toys) = \
                    limits.expectedPlcLimit(fitter_mWW.ws.var(pars_mWW.var[0]),
                                            fitter_mWW.ws.var('r_signal'),
                                            full_pdf, fitter_mWW.ws,
                                            ntoys = opts.doLimit,
                                            binData = pars_mWW.binData)

    upperHist = TH1F('upperHist', 'upper limit hist', 50,
                     fitter_mWW.ws.var('r_signal').getMin(),
                     fitter_mWW.ws.var('r_signal').getMax())
    uppers = []
    for toy in toys:
        #print toy
        if (toy['r_signal']['ok']) and \
           (toy['r_signal']['upper'] < (fitter_mWW.ws.var('r_signal').getMax()-0.02)):
            upperHist.Fill(toy['r_signal']['upper'])
            uppers.append(toy['r_signal']['upper'])
            
    qs = array('d', [0.]*5)
    probs = array('d', [0.022, 0.16, 0.5, 0.84, 0.978])
    uppers.sort()
    uppersArray = array('d', uppers)
    #upperHist.Print()
        
    print 'expected 95%% CL upper limit: %0.4f +/- %0.4f' % \
          (expectedLimit['upper'], expectedLimit['upperErr'])
    c_upper = TCanvas('c_upper', 'toy upper limits')
    upperHist.Draw()
    c_upper.Update()

    TMath.Quantiles(len(uppers), len(qs), uppersArray, qs, probs)
    # nquants = upperHist.GetQuantiles(len(qs), qs, probs)
    print 'sensible expected 95%% CL upper limit quantiles for %i toys: [' % len(uppers),
    for q in qs:
        print '%.4f,' % q,
    print ']'
    print 'expected 95%% CL lower limit: %0.4f +/- %0.4f' % \
          (expectedLimit['lower'], expectedLimit['lowerErr'])

fitter_mWW.ws.var('r_signal').setVal(0.0)
fitter_mWW.ws.var('r_signal').setConstant(True)
if opts.sigInject != None:
    fitter_mWW.ws.var('r_signal').setConstant(False)
fitter_mWW.ws.var('r_signal').setError(0.1)
# fitter_mWW.ws.var('r_signal').setRange(-0.5, 5.)

fr_mWW = None
fr_mWW = fitter_mWW.fit()

# if fr_mWW.statusCodeHistory(0) != 0:
#     print fr_mWW.statusLabelHistory(0), 'failed'
#     fr_mWW = fitter_mWW.fit(True)

if opts.toyOut:
    outFile = open(opts.toyOut, 'a', 1)
    fitter_mWW.ws.loadSnapshot("genPars")

    line = '%s %.6g %.6g %.6g '
    outFile.write('nll %f covQual %i edm %.4g ' % (fr_mWW.minNll(), 
                                                   fr_mWW.covQual(),
                                                   fr_mWW.edm())
                  )
    # outFile.write('chi2Prob %f ' % (TMath.Prob(chi2, ndf)))
    for cycle in range(0, fr_mWW.numStatusHistory()):
        outFile.write('%s %i ' % (fr_mWW.statusLabelHistory(cycle),
                                  fr_mWW.statusCodeHistory(cycle))
                      )
    for i in range(0, fr_mWW.floatParsFinal().getSize()):
        outFile.write(line % (fr_mWW.floatParsFinal().at(i).GetName(),
                              fr_mWW.floatParsFinal().at(i).getVal(),
                              fr_mWW.floatParsFinal().at(i).getError(),
                              fitter_mWW.ws.var(fr_mWW.floatParsFinal().at(i).GetName()).getVal()
                              )
                      )
    outFile.write('\n')
    outFile.close()

plot_mWW = fitter_mWW.stackedPlot(pars_mWW.var[0])
plot_mWW.SetName('%s_plot_stacked' % (pars_mWW.var[0]))
sigPdf = fitter_mWW.ws.pdf('ggH')
n_sig = fitter_mWW.ws.var('n_ggH').getVal() + fitter_mWW.ws.var('n_qqH').getVal()
# sigScaler = 2.
print "N HWW:", n_sig

sigPdf.plotOn(plot_mWW, RooFit.Range('plotRange'),
              RooFit.NormRange('plotRange'),
              RooFit.Normalization(n_sig*sigScaler, RooAbsReal.NumEvent),
              RooFit.Name('signal_HWW'),
              RooFit.LineColor(kBlue+2))
plot_mWW.getCurve().SetTitle("H(%i)#times%.0f" % (opts.mH, sigScaler))

leg_mWW = RooWjj2DFitter.Wjj2DFitter.legend4Plot(plot_mWW)
plot_mWW.addObject(leg_mWW)

plot_mWW_withErrs = fitter_mWW.stackedPlot(pars_mWW.var[0])
plot_mWW_withErrs.SetName('%s_plot_with_errors' % (pars_mWW.var[0]))
nexp = totalPdf_mWW.expectedEvents(fitter_mWW.ws.set('obsSet'))
if fr_mWW:
    totalPdf_mWW.plotOn(plot_mWW_withErrs, 
                        RooFit.VisualizeError(fr_mWW, 1, True),
                        RooFit.Range('plotRange'),
                        RooFit.NormRange('plotRange'),
                        RooFit.Normalization(nexp, RooAbsReal.NumEvent),
                        RooFit.Name('fitErrors'),
                        RooFit.FillColor(kOrange+1),
                        RooFit.FillStyle(3001))
    plot_mWW_withErrs.getCurve().SetTitle('Fit errors')
    # errs = plot_mWW_withErrs.getCurve('fitErrors')
    # (upper, lower) = pulls.splitErrCurve(errs)

sigPdf.plotOn(plot_mWW_withErrs, RooFit.Range('plotRange'),
              RooFit.NormRange('plotRange'),
              RooFit.Normalization(n_sig*sigScaler, RooAbsReal.NumEvent),
              RooFit.Name('signal_HWW'),
              RooFit.LineColor(kBlue+2))
plot_mWW_withErrs.getCurve().SetTitle("H(%i)#times%.0f" % (opts.mH, sigScaler))

leg_mWW_withErrs = RooWjj2DFitter.Wjj2DFitter.legend4Plot(plot_mWW_withErrs)
plot_mWW_withErrs.addObject(leg_mWW_withErrs)

# errs.Print('v')
# upper.Print('v')
# lower.Print('v')

c_mWW = TCanvas('c_mWW', pars_mWW.var[0] + ' plot')
plot_mWW.Draw()
# if blinder:
#     plot_mWW.setInvisible('theData', True)
c_mWW.Update()

c_mWW_err = TCanvas('c_mWW_err', 'with errors')
plot_mWW_withErrs.Draw()
c_mWW_err.Update()

(chi2_mWW, ndf_mWW) = pulls.computeChi2(plot_mWW.getHist('theData'), 
                                        plot_mWW.getObject(1))
nFreeFitParams = 0
if fr_mWW:
    nFreeFitParams = fr_mWW.floatParsFinal().getSize() - \
        fitter_mWW.ws.set('constraintSet').getSize()
    ndf_mWW -= nFreeFitParams
pull_mWW = pulls.createPull(plot_mWW.getHist('theData'), 
                            plot_mWW.getObject(1))
cp_mWW = TCanvas('cp_mWW', pars_mWW.var[0] + ' mWW pull')
pull_mWW.SetName('mWW_pull')
pull_mWW.Draw('ap')
cp_mWW.SetGridy()
cp_mWW.Update()
pull_mWW.GetXaxis().SetLimits(
    pars_mWW.varRanges[pars_mWW.var[0]][1], 
    pars_mWW.varRanges[pars_mWW.var[0]][2]
    )
pull_mWW.GetXaxis().SetTitle(
    fitter_mWW.ws.var(pars_mWW.var[0]).getTitle(True).Data()
    )
pull_mWW.GetYaxis().SetTitle("pull")
cp_mWW.Update()

parIter = params_mWW.createIterator()
p = parIter.Next()
while p:
    if not p.isConstant():
        p.setRange(p.getVal()-p.getError()*15, p.getVal()+p.getError()*15)
    # else:
    #     p.setRange(p.getVal()-p.getError()*3., p.getVal()+p.getError()*3.)
    p = parIter.Next()


c_bkg = TCanvas('c_bkg', 'histograms')
c_bkg.cd()

TH1F.SetDefaultSumw2(True)
bkgHisto = full_pdf.createHistogram("HWW%snujj_bkg" % mode,
                                    fitter_mWW.ws.var(pars_mWW.var[0]))
bkgHisto.Print()
bkgHisto.Scale(full_pdf.expectedEvents(fitter_mWW.ws.set('obsSet'))/bkgHisto.Integral())
bkgHisto.Print()

# for bi in range(1, bkgHisto.GetNbinsX()+1):
#     print bi, 'content:', bkgHisto.GetBinContent(bi), '+/-', bkgHisto.GetBinError(bi)
# print

bkgHisto.SetName("HWW%snujj_bkg" % mode)
# c_debug = TCanvas('c_debug', 'debug')
bkgHisto.Draw()
# if fr_mWW:
#     bkgHisto_up = fitter_mWW.utils.newEmptyHist(
#         'HWW%snujj_bkg_%sbkgshapeUp' % (mode, mode), 1)
#     bkgHisto_up = pulls.curveToHist(upper, bkgHisto_up)
#     bkgHisto_up.Print()
#     bkgHisto_up.Scale(full_pdf.expectedEvents(fitter_mWW.ws.set('obsSet'))/bkgHisto_up.Integral())
#     bkgHisto_up.SetLineColor(kOrange+2)
#     bkgHisto_up.SetLineStyle(kDashed)
#     bkgHisto_dwn = fitter_mWW.utils.newEmptyHist(
#         'HWW%snujj_bkg_%sbkgshapeDown' % (mode, mode), 1)
#     bkgHisto_dwn = pulls.curveToHist(lower, bkgHisto_dwn)
#     bkgHisto_dwn.Print()
#     bkgHisto_dwn.Scale(full_pdf.expectedEvents(fitter_mWW.ws.set('obsSet'))/bkgHisto_dwn.Integral())
#     bkgHisto_dwn.SetLineColor(kOrange+4)
#     bkgHisto_dwn.SetLineStyle(kDashed)
#     bkgHisto_up.Draw('same')
#     bkgHisto_dwn.Draw('same')

c_bkg.Update()
ggHPdf = fitter_mWW.ws.pdf('ggH')
ggHHisto = ggHPdf.createHistogram("HWW%snujj_ggH" % mode,
                                  fitter_mWW.ws.var(pars_mWW.var[0]))
ggHHisto.Scale(fitter_mWW.ws.var('n_ggH').getVal()/ggHHisto.Integral())
ggHHisto.SetName("HWW%snujj_ggH" % mode)
ggHHisto.SetLineColor(kBlue+2)
ggHHisto.Draw('same')

interf = []
if fitter_mWW.ws.pdf('ggH_interf_ggHDown'):
    for direction in ['Up', 'Down']:
        interfPdf = fitter_mWW.ws.pdf('ggH_interf_ggH%s' % direction)
        interf.append(interfPdf.createHistogram(
            "HWW%snujj_ggH_interf_%s" % (mode, direction),
            fitter_mWW.ws.var(pars_mWW.var[0])))
        interf[-1].Scale(
            fitter_mWW.ws.var('n_ggH_interf_ggH%s' % direction).getVal()/ \
            interf[-1].Integral())
        interf[-1].SetName('HWW%snujj_ggH_interf_ggH%s' % (mode,direction))
        interf[-1].SetLineColor(kBlue+2)
        interf[-1].SetLineStyle(kDashed)
        interf[-1].Draw('same')
        
                     
        
qqHPdf = fitter_mWW.ws.pdf('qqH')
qqHHisto = qqHPdf.createHistogram("HWW%snujj_qqH" % mode,
                                  fitter_mWW.ws.var(pars_mWW.var[0]))
qqHHisto.Scale(fitter_mWW.ws.var('n_qqH').getVal()/qqHHisto.Integral())
qqHHisto.SetName("HWW%snujj_qqH" % mode)
qqHHisto.SetLineColor(kRed+2)
qqHHisto.Draw('same')

dataHisto = RooAbsData.createHistogram(fitter_mWW.ws.data('data_obs'),
                                       'HWW%snujj_data_obs' % mode,
                                       fitter_mWW.ws.var(pars_mWW.var[0]))
dataHisto.SetMarkerStyle(20)
dataHisto.SetName('HWW%snujj_data_obs' % mode)
dataHisto.Draw('same')
dataHisto.Print()
# obs = RooArgList(fitter_mWW.ws.var(pars_mWW.var[0]))
# bkgf = full_pdf.asTF(obs,RooArgList(params_mWW),RooArgSet(obs))
# bkgf.SetLineColor(2)
# bkgf.SetLineWidth(3)
# bkgf.Draw('same')
c_bkg.Update()

print 'histogram chi2 test:'
pval = dataHisto.Chi2Test(bkgHisto, 'uwpchi2')
print 'chi2:',pval,'n fit params:',nFreeFitParams,
print 'ndf:',dataHisto.GetNbinsX()-nFreeFitParams-1,
print 'p-value:',TMath.Prob(pval, dataHisto.GetNbinsX()-nFreeFitParams-1)
# chi2f = dataHisto.Chisquare(bkgf)
# print 'chi2f:', chi2f
print

fitter_mWW.ws.saveSnapshot('nullFitSnapshot', fitter_mWW.ws.allVars())

fitter_mWW.ws.loadSnapshot('nullFitSnapshot')

if opts.obsLimit:
    limit = limits.plcLimit(fitter_mWW.ws.var(pars_mWW.var[0]),
                            fitter_mWW.ws.var('r_signal'),
                            full_pdf, fitter_mWW.ws,
                            fitter_mWW.ws.data('data_obs'))
    
    print '%.0f%% CL upper limit' % (95.), limit['r_signal']['ok'],
    print ': %.4f' % (limit['r_signal']['upper'])
    print '%.0f%% CL lower limit' % (95.), limit['r_signal']['ok'],
    print ': %.4f' % (limit['r_signal']['lower'])

fitter_mWW.ws.loadSnapshot('nullFitSnapshot')

#freeze all parameters in place
#params_mWW.setAttribAll('Constant', True)

import SignalShapeSystematic

allVars = fitter_mWW.ws.allVars()
allVars.remove(fitter_mWW.ws.set('obsSet'))

if opts.mH >= 400:
    print "\nsetting systematic kappa parameters"
    SignalShapeSystematic.SetSystematicAlphas(fitter_mWW.ws)
    print

    fitter_mWW.ws.var('interf_ggH').setConstant(False)

    allVars.remove(fitter_mWW.ws.var('interf_ggH'))
    allVars.remove(fitter_mWW.ws.var('interf_ggH_interf_ggHUp'))
    allVars.remove(fitter_mWW.ws.var('interf_ggH_interf_ggHDown'))
    
varIter = allVars.createIterator()
theVar = varIter.Next()
while theVar:   
    if opts.isElectron:
        theVar.SetName('%s_el' % theVar.GetName())
    else:
        theVar.SetName('%s_mu' % theVar.GetName())
    theVar = varIter.Next()

extraTag = ''
if opts.includeSignal:
    extraTag = '_withSignal'
if mvaCutOverride:
    extraTag += '_mvaCut' + str(mvaCutOverride)
output = TFile("HWW%ilnujj_%s_%ijets_1D2Fit%s_output.root" % \
                   (opts.mH,mode, opts.Nj, extraTag),
               "recreate")

if fr:
    fr.SetName('fit_result')
    fr.Write()
if fr_mWW:
    fr_mWW.SetName('fit_result_mWW')
    fr_mWW.Write()

if not opts.toy:
    plot1.Write()
    pull1.Write()
plot_mWW.Write()
pull_mWW.Write()
plot_mWW_withErrs.Write()
bkgHisto.Write()
# if fr_mWW:
#     bkgHisto_up.Write()
#     bkgHisto_dwn.Write()
qqHHisto.Write()
ggHHisto.Write()
for histo in interf:
    histo.Write()
dataHisto.Write()
fitter.ws.Write()
# fitter_low.ws.Write()
fitter_mWW.ws.Write()
#likPlot.Write()
if upperHist:
    upperHist.Write()

#fitter.ws.Print()
output.Close()

print 'Time elapsed: %.1f sec' % timer.RealTime()
print 'CPU time used: %.1f sec' % timer.CpuTime()

if not opts.toy:
    print '***',pars.var[0], 'fit ***'
    print '%i degrees of freedom' % ndf
    print 'chi2: %.2f / %i = %.2f' % (chi2, ndf, (chi2/ndf))
    print 'chi2 probability: %.4g' % (TMath.Prob(chi2, ndf))

# print '*** low sideband fit ***'
# print '%i degrees of freedom' % ndf_low
# print 'chi2: %.2f / %i = %.2f' % (chi2_low, ndf_low, 
#                                   (chi2_low/ndf_low))
# print 'chi2 probability: %.4g' % (TMath.Prob(chi2_low, ndf_low))

# print '*** high sideband fit ***'
# print '%i degrees of freedom' % ndf_high
# print 'chi2: %.2f / %i = %.2f' % (chi2_high, ndf_high, 
#                                   (chi2_high/ndf_high))
# print 'chi2 probability: %.4g' % (TMath.Prob(chi2_high, ndf_high))

print '*** mWW fit ***'
if fr_mWW:
    print 'minimum NLL:', fr_mWW.minNll()
print '%i degrees of freedom' % ndf_mWW
print 'chi2: %.2f / %i = %.2f' % (chi2_mWW, ndf_mWW, 
                                  (chi2_mWW/ndf_mWW))
print 'chi2 probability: %.4g' % (TMath.Prob(chi2_mWW, ndf_mWW))
