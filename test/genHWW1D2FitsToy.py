from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False,
                  help='no X11 windows')
parser.add_option('-m', '--mode', default="HWW1D2FitsConfig", 
                  dest='modeConfig',
                  help='which config to select look at HWW2DConfig.py for ' +\
                      'an example.  Use the file name minus the .py extension.'
                  )
parser.add_option('--mode_mWW', default="HWW1D2FitsConfig_mWW", 
                  dest='modeConfig_mWW',
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
parser.add_option('--toy', dest='toy', action='store_false', default=True,
                  help='use pseudo-data instead of data file')
parser.add_option('--seed', dest='seed', type='int', help='random seed')
parser.add_option('--debug', dest='debug', action='store_true', default=False,
                  help='turn on extra debugging information')
parser.add_option('--mva', dest='mvaCut', type='float',
                  help='override cut value for mva')
parser.add_option('--xrootd', dest='xrootd', action='store_true',
                  help='use xrootd file opening.')
parser.add_option('--injectS', type='float', dest='sigInject',
                  help='amount of signal to inject')

(opts, args) = parser.parse_args()

import pyroot_logon

#import HWW2DConfig
config = __import__(opts.modeConfig)
import RooWjj2DFitter

from ROOT import TCanvas, RooFit, RooLinkedListIter, TMath, RooRandom, TFile, \
    RooDataHist, RooMsgService, TStopwatch, RooAbsPdf, RooAbsData, \
    RooWorkspace, RooArgList, RooAddPdf
import pulls

timer = TStopwatch()
timer.Start()

#RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9)
#RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9)
if not opts.debug:
    RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)
    RooMsgService.instance().addStream(RooFit.ERROR,
                                       RooFit.Prefix(True), 
                                       RooFit.ClassName('RooExpPoly'),
                                       RooFit.OutputFile('/dev/null'))
    RooMsgService.instance().Print('v')

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

mode = 'muon'
if opts.isElectron:
    mode = 'electron'

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

fitter.readParametersFromFile()
fitter.expectedFromPars()
fitter.resetYields()

startpars = totalPdf.getParameters(fitter.ws.set('obsSet'))
fitter.ws.defineSet("params", startpars)
fitter.ws.saveSnapshot("initPars", startpars)

components = list(pars.backgrounds) + list(pars.signals)
var_name = pars.var[0]
combinedWS = RooWorkspace('w_mjj', 'combined ws')
for comp in components:
    compPdf = fitter.ws.pdf(comp)
    getattr(combinedWS, 'import')(compPdf)
    newPdf = combinedWS.pdf(comp)
    newPdf.SetName('%s_%s' % (comp, var_name))
    yieldf = fitter.ws.function('f_%s_norm' % comp)
    getattr(combinedWS, 'import')(yieldf)

fitter.ws.Print()

# if opts.toy:
#     #generate toy dataset
#     print 'Generated parameters'
#     fitter.ws.set('params').Print('v')
#     fitter.ws.saveSnapshot("genPars", startpars)

#     data = totalPdf.generate(fitter.ws.set('obsSet'), RooFit.Name('data_obs'),
#                              RooFit.Extended())
#     if fitter.pars.binData:
#         data.SetName('data_unbinned')
#         getattr(fitter.ws, 'import')(data)
#         data = RooDataHist('data_obs', 'data_obs', fitter.ws.set('obsSet'),
#                            data)
#     data.Print('v')
#     getattr(fitter.ws, 'import')(data)
# else:    
#     data = fitter.loadData()

# data.Print()
startpars.IsA().Destructor(startpars)

print 'Time elapsed: %.1f sec' % timer.RealTime()
print 'CPU time used: %.1f sec' % timer.CpuTime()
timer.Continue()

sigInt = fitter.ws.pdf('WpJ').createIntegral(fitter.ws.set('obsSet'),
                                             fitter.ws.set('obsSet'),
                                             'signalRegion')
fullInt = fitter.ws.pdf('WpJ').createIntegral(fitter.ws.set('obsSet'),
                                              fitter.ws.set('obsSet'))
n_WpJ_sig = sigInt.getVal()/fullInt.getVal()*fitter.ws.function('f_WpJ_norm').getVal()
print 'W+jets integral in signal region:',sigInt.getVal()/fullInt.getVal(),
print 'W+jets yield in the signal region:',n_WpJ_sig

config_mWW = __import__(opts.modeConfig_mWW)
print mWWArgs
pars_mWW = config_mWW.theConfig(Nj = opts.Nj, mH = opts.mH, 
                                isElectron = opts.isElectron, 
                                initFile = mWWArgs,
                                includeSignal = True,
                                MVACutOverride = mvaCutOverride,
                                mjj_config = opts.modeConfig,
                                xrootd = opts.xrootd)
pars_mWW.yieldConstraints['WpJ'] = fitter.ws.var('WpJ_nrm').getError()

if opts.toy:
    pars_mWW.blind = False

# add systematic errors multipliers to ggH signals
if (opts.mH >= 400):
    pars_mWW.ggHdoSystMult = True

fitter_mWW = RooWjj2DFitter.Wjj2DFitter(pars_mWW)
fitter_mWW.ws.SetName("w_mWW")
totalPdf_mWW = fitter_mWW.makeFitter()
fitter_mWW.readParametersFromFile()
# WpJ_mWW = fitter_mWW.ws.var('n_WpJ')
# WpJ_mWW.setVal(n_WpJ_sig)
fitter_mWW.expectedFromPars()

mWWCut = '((%s>%.0f)&&(%s<%.0f))' % \
    (pars.var[0], pars.exclude[pars.var[0]][0],
     pars.var[0], pars.exclude[pars.var[0]][1])
print 'signal region cut:',mWWCut

fitter_mWW.resetYields()
params_mWW = totalPdf_mWW.getParameters(fitter_mWW.ws.set('obsSet'))

predictedPars = params_mWW.snapshot()

if opts.sigInject:
    fitter_mWW.ws.var('r_signal').setVal(opts.sigInject)
fitter_mWW.ws.var('r_signal').setError(0.1)
fitter_mWW.ws.var('r_signal').setRange(-3., 9.)
fitter_mWW.ws.var('r_signal').setConstant(False)

params_mWW.Print("v")
fitter_mWW.ws.defineSet("params", params_mWW)

var_name = pars_mWW.var[0]
other_var = pars_mWW.var[1]
compPdfs = []
for comp in components:
    compPdf = fitter_mWW.ws.pdf(comp)
    getattr(combinedWS, 'import')(compPdf)
    newPdf = combinedWS.pdf(comp)
    newPdf.SetName('%s_%s' % (comp, var_name))
    norm = combinedWS.function('f_%s_norm' % comp)
    combinedWS.factory('PROD::%s(%s_%s, %s)' % (comp, comp, other_var,
                                                newPdf.GetName()))
    compPdfs.append(combinedWS.factory('RooExtendPdf::%s_extended(%s, %s)' %\
                                           (comp, comp, norm.GetName())))

combinedWS.defineSet('obsSet', '%s,%s' % (other_var, var_name))

if opts.sigInject:
    combinedWS.var('r_signal').setVal(opts.sigInject)
combinedWS.var('r_signal').setError(0.1)
combinedWS.var('r_signal').setRange(-3., 9.)
combinedWS.var('r_signal').setConstant(False)

compNames = [ c.GetName() for c in compPdfs ]
compList = RooArgList(combinedWS.argSet(','.join(compNames)))
getattr(combinedWS, 'import')(RooAddPdf('total', 'total', compList))
combinedPdf = combinedWS.pdf('total')
genPars = combinedPdf.getParameters(combinedWS.set('obsSet'))
combinedWS.defineSet('params', genPars)

if opts.toy:
    #generate toy dataset
    print 'Generated parameters'
    combinedWS.set('params').Print('v')
    combinedWS.saveSnapshot("genPars", combinedWS.set('params'))

    data = combinedPdf.generate(combinedWS.set('obsSet'), 
                                RooFit.Name('data_obs'),
                                RooFit.Extended())
    data.Print('v')
    getattr(combinedWS, 'import')(data)

extraTag = ''
if opts.includeSignal:
    extraTag = '_withSignal'
if mvaCutOverride:
    extraTag += '_mvaCut' + str(mvaCutOverride)

combinedWS.Print()
output = TFile("HWW%ilnujj_%s_%ijets_1D2Fit%s_gen.root" % \
                   (opts.mH,mode, opts.Nj, extraTag),
               "recreate")

# fitter.ws.Write()
# fitter_mWW.ws.Write()

combinedWS.Write()

output.Close()

print 'Time elapsed: %.1f sec' % timer.RealTime()
print 'CPU time used: %.1f sec' % timer.CpuTime()

