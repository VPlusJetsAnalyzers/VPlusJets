from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False,
                  help='no X11 windows')
parser.add_option('-m', '--mode', default="HWW2DConfig", dest='modeConfig',
                  help='which config to select look at HWW2DConfig.py for ' +\
                      'an example.  Use the file name minus the .py extension.'
                  )
parser.add_option('--mjj', dest='mjj_config',
                  help='which config to select for the mjj parameters' +\
                      'an example.  Use the file name minus the .py extension.'
                  )
parser.add_option('-H', '--mH', dest='mH', default=350, type='int',
                  help='Higgs Mass Point')
parser.add_option('-j', '--Njets', dest='Nj', default=2, type='int',
                  help='Number of jets.')
parser.add_option('--bn', dest='bn', default='ShapeParameters',
                  help='basis name for the output files.')
parser.add_option('--comp', dest='component', default='diboson',
                  help='name of component to fit')
parser.add_option('--electrons', dest='isElectron', action='store_true',
                  default=False, help='do electrons instead of muons')
parser.add_option('--makeFree', dest='makeConstant', action='store_false',
                  default=True, help='make parameters free in output')
parser.add_option('--sideband', dest='sb', type='int',
                  default=0, help='use sideband dataset and model instead')
parser.add_option('--signal', dest='sig', action='store_true',
                  default=False, help='use signal model instead')
parser.add_option('--alternateModel', dest='altModel', action='store_true',
                  default=False, help='use the alternate model instead')
parser.add_option('-i', '--interference', dest='interference', default=0, 
                  type='int', help='ggH interference to use.  '+\
                      '(0): none [default]  (1): nominal'+\
                      '  (2): +1 sigma  (3): -1 sigma  ')
parser.add_option('--morphingComp', dest='morphComponent', default=0, 
                  type='int', help='which morphing component to use.  '+\
                      '(0): nominal [default]  (+/- 1): matching up/down'+\
                      '  (+/- 2): scale up/down')
parser.add_option('--ws', dest='ws', help='filename to get data from instead' +\
                  ' reading from source ntuples again.')
parser.add_option('--btag', dest='btag', action='store_true',
                  default=False, help='Use b-tagged selection.')
parser.add_option('--mva', dest='mvaCut', type='float',
                  help='override cut value for mva')
parser.add_option('--rateSyst', dest='RateSyst', type='int', 
                  help='do a rate systematic in region of +/- 10% of mH')
parser.add_option('--xrootd', dest='xrootd', action='store_true',
                  help='use xrootd file opening.')
parser.add_option('--nativeBins', dest='nativeBins', action='store_true',
                  help='don\'t use plot bins.')
parser.add_option('--overrideModel', dest='overrideModel', type='int',
                  help='override the model number')
parser.add_option('--overrideAux', dest='overrideAux', type='int',
                  help='override the aux model number')
parser.add_option('--fixToZero', dest='fixToZero', action='store_true',
                  help='fixLowParametersToZero')
parser.add_option('--Cprime', dest='Cprime', type='float',
                  help='Cprime value')
parser.add_option('--BRnew', dest='BRnew', type='float',
                  help='BRnew value')
parser.add_option('--refit', dest='refit', action='store_true',
                  help='refit parameters')

(opts, args) = parser.parse_args()


import pyroot_logon

#import HWW2DConfig
config = __import__(opts.modeConfig)
#import VBF2DConfig
import RooWjj2DFitter
import HWWSignalShapes

from ROOT import RooFit, TCanvas, RooArgSet, TFile, RooAbsReal, RooAbsData, \
    RooHist, TMath, kRed, kDashed, kOrange, RooMsgService, RooDataSet, \
    RooCmdArg, TH1F

from array import array

import pulls

RooMsgService.instance().setGlobalKillBelow(RooFit.ERROR)

configArgs = {'Nj': opts.Nj, 'mH': opts.mH, 'isElectron': opts.isElectron,
              'initFile' : args}
mvaCutOverride = None
if hasattr(opts, "mvaCut") and (opts.mvaCut != None):
    mvaCutOverride = opts.mvaCut
    configArgs['MVACutOverride'] = opts.mvaCut

if opts.mjj_config:
    configArgs['mjj_config'] = opts.mjj_config

sb = None
if opts.sb > 0:
    sb = 'high'
    configArgs['sideband'] = sb
elif opts.sb < 0:
    sb = 'low'
    configArgs['sideband'] = sb

# pars = config.theConfig(Nj = opts.Nj, mH = opts.mH, 
#                         isElectron = opts.isElectron, initFile = args)
if opts.btag:
    # pars = config.theConfig(Nj = opts.Nj, mH = opts.mH, 
    #                         isElectron = opts.isElectron, initFile = args,
    #                         btagged = opts.btag)
    configArgs['btagged'] = opts.btag
# else:
#     pars = config.theConfig(Nj = opts.Nj, mH = opts.mH, 
#                             isElectron = opts.isElectron, initFile = args)
# if mvaCutOverride:
#     pars = config.theConfig(Nj = opts.Nj, mH = opts.mH, 
#                             isElectron = opts.isElectron, initFile = args,
#                             MVACutOverride = mvaCutOverride)
if opts.xrootd:
    configArgs['xrootd'] = opts.xrootd

pars = config.theConfig(**configArgs)


files = getattr(pars, '%sFiles' % opts.component)
models = getattr(pars, '%sModels' % opts.component)
if hasattr(pars, '%sConvModels' % opts.component):
    convModels  = getattr(pars, '%sConvModels' % opts.component)
else:
    convModels = None

if opts.altModel:
    print 'will fit the alternate model'
    models = getattr(pars, '%sModelsAlt' % opts.component)
    convModels  = getattr(pars, '%sConvModelsAlt' % opts.component)
# if opts.sb:
#     pars.cuts = pars.SidebandCuts

if opts.overrideModel:
    print 'overriding models', models
    models = [ opts.overrideModel ]
if opts.overrideAux:
    print 'overriding aux model'
    setattr(pars, '%sAuxModels' % opts.component, [opts.overrideAux])

compName = opts.component
morphingPdf = False
if models[0] == -2:
    print 'pdf morphing'
    morphingPdf = True
    if opts.morphComponent == 0:
        models = getattr(pars, '%s_NomModels' % opts.component)
        files = getattr(pars, '%s_NomFiles' % opts.component)
        compName = opts.component + '_Nom'
    elif opts.morphComponent == -1:
        models = getattr(pars, '%s_MDModels' % opts.component)
        files = getattr(pars, '%s_MDFiles' % opts.component)
        compName = opts.component + '_MD'
    elif opts.morphComponent == 1:
        models = getattr(pars, '%s_MUModels' % opts.component)
        files = getattr(pars, '%s_MUFiles' % opts.component)
        compName = opts.component + '_MU'
    elif opts.morphComponent == -2:
        models = getattr(pars, '%s_SDModels' % opts.component)
        files = getattr(pars, '%s_SDFiles' % opts.component)
        compName = opts.component + '_SD'
    elif opts.morphComponent == 2:
        models = getattr(pars, '%s_SUModels' % opts.component)
        files = getattr(pars, '%s_SUFiles' % opts.component)
        compName = opts.component + '_SU'
print 'component', compName
print 'models', models
print 'files', files

# print 'eff lumi WW: %.1f invpb' % (pars.dibosonFiles[0][1]/pars.dibosonFiles[0][2])
# print 'eff lumi WZ: %.1f invpb' % (pars.dibosonFiles[1][1]/pars.dibosonFiles[1][2])

#pars.DataFile = files[0][0]
fitter = RooWjj2DFitter.Wjj2DFitter(pars)

#fitter.makeFitter()
#data = fitter.loadData(True)
#data.Print()
data = None
sumNExp = 0.
sumxsec = 0.
weighted = True
cutOverride = None
if compName == 'QCD':
    weighted = False
    cutOverride = pars.QCDcuts
cpw = False
if (compName == "ggH") or (compName == 'qqH'):
    cpw = True
if (compName == 'qqH') and (opts.mH == 170):
    cpw = False
#print 'interference option:',opts.interference
refEffLumi = -1.
print 'integrated lumi:',pars.integratedLumi
in_ws = None
loadIt = True
if opts.ws:
    in_ws = TFile(opts.ws).Get('wjj2dfitter')
    if in_ws.data('data'):
        getattr(fitter.ws, 'import')(in_ws.data('data'))
        data = fitter.ws.data('data')
        loadIt = False
for (ifile, (filename, ngen, xsec)) in enumerate(files):
    if refEffLumi > 0.:
        scale = refEffLumi/(ngen/xsec)
    else:
        scale = 1.0
    if in_ws and in_ws.data('data%i' % ifile):
        getattr(fitter.ws, 'import')(in_ws.data('data%i' % ifile))
    dataArgs  = {'fnames': filename, 
                 'dsName': 'data%i' % ifile, 
                 'ws': fitter.ws,
                 'weighted': weighted, 
                 'CPweight': cpw,
                 'cutOverride': cutOverride,
                 'interference': opts.interference,
                 'additionalWgt': scale}
    if opts.Cprime != None:
        dataArgs['Cprime'] = opts.Cprime
    if opts.BRnew != None:
        dataArgs['BRnew'] = opts.BRnew
    tmpData = fitter.utils.File2Dataset(**dataArgs)

    tmpData.Print()
    print filename,'effective integrated lumi:',ngen/xsec
    expectedYield = xsec*pars.integratedLumi*tmpData.sumEntries()/scale/ngen
    print filename,'A x eff: %.3g' % (tmpData.sumEntries()/scale/ngen)
    print filename,'expected yield: %.1f' % expectedYield
    print filename,'scale:',scale
    sumNExp += expectedYield
    sumxsec += xsec
    if refEffLumi < 0:
        refEffLumi = ngen/xsec
    if loadIt:
        if not data:
            data = tmpData.Clone('data')
        else:
            data.append(tmpData)

print compName,'total expected yield: %.1f' % sumNExp
print compName,'overall A x eff: %.3g' % (sumNExp/sumxsec/pars.integratedLumi)
data.Print()
if not fitter.ws.data('data'):
    getattr(fitter.ws, 'import')(data)

hist2d = None
try:
    obs = [ pars.varNames[x] for x in pars.var ]
except AttributeError:
    obs = pars.var

if (len(pars.var) > 1):
    hist2d = data.createHistogram(fitter.ws.var(obs[0]),
                                  fitter.ws.var(obs[1]), '',
                                  'hist2d')
    print 'correlation between the first two observables:', hist2d.GetCorrelationFactor()

if opts.interference in [1,2,3]:
    fitter.makeFitter()
    # fitter.ws.Print()
    if opts.interference == 1:
        pdfName = '%s' % compName
    elif opts.interference == 2:
        pdfName = '%s_interf_%sUp' % (compName, compName)
    elif opts.interference == 3:
        pdfName = '%s_interf_%sDown' % (compName, compName)
    print pdfName
    sigPdf = fitter.ws.pdf(pdfName)
    if sigPdf:
        sigPdf.Print()
    else:
        print 'failed to fined pdf',pdfName
else:
    sigPdf = fitter.makeComponentPdf(compName, files, models, opts.altModel, convModels)
    if sb:
        sbPdf = fitter.ws.pdf('%s_%s_%s' % (compName,pars.var[0],sb))
        if sbPdf:
            sigPdf = sbPdf
            sigPdf.getParameters(data).setAttribAll('Constant', False)
        sigPdf.Print()
    elif opts.sig and fitter.ws.pdf('%s_%s_sig' % (compName,pars.var[0])):
        sigPdf = fitter.ws.pdf('%s_%s_sig' % (compName,pars.var[0]))


extraTag = ''
if opts.interference == 2:
    extraTag = '_interf_%sUp' % compName
if opts.interference == 3:
    extraTag = '_interf_%sDown' % compName

if fitter.ws.var('mean_%s_fit_mlvjj_core%s' % (compName,extraTag)):
    fitter.ws.var('mean_%s_fit_mlvjj_core%s' % \
                      (compName,extraTag)).setVal(opts.mH)
if fitter.ws.var('mean_%s_fit_mlvjj%s' % (compName,extraTag)):
    fitter.ws.var('mean_%s_fit_mlvjj%s' % \
                      (compName,extraTag)).setVal(opts.mH)
if fitter.ws.var('mean_%s_fit_mlvjj_tail%s' % (compName,extraTag)):
    fitter.ws.var('mean_%s_fit_mlvjj_tail%s' % \
                      (compName,extraTag)).setVal(opts.mH)
if fitter.ws.var('sigma_%s_fit_mlvjj_tail%s' % (compName,extraTag)):
    fitter.ws.var('sigma_%s_fit_mlvjj_tail%s' % \
                      (compName,extraTag)).setVal(opts.mH*0.25)
if fitter.ws.var('sigma_%s_fit_mlvjj_core%s' % (compName,extraTag)):
    fitter.ws.var('sigma_%s_fit_mlvjj_core%s' % \
                      (compName,extraTag)).setVal(opts.mH*0.1)
if (compName in ['ggH','qqH']) and \
        fitter.ws.var('width_%s_fit_mlvjj%s' % (compName,extraTag)):
    fitter.ws.var('width_%s_fit_mlvjj%s' % \
                      (compName,extraTag)).setVal(HWWSignalShapes.HiggsWidth[int(opts.mH)])

params = sigPdf.getParameters(data)
parCopy = params.snapshot()
for filename in args:
    parCopy.readFromFile(filename)

if opts.refit or opts.makeConstant:
    params.assignValueOnly(parCopy)
else:
    params.__assign__(parCopy)

params.Print('v')
parCopy.IsA().Destructor(parCopy)
    
if fitter.ws.var('npow_ggH_Mass2j_PFCor'):
    # fitter.ws.var('npow_ggH_Mass2j_PFCor').setConstant(False)
    # fitter.ws.var('alpha_ggH_Mass2j_PFCor').setVal(1.0)
    # fitter.ws.var('alpha_ggH_Mass2j_PFCor').setConstant(True)
    pass

if fitter.ws.var('npow_diboson_Mass2j_PFCor'):
    fitter.ws.var('npow_diboson_Mass2j_PFCor').setConstant(True)
    fitter.ws.var('npow_diboson_Mass2j_PFCor').setError(0.5)

if fitter.ws.var('npow_WpJ_Mass2j_PFCor'):
    fitter.ws.var('npow_WpJ_Mass2j_PFCor').setConstant(False)

fitter.ws.Print()


fr = None
if models[0] >= 0:
    fr = sigPdf.fitTo(data, RooFit.Save(), 
                      RooFit.SumW2Error(False),
                      # RooFit.Hesse(False),
                      RooFit.Minimizer("Minuit2", "minimize"),
                      # RooFit.Minos(True),
                      # RooFit.InitialHesse(True)
                      )


mode = 'muon'
if opts.isElectron:
    mode = 'electron'

cans = []
plots = []
chi2s = []
ndfs = []

dhists = []
fhists = []
hchi2s = []

TH1F.SetDefaultSumw2(True)

for (i,m) in enumerate(models):
    par = obs[i]
    c1 = TCanvas('c%i' % i, par)
    binArg = RooFit.Bins(fitter.ws.var(par).getBins('plotBins'))
    if opts.nativeBins:
        binArg = RooCmdArg.none()
    sigPlot = fitter.ws.var(par).frame(RooFit.Name('%s_Plot' % par),
                                       RooFit.Range('plotRange'),
                                       binArg
                                       )
    dataHist = RooAbsData.createHistogram(data,'dataHist_%s' % par,
                                          fitter.ws.var(par))
    fitHist = sigPdf.createHistogram('sigHist_%s' % par,
                                     fitter.ws.var(par))
    fitHist.Scale(dataHist.Integral()/fitHist.Integral())
    #this will let it ignore empty bins
    for bin in range(1, dataHist.GetNbinsX()+1):
        if dataHist.GetBinContent(bin) < 1e-9:
            dataHist.SetBinContent(bin, fitHist.GetBinContent(bin))
            dataHist.SetBinError(bin, 1.0)
    h_chi2 = dataHist.Chi2Test(fitHist, 'wwpchi2')
    fhists.append(fitHist)
    dhists.append(dataHist)
    hchi2s.append(h_chi2)
    # theData = RooHist(dataHist, 1., 1, RooAbsData.SumW2, 1.0, True)
    # theData.SetName('theData')
    # theData.SetTitle('data')
    # sigPlot.addPlotable(theData, 'pe', False, True)
    data.plotOn(sigPlot, RooFit.Name('theData'),
                RooFit.DataError(RooAbsData.SumW2))
    sigPlot.getHist('theData').SetTitle('data')
    if fr:
        sigPdf.plotOn(sigPlot, 
                      RooFit.VisualizeError(fr, 1, True),
                      RooFit.Name('fitErrors'),
                      RooFit.FillColor(kOrange+1),
                      RooFit.FillStyle(3001))
        sigPlot.getCurve().SetTitle('Fit errors')
    sigPdf.plotOn(sigPlot, RooFit.Name('fitCurve'))
    sigPlot.getCurve().SetTitle('%s fit' % compName)
    ret = sigPdf.plotOn(sigPlot, RooFit.Name('fitTails'),
                        RooFit.Components('*tail'),
                        RooFit.LineColor(kRed),
                        RooFit.LineStyle(kDashed))
    if ret:
        sigPlot.getCurve().SetTitle('tail')

    # sigPlot.GetYaxis().SetTitle('Events / GeV')

    sigPlot.getHist('theData').SetTitle('%s MC' % compName)

    leg = RooWjj2DFitter.Wjj2DFitter.legend4Plot(sigPlot)
    sigPlot.addObject(leg)
    sigPlot.Draw()
    c1.Update()
    sigPlot.SetMinimum(1e-6)
    cans.append(c1)
    plots.append(sigPlot)
    (chi2_1, ndf_1) = pulls.computeChi2(sigPlot.getHist('theData'),
                                        sigPlot.getCurve('fitCurve'))
    chi2s.append(chi2_1)
    ndfs.append(ndf_1)

    residuals = pulls.createResid(sigPlot.getHist('theData'),
                                  sigPlot.getCurve('fitCurve'))
    plots.append(residuals)
    residuals.SetName('%s_Residuals' % par)
    c3 = TCanvas('c%i_residual' % i, par + ' residual')
    c3.SetGridy()
    residuals.Draw('ap')
    residuals.GetXaxis().SetLimits(fitter.ws.var(par).getMin(),
                                   fitter.ws.var(par).getMax())
    c3.Update()
    cans.append(c3)

    pull = pulls.createPull(sigPlot.getHist('theData'),
                            sigPlot.getCurve('fitCurve'))
    plots.append(pull)
    pull.SetName('%s_Pulls' % par)
    c2 = TCanvas('c%i_pull' % i, par + ' pull')
    c2.SetGridy()
    pull.Draw('ap')
    pull.GetXaxis().SetLimits(fitter.ws.var(par).getMin(),
                              fitter.ws.var(par).getMax())
    c2.Update()
    cans.append(c2)

    c1.Print('%s_compFit.png' % opts.bn)
#    c1.Print('%s_compFit.pdf' % opts.bn)

mode = 'muon'
if opts.isElectron:
    mode = 'electron'

ndf = 0

finalPars = params.snapshot()

if fr and (opts.RateSyst):
    var = fitter.ws.var(obs[0])
    var.setRange("SystSignal", opts.mH*0.9, opts.mH*1.1)
    sigPdfInt = sigPdf.createIntegral(fitter.ws.set('obsSet'),
                                      fitter.ws.set('obsSet'))
    sigPdfInt_sig = sigPdf.createIntegral(fitter.ws.set('obsSet'),
                                          fitter.ws.set('obsSet'),
                                          'SystSignal')
    fracs = [sigPdfInt_sig.getVal()/sigPdfInt.getVal()]
    while len(fracs) < opts.RateSyst:
        params.assignValueOnly(fr.randomizePars())
        frac = sigPdfInt_sig.getVal()/sigPdfInt.getVal()
        if (frac > 0) and (frac < 1):
            fracs.append(sigPdfInt_sig.getVal()/sigPdfInt.getVal())

    farray = array('d', fracs)
    print 'yield in signal region around %i: %.4f +/- %.4f' % \
        (opts.mH, TMath.Mean(len(farray),farray), 
         TMath.RMS(len(farray), farray))

sigFile = TFile('%s.root' % opts.bn, 'recreate')
if fr:
    fr.Print('v')
    fr.Write('fr')
for plot in plots:
    plot.Write()
fitter.ws.Write()
if hist2d:
    hist2d.Write()

sigFile.Close()

parIter = finalPars.createIterator()
p = parIter.Next()
fixCnt = 0
while p:
    if not p.isConstant():
        p.setRange(max(p.getVal()-p.getError()*20., p.getMin()),
                   min(p.getVal()+p.getError()*20., p.getMax()))
        if (p.getError() > abs(p.getVal()*10)):
            fixCnt += 1
            print p.GetName(), 'is not significantly deviated from zero',
            print '(%f' % (p.getVal()/p.getError()), 'sigma),',
            print 'consider fixing it to zero.'
            if opts.fixToZero:
                print 'fixing',p.GetName(),'to zero'
                p.setVal(0.0)
                p.setConstant()
        # elif (fixCnt > 1) and opts.fixToZero:
        #     print 'fixing',p.GetName(),'to zero'
        #     p.setVal(0.0)
        #     p.setConstant()
        else:
            fixCnt = 0
            
    p = parIter.Next()

if opts.makeConstant:
    finalPars.setAttribAll('Constant', True)
finalPars.writeToFile("%s.txt" % opts.bn)

if compName != 'QCD':
    paramsFile = open('%s.txt' % opts.bn)
    lines = paramsFile.readlines()
    paramsFile.close()
    outfile = open('%s.txt' % opts.bn, 'w')
    for line in lines:
        tag = ' '
        if opts.sb == 2:
            tag = '_low '
        elif opts.sb == 3:
            tag = '_high '
        words = line.split()
        outfile.write(words[0] + tag + ' '.join(words[1:]) + '\n')
    yieldName = compName
    if opts.interference == 2:
        yieldName += '_interf_%sUp' % compName
    elif opts.interference == 3:
        yieldName += '_interf_%sDown' % compName
    # if opts.sb:
    #     yieldName = 'dummy'
    if morphingPdf and (opts.morphComponent == 0):
        yieldName = opts.component
        outfile.write('fMU_%s = 0.0 +/- 10. L(-1.0 - 1.0)\n' % (yieldName))
        outfile.write('fSU_%s = 0.0 +/- 10. L(-1.0 - 1.0)\n' % (yieldName))
    outfile.write('n_%s = %.1f +/- %.1f C\n' % (yieldName, sumNExp, 
                                                TMath.Sqrt(sumNExp))
                  )
    outfile.close()

ndf = 1
if fr:
    ndf = fr.floatParsFinal().getSize()
finalPars.IsA().Destructor(finalPars)

params.IsA().Destructor(params)

print '%i free parameters in the fit' % ndf

for (i,chi2) in enumerate(hchi2s):
    print 'histogram chi2:',chi2,'ndf:',dhists[i].GetNbinsX()-ndf-1,
    print 'p-value:', TMath.Prob(chi2, dhists[i].GetNbinsX()-ndf-1)

chi2 = 0
print 'chi2: (',
for c in chi2s:
    chi2 += c
    if c != chi2s[-1]:
        print '%.2f +' % c,
    else:
        print '%.2f )' % c,
    ndf += 1

bins = 0
for b in ndfs:
    bins += b
ndf = bins-ndf

print '/%i = %.2f' % (ndf, (chi2/ndf))
print 'chi2 probability: %.4g' % (TMath.Prob(chi2, ndf))
