from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False,
                  help='no X11 windows')
parser.add_option('-m', '--mode', default="HWW2DConfig", dest='modeConfig',
                  help='which config to select look at HWW2DConfig.py for ' +\
                      'an example.  Use the file name minus the .py extension.'
                  )
parser.add_option('-H', '--mH', dest='mH', default=350, type='int',
                  help='Higgs Mass Point')
parser.add_option('-j', '--Njets', dest='Nj', default=2, type='int',
                  help='Number of jets.')
parser.add_option('--comp', dest='component', default='diboson',
                  help='name of component to fit')
parser.add_option('--electrons', dest='isElectron', action='store_true',
                  default=False, help='do electrons instead of muons')
parser.add_option('--nosig', dest='includeSignal', action='store_false',
                  default=True, help='include signal shape in fit')
parser.add_option('--fitparfn', dest='fitparfn', default='', type='string', help='the file name to be used to store fit output parameters')
parser.add_option('--useAlternateModelsToFit', dest='useAlternateModelsToFit', action='store_true', default=False, help='use alternate models to fit the data')
parser.add_option('--toy', dest='toy', action='store_true',
                  help='use pseudo-data instead of data file')
parser.add_option('--toyOut', dest='toyOut', help='filename for toy output')
parser.add_option('--runMCGenToySim', dest='runMCGenToySim', action='store_true', default=False, help='generate toys using MC and fit the combined distribution')
parser.add_option('--runPdfGenToySim', dest='runPdfGenToySim', action='store_true', default=False, help='generate toys using pdfs and fit the combined distribution')
parser.add_option('--useAlternateModelsToGen', dest='useAlternateModelsToGen', action='store_true', default=False, help='use alternate models when generating toy pdfs (and fit with the standard ones)')
parser.add_option('--genParamFile', dest='genParamFile', default='NONE', type='string', help='the parameters files to be used for generation')
parser.add_option('--extdiboson', dest='extdiboson', default=-1, type='int',
                  help='Set the number of diboson events to generate for toy validation')
parser.add_option('--extendedGen', dest='extendedGen', action='store_true',
                  default=False, help='use the extended option when generating')
parser.add_option('--extWpJ', dest='extWpJ', default=-1, type='int',
                  help='Set the number of WpJ events to generate for toy validation')
parser.add_option('--extWHbb', dest='extWHbb', default=-1, type='int',
                  help='Set the number of WHbb events to generate for toy validation')
parser.add_option('--exttop', dest='exttop', default=-1, type='int',
                  help='Set the number of top events to generate for toy validation')
parser.add_option('--extZpJ', dest='extZpJ', default=-1, type='int',
                  help='Set the number of ZpJ events to generate for toy validation')
parser.add_option('--extQCD', dest='extQCD', default=-1, type='int',
                  help='Set the number of QCD events to generate for toy validation')
parser.add_option('--seed', dest='seed', type='int', help='random seed')
parser.add_option('--ws', dest='ws', 
                  help='filename that contains workspace to be used cloned ' +\
                      'for use')
parser.add_option('--debug', dest='debug', action='store_true', default=False,
                  help='turn on extra debugging information')
parser.add_option('--noNull', dest='nullFit', action='store_false',
                  default=True, help='do not perform the null hypothesis fit.')
parser.add_option('--btag', dest='btag', action='store_true',
                  default=False, help='Use b-tagged selection.')
parser.add_option('--plotDibosonContributions', dest='plotDibosonContributions', action='store_true',
                  default=False, help='make an extra canvas with WW and WZ contributions shown separately')

(opts, args) = parser.parse_args()

import pyroot_logon
import CMS_lumi
CMS_lumi.lumi_8TeV = "19 fb^{-1}"
CMS_lumi.relPosX = 0.135
#CMS_lumi.writeExtraText = 1
CMS_lumi.writeExtraText = 0
CMS_lumi.extraText = "Preliminary"

#import HWW2DConfig
config = __import__(opts.modeConfig)
import RooWjj2DFitter

from ROOT import TCanvas, TLegend, TLine, RooFit, RooStats, RooArgSet, RooDataSet, RooLinkedListIter, TMath, RooRandom, TFile, TTree, \
         RooDataHist, RooMsgService, TStopwatch, RooAbsPdf, RooAbsReal, TBox, kBlack, kRed, \
         kBlue, kGreen, kOrange, kGray, kDashed, RooAbsCollection
import pulls

timer = TStopwatch()
timer.Start()

if not opts.debug:
    RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)

if hasattr(opts, "seed") and (opts.seed >= 0):
    print "random seed:", opts.seed
    RooRandom.randomGenerator().SetSeed(opts.seed)
pars = config.theConfig(Nj = opts.Nj, mH = opts.mH, 
                        isElectron = opts.isElectron, initFile = args,
                        includeSignal = opts.includeSignal,
                        btagged = opts.btag, includeExtras=opts.plotDibosonContributions)


fitter = RooWjj2DFitter.Wjj2DFitter(pars)
if opts.ws:
    fitter.loadWorkspaceFromFile(opts.ws, getFloatPars = False)    

totalPdf = fitter.makeFitter(opts.useAlternateModelsToFit)


fitter.readParametersFromFile()
fitter.expectedFromPars()
fitter.resetYields()

startpars = totalPdf.getParameters(fitter.ws.set('obsSet'))
print 'startpars:'
startpars.Print('v')
fitter.ws.defineSet("params", startpars)
fitter.ws.saveSnapshot("initPars", startpars)


if opts.runPdfGenToySim:
    parFiles = args
    if opts.genParamFile=='NONE':
        print "Toy generation parameters taken from ", args
    else:
        print "Toy generation parameters taken from ", opts.genParamFile
        parFiles = opts.genParamFile
    #print "performing Toy Generation from param files: ", parFiles
    genPars = config.theConfig(Nj = opts.Nj, mH = opts.mH, isElectron = opts.isElectron, initFile = parFiles, includeSignal = opts.includeSignal)
    genFitter = RooWjj2DFitter.Wjj2DFitter(genPars)

    


    genPdf = genFitter.makeFitter(opts.useAlternateModelsToGen)
    genFitter.readParametersFromFile(parFiles)
    genFitter.initFromExplicitVals(opts)
    genFitter.expectedFromPars()
    genFitter.resetYields()
    
    genstartpars = genPdf.getParameters(genFitter.ws.set('obsSet'))
    genFitter.ws.defineSet("params", genstartpars)
    genFitter.ws.saveSnapshot("genInitPars", genstartpars)
    print "genPars:"
    genFitter.ws.set('params').Print('v')

    if opts.extendedGen:
        data = genPdf.generate(genFitter.ws.set('obsSet'), RooFit.Name('data_obs'),RooFit.Extended())
    else:
        data = genPdf.generate(genFitter.ws.set('obsSet'), RooFit.Name('data_obs'))
    
    if fitter.pars.binData:
        data = RooDataHist('data_obs', 'data_obs', genFitter.ws.set('obsSet'), data)
        data.Print('v')
    getattr(fitter.ws, 'import')(data)
    
    print 'Generated Data Events: %.0f' % (data.sumEntries())
    
elif opts.runMCGenToySim:
    print "performing ToyMC Generation:"
    print "Generating ", opts.extdiboson, " diboson events"
    dibosonfiles = getattr(pars, 'dibosonFiles')
    dibosontoymodels = [-1]
    dibosonPdf = fitter.makeComponentPdf('dibosontoy',dibosonfiles,dibosontoymodels)
    dibosontoymc = dibosonPdf.generate(fitter.ws.set('obsSet'),opts.extdiboson,RooFit.Name('data_obs'))
    dibosontoymc.Print()
    gentoymc=dibosontoymc
    print "Generating ", opts.extWpJ, " WpJ events"
    WpJfiles = getattr(pars, 'WpJFiles')
    WpJtoymodels = [-1]
    WpJPdf = fitter.makeComponentPdf('WpJtoy',WpJfiles,WpJtoymodels)
    WpJtoymc = RooDataSet()   
    WpJtoymc = WpJPdf.generate(fitter.ws.set('obsSet'),opts.extWpJ,RooFit.Name('data_obs'))
    WpJtoymc.Print()
    gentoymc.append(WpJtoymc)
    print "Generating ", opts.exttop, " top events"
    topfiles = getattr(pars, 'topFiles')
    toptoymodels = [-1]
    topPdf = fitter.makeComponentPdf('toptoy',topfiles,toptoymodels)
    toptoymc = RooDataSet()
    toptoymc = topPdf.generate(fitter.ws.set('obsSet'),opts.exttop,RooFit.Name('data_obs'))
    toptoymc.Print()
    gentoymc.append(toptoymc)
    data = gentoymc
    if fitter.pars.binData:
        data = RooDataHist('data_obs', 'data_obs', fitter.ws.set('obsSet'),
                           data)
        data.Print('v')
    getattr(fitter.ws, 'import')(data)
elif opts.toy:
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

#fitter.setMultijetYield()
data.Print()
fitter.resetYields()
startpars.IsA().Destructor(startpars)

print 'Time elapsed: %.1f sec' % timer.RealTime()
print 'CPU time used: %.1f sec' % timer.CpuTime()
print 'starting fitting routine'
timer.Continue()
#fitter.ws.var('top_nrm').setConstant()
fitter.ws.var('r_signal').setVal(1.0)
fitter.ws.var('r_signal').setError(0.04)
fr = None
fr = fitter.fit()

try:
    vName = pars.varNames[pars.var[0]]
except AttributeError:
    vName = pars.var[0]

xvar = fitter.ws.var(vName)

plot1 = fitter.stackedPlot(pars.varNames[pars.var[0]])
chanlabel = '#mu, Data'
if opts.isElectron:
    chanlabel = 'e, Data'
leg1 = RooWjj2DFitter.Wjj2DFitter.legend4Plot(plot1,False,chanlabel)
plot1.addObject(leg1)
minline1 = RooWjj2DFitter.Wjj2DFitter.vline4Plot(plot1.GetMinimum(),plot1.GetMaximum(),pars.sigRegionMin,yScale=0.695)
maxline1 = RooWjj2DFitter.Wjj2DFitter.vline4Plot(plot1.GetMinimum(),plot1.GetMaximum(),pars.sigRegionMax,yScale=0.695)
plot1.addObject(minline1)
plot1.addObject(maxline1)
c1 = TCanvas('c1', fitter.ws.var(pars.varNames[pars.var[0]]).GetTitle() + ' plot')
plot1.Draw()
CMS_lumi.CMS_lumi(c1,2,0)
c1.Update()

print 'Contributions in the diboson signal region of', fitter.pars.sigRegionMin, '-', fitter.pars.sigRegionMax, 'GeV :'
fitter.ws.var(vName).setRange('signalRegion', fitter.pars.sigRegionMin , fitter.pars.sigRegionMax)
for comp in fitter.pars.backgrounds:
    sigInt = fitter.ws.pdf('%s' % (comp)).createIntegral(fitter.ws.set('obsSet'), fitter.ws.set('obsSet'), 'signalRegion')
    fullInt = fitter.ws.pdf('%s' % (comp)).createIntegral(fitter.ws.set('obsSet'), fitter.ws.set('obsSet'))
    n_Exp_sigRegion = sigInt.getVal()/fullInt.getVal()*fitter.ws.function('f_%s_norm' % (comp) ).getVal()
    n_Exp_sigRegion_Err = sigInt.getVal()/fullInt.getVal()*fitter.ws.function('f_%s_norm' % (comp) ).getPropagatedError(fr)    
    print '%s integral in signal region:' % (comp) , sigInt.getVal()/fullInt.getVal()
    print '%s yield in the signal region:' % (comp) , n_Exp_sigRegion, ' +/- ', n_Exp_sigRegion_Err


#Make the Data-NonDiboson subtracted plot
#print 'plot1 : '
#plot1.Print()
xvar = fitter.ws.var(pars.varNames[pars.var[0]])
xvar.setRange('plotRange', xvar.getMin(), xvar.getMax())
dibosonSubtractedFrame = xvar.frame()
dibosonSubtractedFrame.SetName("%s_subtracted" % pars.varNames[pars.var[0]])
dibosonResidual = plot1.residHist('theData', pars.backgrounds[1], False, True)#The first background is the diboson
dibosonResidual.SetTitle('#mu, Data-Bkg')
if opts.isElectron:
    dibosonResidual.SetTitle('e, Data-Bkg')
dibosonSubtractedFrame.addPlotable(dibosonResidual, 'p', False, True)
fitter.ws.pdf('diboson').plotOn(dibosonSubtractedFrame)
dibosonSubtractedFrame.getCurve().SetTitle(pars.dibosonPlotting['title'])
dibosonSubtractedLegend = RooWjj2DFitter.Wjj2DFitter.legend4Plot(dibosonSubtractedFrame)
c2 = TCanvas('c2', xvar.GetTitle() + ' Subtracted')
dibosonSubtractedFrame.GetYaxis().SetTitle(plot1.GetYaxis().GetTitle())
dibosonSubtractedFrame.addObject(dibosonSubtractedLegend)
minline2 = RooWjj2DFitter.Wjj2DFitter.vline4Plot(dibosonSubtractedFrame.GetMinimum(),dibosonSubtractedFrame.GetMaximum(),pars.sigRegionMin,yScale=0.77)
maxline2 = RooWjj2DFitter.Wjj2DFitter.vline4Plot(dibosonSubtractedFrame.GetMinimum(),dibosonSubtractedFrame.GetMaximum(),pars.sigRegionMax,yScale=0.77)
dibosonSubtractedFrame.addObject(minline2)
dibosonSubtractedFrame.addObject(maxline2)
dibosonSubtractedFrame.Draw()
CMS_lumi.CMS_lumi(c2,2,0)
c2.Update()
if opts.plotDibosonContributions:
    dibosonSplitFrame = xvar.frame()
    dibosonSplitFrame.SetName("%s_splitDiboson" % vName)
    dibosonSplitFrame.addPlotable(dibosonResidual, 'p', False, True)
    fitter.ws.pdf('diboson').plotOn(dibosonSplitFrame,RooFit.DrawOption('LF'),RooFit.FillStyle(1001),RooFit.FillColor(kBlue+2),RooFit.LineColor(kBlue+2),RooFit.VLines(),RooFit.Range('plotRange'),RooFit.NormRange('plotRange'))
    dibosonSplitFrame.getCurve().SetTitle('WW')
    WZfiles = getattr(pars, 'WZFiles')
    WZplotmodels = getattr(pars, 'WZModels')
    WZplotPdf = fitter.makeComponentPdf('WZ',WZfiles,WZplotmodels,False,[-1])
##         WZplotPdf.plotOn(dibosonSplitFrame,RooFit.DrawOption('LF'), RooFit.FillStyle(1001),RooFit.FillColor(kGreen+3),RooFit.LineColor(kGreen+3),RooFit.VLines(),RooFit.Range('plotRange'),RooFit.NormRange('plotRange'),RooFit.Normalization(fitter.WZExpected*fitter.ws.var('diboson_nrm').getVal(), RooAbsReal.NumEvent))
    WZplotPdf.plotOn(dibosonSplitFrame,RooFit.DrawOption('LF'), RooFit.FillStyle(1001),RooFit.FillColor(kGreen+3),RooFit.LineColor(kGreen+3),RooFit.VLines(),RooFit.Range('plotRange'),RooFit.NormRange('plotRange'),RooFit.Normalization(fitter.ws.var('n_WZ').getVal()*fitter.ws.var('diboson_nrm').getVal(), RooAbsReal.NumEvent))
    dibosonSplitFrame.getCurve().SetTitle('WZ')
    dibosonSplitLegend = RooWjj2DFitter.Wjj2DFitter.legend4Plot(dibosonSplitFrame,True)
    dibosonSplitFrame.addPlotable(dibosonResidual, 'p', False, True)
    cWWvsWZ = TCanvas('cWWvsWZ', xvar.GetTitle() + ' Split')
    dibosonSplitFrame.GetYaxis().SetTitle(plot1.GetYaxis().GetTitle())
    dibosonSplitFrame.addObject(dibosonSplitLegend)
    dibosonSplitFrame.Draw()
    CMS_lumi.CMS_lumi(cWWvsWZ,2,0)
    cWWvsWZ.Update()

ndf = 0

if fr:
    floatVars = [ fr.floatParsFinal().at(i).GetName() \
                      for i in range(0, fr.floatParsFinal().getSize()) ]
    fitter.ws.defineSet('floatingParams', ','.join(floatVars))
    fitter.ws.saveSnapshot("fitPars", ','.join(floatVars))
    ndf = fr.floatParsFinal().getSize() - \
        fitter.ws.set('constraintSet').getSize()
    fr.Print('v')

print '%i free parameters in the fit' % ndf

firstCurve1 = plot1.getObject(1)
# firstCurve1.Print()

(chi2_1, ndf_1) = pulls.computeChi2(plot1.getHist('theData'), firstCurve1)
pull1 = pulls.createPull(plot1.getHist('theData'), firstCurve1)

chi2 = chi2_1
ndf = ndf_1-ndf

cp1 = TCanvas("cp1", fitter.ws.var(pars.varNames[pars.var[0]]).GetTitle() + ' pull')

pull1.Draw('ap')
minline3 = RooWjj2DFitter.Wjj2DFitter.vline4Plot(-1.9,1.9,pars.sigRegionMin) #Provide the y coordinates of the lines manually
maxline3 = RooWjj2DFitter.Wjj2DFitter.vline4Plot(-1.9,1.9,pars.sigRegionMax)
minline3.Draw()
maxline3.Draw()

## Make the legend "by hand"
legPull=leg1.Clone()
legPull.Clear()
if opts.isElectron:
    legPull.AddEntry(dibosonResidual,"e, Pull","p")
else:
    legPull.AddEntry(dibosonResidual,"#mu, Pull","p")
legPull.SetX1NDC(0.74)
legPull.SetY1NDC(0.8)
legPull.Draw()

pull1.SetName(pars.varNames[pars.var[0]] + "_pull")
cp1.SetGridy()
cp1.Update()
pull1.GetXaxis().SetLimits(pars.varRanges[pars.var[0]][1], 
                           pars.varRanges[pars.var[0]][2])
pull1.GetXaxis().SetTitle(fitter.ws.var(pars.varNames[pars.var[0]]).getTitle(True).Data())
pull1.GetYaxis().SetTitle("pull (#sigma)")
CMS_lumi.CMS_lumi(cp1,2,0)
cp1.Update()

if opts.toyOut:
    outFile = open(opts.toyOut, 'a', 1)
    fitter.ws.loadSnapshot("genPars")

    line = '%s %.6g %.6g %.6g '
    outFile.write('nll %f covQual %i edm %.4g ' % (fr.minNll(), fr.covQual(),
                                                   fr.edm())
                  )
    outFile.write('chi2Prob %f ' % (TMath.Prob(chi2, ndf)))
    for cycle in range(0, fr.numStatusHistory()):
        outFile.write('%s %i ' % (fr.statusLabelHistory(cycle),
                                  fr.statusCodeHistory(cycle))
                      )
    for i in range(0, fr.floatParsFinal().getSize()):
        outFile.write(line % (fr.floatParsFinal().at(i).GetName(),
                              fr.floatParsFinal().at(i).getVal(),
                              fr.floatParsFinal().at(i).getError(),
                              fitter.ws.var(fr.floatParsFinal().at(i).GetName()).getVal()
                              )
                      )
    outFile.write('\n')
    outFile.close()

##write the output parameters to one file
if hasattr(opts, "fitparfn"):
    outparams = totalPdf.getParameters(data)
    finalPars = outparams.snapshot()
    finalPars.writeToFile(opts.fitparfn)


if opts.nullFit:
    print 'doing fit under null hypothesis'

    fitter_null = RooWjj2DFitter.Wjj2DFitter(pars)
    fitter_null.loadHistogramsFromWorkspace(fitter.ws)
    totalPdf_null = fitter_null.makeFitter(opts.useAlternateModelsToFit)
    fitter_null.readParametersFromFile()
    fitter_null.expectedFromPars()
    fitter_null.resetYields()

    fitter_null.loadDataFromWorkspace(fitter.ws)

    fitter_null.ws.var('diboson_nrm').setVal(0.)
    fitter_null.ws.var('diboson_nrm').setConstant()

    fr_null = fitter_null.fit()
    fr_null.SetName('fitResult_null')
    fr_null.Print('v')


mode = 'muon'
if opts.isElectron:
    mode = 'electron'

output = TFile("DibosonBoostedlnuJ_%s_output.root" % mode ,"recreate")


plot1.Write()
dibosonSubtractedFrame.Write()
if opts.plotDibosonContributions:
    dibosonSplitFrame.Write()

pull1.Write()
fitter.ws.SetName("w")
fitter.ws.Write()
if opts.nullFit:
    fr_null.Write()
#fitter.ws.Print()
output.Close()

c1.SaveAs("DibosonBoostedlnuJ_%s_Stacked.png" % mode)
c2.SaveAs("DibosonBoostedlnuJ_%s_Subtracted.png" % mode)
cp1.SaveAs("DibosonBoostedlnuJ_%s_Pull.png" % mode)
if opts.plotDibosonContributions:
    cWWvsWZ.SaveAs("DibosonBoostedlnuJ_%s_Split.png" % mode)
fitter.ws.var('n_diboson').Print()
fitter.ws.var('diboson_nrm').Print()

print 'Time elapsed: %.1f sec' % timer.RealTime()
print 'CPU time used: %.1f sec' % timer.CpuTime()

print '%i degrees of freedom' % ndf
print 'chi2: %.2f / %i = %.2f' % (chi2_1, ndf, (chi2/ndf))
print 'chi2 probability: %.4g' % (TMath.Prob(chi2, ndf))

if opts.nullFit:
    likelihoodRatio = 2.*fr_null.minNll()-2.*fr.minNll()
    print '2*nll_null - 2*nll: %.4f - %.4f = %.4f' % (2.*fr_null.minNll(),
                                                      2.*fr.minNll(),
                                                      likelihoodRatio)

    pval = TMath.Prob(likelihoodRatio, 1)
    print 'p-value: %.4g' % pval
    print 'Gaussian significance: %.3g' % RooStats.PValueToSignificance(pval)

