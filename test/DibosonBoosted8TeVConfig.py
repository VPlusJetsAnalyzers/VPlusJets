from RooWjj2DFitterPars import Wjj2DFitterPars
from ROOT import kRed, kAzure, kGreen, kBlue, kCyan, kViolet, kGray

def theConfig(Nj, mH, isElectron = False, initFile = [], includeSignal = True,
              btagged = False, includeExtras = False):
    pars = Wjj2DFitterPars()
    # pars.MCDirectory = '/uscms_data/d2/andersj/Wjj/2012/data/Moriond2013/ReducedTrees/'
    #pars.MCDirectory = '/uscms_data/d1/lnujj/RDTrees_BoostedW_2013_1_29/'
    #pars.MCDirectory = 'root://cmseos:1094//eos/uscms/store/user/lnujj/BoostedDibosonFitPostMoriond2013/'
    #pars.MCDirectory = 'root://cmseos:1094//eos/uscms/store/user/pdudero/lnujj/JetLepSepTrees4WW/'
    #pars.MCDirectory = "root://cmseos:1094//eos/uscms/store/user/lnujj/DibosonFitPostMoriond2013/"
    pars.MCDirectory = 'root://cmseos:1094//eos/uscms/store/user/lnujj/BoostedDibosonFitPostMoriond2014_ilyao/'
    pars.isElectron = isElectron
    pars.btagSelection = btagged
    pars.boostedSelection = True
    pars.useTopSideband = False
    pars.useTopMC = True
    pars.initialParametersFile = initFile
    pars.extras = []

    if includeExtras:
        pars.extras = ['WZ']
    pars.backgrounds = ['diboson', 'top', 'WpJ']
    pars.includeSignal = includeSignal
    pars.signals = []
    pars.constrainShapes = []
    if pars.btagSelection:
        pars.yieldConstraints = {'top' : 0.50, 'WpJ' : 0.50 }
    else:
        #        pars.yieldConstraints = {'top' : 0.07, 'WpJ' : 0.05 }
        if isElectron:
            pars.yieldConstraints = {'top' : 0.10, 'WpJ' : 0.50 }
        else:
            pars.yieldConstraints = {'top' : 0.08, 'WpJ' : 0.50 }
        #pars.constrainShapes = ['WpJ', 'top', 'diboson']
        pars.constrainShapes = ['top', 'diboson']
        
    #pars.yieldConstraints = {}

#    pars.yieldConstraints = {'top' : 0.50, 'WpJ' : 0.50 }
    #pars.constrainShapes = ['WpJ']

    pars.Njets = Nj
    pars.mHiggs = mH

    if isElectron:
        flavorString = 'el'
    else:
        flavorString = 'mu'

    ##Second Jet Veto: &&(GroomedJet_CA8_pt[1]<150)
    ##Second Jet Separation from the lepton: '&&( (GroomedJet_CA8_deltaR_lca8jet[1]<-900.0)||(GroomedJet_CA8_deltaR_lca8jet[1]>7.0) )'
    pars.cuts = \
              '(W_pt>200.)&&(GroomedJet_CA8_pt[0]>200)' +\
              '&&(abs(GroomedJet_CA8_eta[0])<2.4)' +\
              '&&(GroomedJet_CA8_pt[1]<80)' +\
              '&&(GroomedJet_CA8_mass_pr[0]>40)' +\
              '&&(GroomedJet_CA8_tau2tau1[0]<0.55)'
    

    pars.btagVeto = False

    if pars.useTopSideband:
        pars.cuts += '&&(GroomedJet_numberbjets_csvm>=1)'
    # elif pars.btagSelection:
    #     pars.cuts += '&&(GroomedJet_numberbjets==1)'
    else:
        pars.cuts += '&&(ggdboostedWevt==1)' +\
                     '&&(GroomedJet_CA8_deltaphi_METca8jet[0]>2.0)' +\
                     '&&(GroomedJet_CA8_deltaR_lca8jet[0]>1.57)'
        pars.cuts += '&&(numPFCorJetBTags<1)'

        

    # veto vbf
    pars.cuts += '&&(vbf_event==0)'

    # you need a files entry and a models entry for each of the fit 
    # compoents in backgrounds and signals
    # the files should a list with entries like (filename, Ngen, xsec)
    
#################### Global Convolution Models ########################
    #pars.GlobalConvModels=[27]
    pars.GlobalConvModels=[-1]
    pars.GlobalConvModelsAlt=pars.GlobalConvModels
    
#####################  diboson: #######################################
##     pars.dibosonFiles = [
##         (pars.MCDirectory + 'RD_%s_WW_CMSSW532.root' % (flavorString),
##          9450414, 57.25),
##         (pars.MCDirectory + 'RD_%s_WZ_CMSSW532.root' % (flavorString),
##          10000267, 22.88),
##         ]
    
    ### aMC@NLO Cross-Check ###
    ### the cross-sections used are 'absolute value' ones (i.e. the negative weight is propagated into effwt and will be taken into account when producing the distributions)
    ### multiple samples for WZ are used
    pars.dibosonFiles = [
        (pars.MCDirectory + 'RD_%s_WW_minPt150_amcnlo_CMSSW532.root' % (flavorString),
         119692, 2*0.2222038*1.883),
        (pars.MCDirectory + 'RD_%s_WZ_minPt150_amcnlo_CMSSW532.root' % (flavorString),
         95230+55805, 0.230996*1.00125782),
        (pars.MCDirectory + 'RD_%s_WZ_minPt150_amcnlo_add_CMSSW532.root' % (flavorString),
         95230+55805, 0.230996*1.00125782),
        ]

##     ### To get AxEff into the inclusive final state:
##     pars.dibosonFiles = [
##         (pars.MCDirectory + 'RD_%s_WW_minPt150_amcnlo_CMSSW532.root' % (flavorString),
##          119692*57.25/(2*0.2222038*1.883), 57.25),
## ##         (pars.MCDirectory + 'RD_%s_WZ_minPt150_amcnlo_CMSSW532.root' % (flavorString),
## ##          (55805)*22.88/(0.230996*1.00125782), 22.88),
##         (pars.MCDirectory + 'RD_%s_WZ_minPt150_amcnlo_add_CMSSW532.root' % (flavorString),
##          (95230)*22.88/(0.230996*1.00125782), 22.88),
##         ]



##     ### High-Stat aMC@NLO Samples
##     WWfileName = pars.MCDirectory + 'RD_%s_WWtoLNuQQ_amcnlo_800k_WWptResumWts_CMSSW532.root' % (flavorString)
##     pars.dibosonFiles = [
##         (WWfileName, 799899/0.444408, 57.52),
##         (pars.MCDirectory + 'RD_%s_WZtoLNuQQ_amcnlo_100k_CMSSW532.root' % (flavorString), 99791/0.230996, 24.20),
##         ]
    
    pars.dibosonFracOfData = -1
    #pars.dibosonModels = [5]
    pars.dibosonModels = [22]
    pars.dibosonModelsAlt = pars.dibosonModels
    pars.dibosonConvModels = pars.GlobalConvModels
    pars.dibosonConvModelsAlt = pars.dibosonConvModels

    
### Test with aMC@NLO
##     pars.WWaMCNLOFiles = [
##         (pars.MCDirectory + 'RD_%s_WW_amcnlo_CMSSW532.root' % (flavorString),
##          39400, 3*0.727),
## ##         (pars.MCDirectory + 'RD_%s_WW_CMSSW532.root' % (flavorString),
## ##          9450414, 57.25),
##         ]
    
##     pars.WWaMCNLOFracOfData = -1
##     pars.WWaMCNLOModels = [13]
##     pars.WWaMCNLOModelsAlt = pars.WWaMCNLOModels
##     pars.WWaMCNLOConvModels = pars.GlobalConvModels
##     pars.WWaMCNLOConvModelsAlt = pars.WWaMCNLOConvModels

## ### default WW separately: ###
##     pars.WWdefaultFiles = [
##         (pars.MCDirectory + 'RD_%s_WW_CMSSW532.root' % (flavorString),
##          9450414, 57.25),
##         ]
    
##     pars.WWdefaultFracOfData = -1
##     pars.WWdefaultModels = [13]
##     pars.WWdefaultModelsAlt = pars.WWdefaultModels
##     pars.WWdefaultConvModels = pars.GlobalConvModels
##     pars.WWdefaultConvModelsAlt = pars.WWdefaultConvModels


### WZ separately: ###
##     pars.WZFiles = [(pars.MCDirectory + 'RD_%s_WZ_CMSSW532.root' % (flavorString),
##                      10000267, 22.88),
##                     ]
    pars.WZFiles = [(pars.MCDirectory + 'RD_%s_WZ_minPt150_amcnlo_CMSSW532.root' % (flavorString), 95230+55805, 0.230996*1.00125782),
                    (pars.MCDirectory + 'RD_%s_WZ_minPt150_amcnlo_add_CMSSW532.root' % (flavorString), 95230+55805, 0.230996*1.00125782),]

    pars.WZFracOfData = -1
    pars.WZModels = [13]
    pars.WZModelsAlt = pars.WZModels
    pars.WZConvModels = pars.GlobalConvModels
    pars.WZConvModelsAlt = pars.WZConvModels    

#####################  WpJ: ###########################################
##     pars.WpJFiles = [
##         (pars.MCDirectory + 'RD_%s_WJets_madgraph_CMSSW532.root' % (flavorString),
##          8955318, 1.3*228.9),
##         ]

    pars.WpJFiles = [
        (pars.MCDirectory + 'RD_%s_WpJ_PT180_Madgraph_CMSSW532.root' % (flavorString),
         9492452, 1.3*23.5),
        ]
    
##     ### HERWIG Cross-check ###
##     pars.WpJFiles = [
##         (pars.MCDirectory + 'RD_%s_WpJ_HERWIG_CMSSW532.root' % (flavorString),
##          13256090, 1.3*228.9),
##         ] #update crossX info
    
    if pars.btagSelection:
        pars.WpJFracOfData = 0.332
    else:
        pars.WpJFracOfData = -1
##         if isElectron:
##             #pars.WpJFracOfData = 0.733
##             pars.WpJFracOfData = 0.800
##         else:
##             #pars.WpJFracOfData = 0.737
##             pars.WpJFracOfData = 0.782

    pars.WpJModels = [8]
    #    pars.WpJModelsAlt = [8]
    #    pars.WpJModelsAlt = [308] ##Alt2 Model
    pars.WpJModelsAlt = [10] ##Alt3 Model
    pars.WpJAuxModelsAlt = [5]
    pars.WpJConvModels = pars.GlobalConvModels
    pars.WpJConvModelsAlt = pars.WpJConvModels

#####################  top: #######################################  
    ttkfactor = 0.95
    if isElectron:
        ttkfactor = 0.92

    if pars.useTopSideband:
        if isElectron:
            pars.topFiles = [(pars.MCDirectory + 'RD_WenuJets_DataAllSingleElectronTrigger_GoldenJSON_19p2invfb.root',1,1), ]
        else:
            pars.topFiles = [(pars.MCDirectory + 'RD_WmunuJets_DataAll_GoldenJSON_19p3invfb.root',1,1), ]
    else:
        pars.topFiles = [
        (pars.MCDirectory + 'RD_%s_STopTW_Tbar_CMSSW532.root' % (flavorString),
         493458, 11.1773*ttkfactor),
        (pars.MCDirectory + 'RD_%s_STopTW_T_CMSSW532.root' % (flavorString),
         497657, 11.1773*ttkfactor),
##         (pars.MCDirectory + 'RD_%s_TTJetsPoheg_CMSSW532.root' % (flavorString),
##          20975917, 225.197*ttkfactor),
##         (pars.MCDirectory + 'RD_%s_TTJets_poheg_CMSSW532_v2.root' % (flavorString),
##          20975917, 225.197*ttkfactor),
        ('root://cmseos:1094//eos/uscms/store/user/lnujj/VBF_Higgs_v1/VBF_Higgs_5Aug_v1/RD_%s_TTJets_poheg_CMSSW532.root' % (flavorString),
         20975917, 225.197*ttkfactor),
        (pars.MCDirectory + 'RD_%s_STopS_Tbar_CMSSW532.root' % (flavorString),
         139974, 1.75776*ttkfactor),
        (pars.MCDirectory + 'RD_%s_STopS_T_CMSSW532.root' % (flavorString),
         259960, 3.89394*ttkfactor),
        (pars.MCDirectory + 'RD_%s_STopT_Tbar_CMSSW532.root' % (flavorString),
         1935066, 30.0042*ttkfactor),
        (pars.MCDirectory + 'RD_%s_STopT_T_CMSSW532.root' % (flavorString),
         3758221, 55.531*ttkfactor),
        ]
    
    if pars.btagSelection:
        pars.topFracOfData = 0.616
##         if isElectron:
##             pars.topModels = [5]
##         else:
##             pars.topModels = [13]
    else:
        pars.topFracOfData = -1
##         if isElectron:
##             pars.topFracOfData = 0.201
##         else:
##             pars.topFracOfData = 0.199


    pars.topModels = [30]
    pars.topModelsAlt = pars.topModels
#    pars.topModelsAlt = [330]
    pars.topConvModels = pars.GlobalConvModels
    pars.topConvModelsAlt = pars.topConvModels

################################################################### 
    pars.WZPlotting = {'color' : kGreen+3, 'title' : 'WZ'}
    pars.dibosonPlotting = {'color' : kAzure+8, 'title' : 'WW+WZ'}
    pars.WpJPlotting = { 'color' : kRed, 'title' : 'W+jets'}
    pars.topPlotting = {'color' : kGreen+2, 'title' : 'Top'}
    pars.ggHPlotting = {'color' : kBlue, 'title' : "ggH(%i) #rightarrow WW" % mH}

    pars.var = ['GroomedJet_CA8_mass_pr[0]']
    pars.varRanges = {'GroomedJet_CA8_mass_pr[0]': (10, 40., 140., []),}
    pars.sigRegionMin = 70.0
    pars.sigRegionMax = 100.0
    pars.varTitles = {'GroomedJet_CA8_mass_pr[0]': 'm_{J}',
                      }
    pars.varNames = {'GroomedJet_CA8_mass_pr[0]': 'GroomedJet_CA8_mass_pr' }

    
    pars.exclude = {}
    pars.blind = False

    pars.integratedLumi = 19300.

    pars.binData = False
    # pars.binData = True

    return customizeElectrons(pars) if isElectron else \
        customizeMuons(pars)

def customizeElectrons(pars):
    pars.DataFile = pars.MCDirectory + 'RD_WenuJets_DataAllSingleElectronTrigger_GoldenJSON_19p2invfb.root'

    if pars.useTopSideband and not pars.useTopMC:
        pars.topFiles = [(pars.DataFile,1,1),]
    
    pars.integratedLumi = 19200.
    pars.doEffCorrections = True
    pars.effToDo = ['lepton']
    pars.leptonEffFiles = {
        'id': ["EffTable2012/scaleFactor-Run2012ABC-GsfElectronToId.txt"],
        'reco': ["EffTable2012/scaleFactor-Run2012ABC-SCToElectron.txt"],
        'HLT': ["EffTable2012/efficiency-Run2012ABC-WP80ToHLTEle.txt"]
        }
    pars.lumiPerEpoch = [pars.integratedLumi]

##     pars.cuts += '&&(W_electron_pt>30)'
    pars.cuts += '&&(event_met_pfmet >70)&&(W_electron_pt>35)'
    return pars

def customizeMuons(pars):
    pars.DataFile = pars.MCDirectory + 'RD_WmunuJets_DataAll_GoldenJSON_19p3invfb.root'
    
    if pars.useTopSideband and not pars.useTopMC:
        pars.topFiles = [(pars.DataFile,1,1),]
    
    pars.doEffCorrections = True
    pars.effToDo = ['lepton']
    pars.leptonEffFiles = {
        'id': ["EffTable2012/scaleFactor-Run2012ABC-RecoToIso.txt"],
        'HLT': ["EffTable2012/efficiency-Run2012ABC-IsoToIsoMuHLT.txt"]
        }
    pars.lumiPerEpoch = [pars.integratedLumi]

    pars.cuts += '&&(event_met_pfmet >50)&&(abs(W_muon_eta)<2.1)&&(W_muon_pt>30.)'
    
    return pars
