from RooWjj2DFitterPars import Wjj2DFitterPars
from ROOT import kRed, kAzure, kGreen, kBlue, kCyan, kViolet, kGray, kYellow

def theConfig(Nj, mH, isElectron = False, initFile = [], includeSignal = True,
              btagged = False, includeExtras = False):
    pars = Wjj2DFitterPars()

    #pars.MCDirectory = '/uscms_data/d2/andersj/Wjj/2012/data/Moriond2013/ReducedTrees/'
    # pars.MCDirectory = "root://cmseos:1094//eos/uscms/store/user/lnujj/HCP2012METfix/ReducedTrees/"
    #    fitterPars.MCDirectory = "root://cmssrv32.fnal.gov//store/user/lnujj/WjjTrees/Full2011DataFall11MC/OriginalTree_PAT/"
    pars.MCDirectory = "root://cmseos:1094//eos/uscms/store/user/lnujj/DibosonFitPostMoriond2013/"
    #pars.SecondaryDirectory = "root://cmseos:1094//eos/uscms/store/user/ilyao/Diboson8TeV_DibosonFitPostMoriond2013/"
    pars.SecondaryDirectory = pars.MCDirectory
    pars.isElectron = isElectron
    pars.btagSelection = btagged
    pars.boostedSelection = False
    pars.useTopSideband = False
    pars.initialParametersFile = initFile
    pars.extras = []
    pars.Njets = Nj
    pars.mHiggs = 126.
    includeWWpTCorr = True
    WWpTCorrOpt = 0
    

    if isElectron:
        flavorString = 'el'
    else:
        flavorString = 'mu'


    pars.cuts = '(sqrt(JetPFCor_Pt[0]**2+JetPFCor_Pt[1]**2+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))>70.)' + \
        '&&(abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5)' + \
        '&&(abs(JetPFCor_dphiMET[0])>0.4)' + \
        '&&(W_mt>30.)'

    pars.cuts += '&&(JetPFCor_Pt[1]>35.)' + \
        '&&(JetPFCor_Pt[0]>40.)'
##        '&&(JetPFCor_Pt[2]<30.)'

##     ### To get the mass peak in the higgs samples
##     pars.cuts = '(abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5)' + \
##         '&&(abs(JetPFCor_dphiMET[0])>0.4)' + \
##         '&&(W_mt>30.)'

##     pars.cuts += '&&(JetPFCor_Pt[1]>25.)' + \
##         '&&(JetPFCor_Pt[0]>30.)'

    
##    pars.cuts += '&&(JetPFCor_Pt[1]>35.)' + \
##                 '&&(JetPFCor_Pt[0]>40.)'

##    pars.cuts += '&&(JetPFCor_Pt[2]<30.)'

### Cross-checks explicitly changing the top sideband selection (useTopSideband should be set to true)
##    pars.cuts += '&&( ((JetPFCor_Pt[0]>40.)&&(JetPFCor_bDiscriminatorCSV[0]<0.679)&&(JetPFCor_Pt[1]>35.)&&(JetPFCor_bDiscriminatorCSV[1]<0.679)&&(JetPFCor_Pt[2]<30.)&&(JetPFCor_Pt[2]>20##.)&&(JetPFCor_bDiscriminatorCSV[2]>0.679)&&(JetPFCor_Pt[3]<30.)&&(JetPFCor_Pt[3]>20.)&&(JetPFCor_bDiscriminatorCSV[3]>0.679)) )'
##       '||((JetPFCor_Pt[0]>40.)&&(JetPFCor_bDiscriminatorCSV[0]<0.244)&&(JetPFCor_Pt[1]>35.)&&(JetPFCor_bDiscriminatorCSV[1]>0.679)&&(JetPFCor_Pt[2]>30.)&&(JetPFCor_bDiscriminatorCSV[2]<0.244)&&(JetPFCor_Pt[3]>30.)&&(JetPFCor_bDiscriminatorCSV[3]>0.679))' + \
##       '||((JetPFCor_Pt[0]>40.)&&(JetPFCor_bDiscriminatorCSV[0]>0.679)&&(JetPFCor_Pt[1]>35.)&&(JetPFCor_bDiscriminatorCSV[1]<0.244)&&(JetPFCor_Pt[2]>30.)&&(JetPFCor_bDiscriminatorCSV[2]<0.244)&&(JetPFCor_Pt[3]>30.)&&(JetPFCor_bDiscriminatorCSV[3]>0.679))' + \
##        '||((JetPFCor_Pt[0]>40.)&&(JetPFCor_bDiscriminatorCSV[0]<0.244)&&(JetPFCor_Pt[1]>35.)&&(JetPFCor_bDiscriminatorCSV[1]>0.679)&&(JetPFCor_Pt[2]>30.)&&(JetPFCor_bDiscriminatorCSV[2]>0.679)&&(JetPFCor_Pt[3]>30.)&&(JetPFCor_bDiscriminatorCSV[3]<0.244))' + \
##        '||((JetPFCor_Pt[0]>40.)&&(JetPFCor_bDiscriminatorCSV[0]>0.679)&&(JetPFCor_Pt[1]>35.)&&(JetPFCor_bDiscriminatorCSV[1]<0.244)&&(JetPFCor_Pt[2]>30.)&&(JetPFCor_bDiscriminatorCSV[2]>0.679)&&(JetPFCor_Pt[3]>30.)&&(JetPFCor_bDiscriminatorCSV[3]<0.244))' + \
##        '||((JetPFCor_Pt[0]>40.)&&(JetPFCor_bDiscriminatorCSV[0]>0.679)&&(JetPFCor_Pt[1]>35.)&&(JetPFCor_bDiscriminatorCSV[1]>0.679)&&(JetPFCor_Pt[2]>30.)&&(JetPFCor_bDiscriminatorCSV[2]<0.244)&&(JetPFCor_Pt[3]>30.)&&(JetPFCor_bDiscriminatorCSV[3]<0.244))' + \
##        ' )'

### Cross-checks fitting the W+ and W- regions separately
##    pars.cuts += '&&(W_muon_charge>0.0)'
##    pars.cuts += '&&(W_muon_charge<0.0)'    
    
 
    #implement topselection, btagged or anti-btagged cuts
    pars.btagVeto = False
    if pars.useTopSideband:
        pars.cuts += '&&(TopWm>0)&&(JetPFCor_Pt[4]<30.)'
##        pars.cuts +='&&(JetPFCor_Pt[4]<30.)'
    else:
        pars.cuts += '&&(JetPFCor_Pt[2]<30.)'
        
        startRange = 0
        ##startRange = 2 ##Test with no b-tagging
        if pars.btagSelection:
            startRange = 2
            for i in range(0, startRange):
                pars.cuts += '&&(JetPFCor_bDiscriminatorCSV[%i]>0.244)' % i
                #pars.cuts += '&&(JetPFCor_bDiscriminatorCSV[%i]>0.679)' % i
        for i in range(startRange, 6):
            pars.cuts += '&&((abs(JetPFCor_Eta[%i])>2.4)||' % i + \
                         '(JetPFCor_Pt[%i]<30.)||' % i + \
                         '(JetPFCor_bDiscriminatorCSV[%i]<0.244))' % i  
            

## #### Attempts to select top events near the signal region by looking for additional b-tagged jets (3rd jet and above with 20 < pt < 30). ####
##        pars.cuts += '&&((abs(JetPFCor_Eta[0])>2.4)||(JetPFCor_bDiscriminatorCSV[0]<0.244))&&((abs(JetPFCor_Eta[1])>2.4)||(JetPFCor_bDiscriminatorCSV[1]<0.244))'
##        pars.cuts += '&&(  ' +\
##            '( (JetPFCor_bDiscriminatorCSV[2]>0.898)&&(JetPFCor_Pt[2]<30.)&&(JetPFCor_Pt[2]>20.) &&((abs(JetPFCor_Eta[3])>2.4)||(JetPFCor_Pt[3]<30.)||(JetPFCor_bDiscriminatorCSV[3]<0.244))&&((abs(JetPFCor_Eta[4])>2.4)||(JetPFCor_Pt[4##]<30.)||(JetPFCor_bDiscriminatorCSV[4]<0.244))&&((abs(JetPFCor_Eta[5])>2.4)||(JetPFCor_Pt[5]<30.)||(JetPFCor_bDiscriminatorCSV[5]<0.244)) )' +\
##            '||( (JetPFCor_bDiscriminatorCSV[3]>0.898)&&(JetPFCor_Pt[3]<30.)&&(JetPFCor_Pt[3]>20.) &&((abs(JetPFCor_Eta[2])>2.4)||(JetPFCor_Pt[2]<30.)||(JetPFCor_bDiscriminatorCSV[2]<0.244))&&((abs(JetPFCor_Eta[4])>2.4)||(JetPFCor_Pt##[4]<30.)||(JetPFCor_bDiscriminatorCSV[4]<0.244))&&((abs(JetPFCor_Eta[5])>2.4)||(JetPFCor_Pt[5]<30.)||(JetPFCor_bDiscriminatorCSV[5]<0.244)) )' +\
##            '||( (JetPFCor_bDiscriminatorCSV[4]>0.898)&&(JetPFCor_Pt[4]<30.)&&(JetPFCor_Pt[4]>20.) &&((abs(JetPFCor_Eta[3])>2.4)||(JetPFCor_Pt[3]<30.)||(JetPFCor_bDiscriminatorCSV[3]<0.244))&&((abs(JetPFCor_Eta[2])>2.4)||(JetPFCor_Pt##[2]<30.)||(JetPFCor_bDiscriminatorCSV[2]<0.244))&&((abs(JetPFCor_Eta[5])>2.4)||(JetPFCor_Pt[5]<30.)||(JetPFCor_bDiscriminatorCSV[5]<0.244)) )' +\
##            '||( (JetPFCor_bDiscriminatorCSV[5]>0.898)&&(JetPFCor_Pt[5]<30.)&&(JetPFCor_Pt[5]>20.) &&((abs(JetPFCor_Eta[3])>2.4)||(JetPFCor_Pt[3]<30.)||(JetPFCor_bDiscriminatorCSV[3]<0.244))&&((abs(JetPFCor_Eta[4])>2.4)||(JetPFCor_Pt##[4]<30.)||(JetPFCor_bDiscriminatorCSV[4]<0.244))&&((abs(JetPFCor_Eta[2])>2.4)||(JetPFCor_Pt[2]<30.)||(JetPFCor_bDiscriminatorCSV[2]<0.244)) )' +\
##            '  )'



    # veto boosted topology
    # pars.cuts += '&&(ggdboostedWevt==0)&&(W_pt<200.)'
    pars.cuts += '&&(W_pt<200.)'

    # veto vbf
    pars.cuts += '&&(vbf_event==0)'


    if includeExtras:
        pars.extras = ['WZ']
    pars.backgrounds = ['diboson', 'top', 'WpJ']
    pars.yieldConstraints = {'top' : 0.06, 'WpJ' : 0.05}
##    pars.yieldConstraints = {'top' : 1.0, 'WpJ' : 0.0001}##Saschas Test 
##    pars.yieldConstraints = {'top' : 0.20, 'WpJ' : 0.05} ##Relaxed Constraint Cross-Check
    pars.constrainShapes = ['WpJ']
    
    if pars.useTopSideband:
        pars.backgrounds = ['top', 'WpJ']
        pars.yieldConstraints = {}

    if pars.btagSelection:
        if includeExtras:
            pars.extras = ['WZ']
        pars.backgrounds = ['diboson', 'WHbb', 'top', 'WpJ']
        pars.yieldConstraints = {'WHbb' : 0.000001, 'top' : 0.06, 'WpJ' : 0.05}
##        pars.yieldConstraints = {'WHbb' : 0.000001, 'top' : 0.20, 'WpJ' : 0.05} ##Relaxed Constraint Cross-Check
        pars.constrainShapes = ['WpJ']
        

    # you need a files entry and a models entry for each of the fit 
    # compoents in backgrounds and signals
    # the files should a list with entries like (filename, Ngen, xsec)
#################### Global Convolution Models ####################
    #pars.GlobalConvModels=[27]
    pars.GlobalConvModels=[-1]
    pars.GlobalConvModelsAlt=pars.GlobalConvModels
    
#####################  diboson: #######################################
##     ### PYTHIA Samples:
##     pars.dibosonFiles = [
##         (pars.MCDirectory + 'RD_%s_WW_CMSSW532.root' % (flavorString),
##          9450414, 57.25),
##         (pars.MCDirectory + 'RD_%s_WZ_CMSSW532.root' % (flavorString),
##          10000267, 22.88),
##         ]

##     ### PYTHIA Samples with Scale Up and Scale Down fluctuations
##     pars.dibosonFiles = [
##         (pars.MCDirectory + 'RD_%s_WW_CMSSW532.root' % (flavorString), 9450414, 57.25),
##         (pars.MCDirectory + 'RD_%s_WZ_CMSSW532.root' % (flavorString),
##          10000267, 22.88),
##         ]

    ### High-Stat aMC@NLO Samples
    WWfileName = pars.MCDirectory + 'RD_%s_WWtoLNuQQ_amcnlo_800k_WWptResumWts_CMSSW532.root' % (flavorString)

    ## with WWpT corrections included
    if includeWWpTCorr:
        if WWpTCorrOpt == 0:
            WWfileName = pars.MCDirectory + 'FT_%s_WWtoLNuQQ_amcnlo_800k_WWptResumWts_CMSSW532_NominalWWpTwt.root' % (flavorString)
        elif WWpTCorrOpt == 1:
            WWfileName = pars.MCDirectory + 'FT_%s_WWtoLNuQQ_amcnlo_800k_WWptResumWts_CMSSW532_MUWWpTwt.root' % (flavorString)
        elif WWpTCorrOpt == -1:
            WWfileName = pars.MCDirectory + 'FT_%s_WWtoLNuQQ_amcnlo_800k_WWptResumWts_CMSSW532_MDWWpTwt.root' % (flavorString)
        elif WWpTCorrOpt == 2:
            WWfileName = pars.MCDirectory + 'FT_%s_WWtoLNuQQ_amcnlo_800k_WWptResumWts_CMSSW532_SUWWpTwt.root' % (flavorString)
        elif WWpTCorrOpt == -2:
            WWfileName = pars.MCDirectory + 'FT_%s_WWtoLNuQQ_amcnlo_800k_WWptResumWts_CMSSW532_SDWWpTwt.root' % (flavorString)           


##     pars.dibosonFiles = [
##         ('/uscms_data/d3/ilyao/Diboson8TeV/FilteredTrees/RD_mu_WWtoLNuQQ_amcnlo_800k_WWptResumWts_CMSSW532_NoWWpTwt.root', 799899, 57.52*0.444408),
##         (pars.MCDirectory + 'RD_%s_WZtoLNuQQ_amcnlo_100k_CMSSW532.root' % (flavorString), 99791, 24.20*0.230996),
##         ]

##     pars.dibosonFiles = [
##         (WWfileName, 799899, 57.52*0.444408),
##         (pars.MCDirectory + 'RD_%s_WZtoLNuQQ_amcnlo_100k_CMSSW532.root' % (flavorString), 99791, 24.20*0.230996),
##         ]

### To Compute AxEff properly: need to take into account the fact that a) we want AxEff for the final state covering the full phase space and that b) the code adds the cross-sections when computing AxEff. I.e. Nexpected = AxEff_total*L*CrossX_total=L*Sum_i(AxEff_i*CrossX_i). CrossX_i=fracPhaseSpace_i*CrossX_WW, AxEff_i=N_i_passing/N_i_generated. Take CrossX_i->CrossX_i/Sum_i(fracPhaseSpace_i), N_i_generated->N_i_generated/Sum_i(fracPhaseSpace_i). Then CrossX_total=CrossX_WW, while Nexpected is unchanged. Similarly for WZ.
    CrossX_WW_aMCNLO=57.52
    fracPh_WWmain=0.444408
    fracPh_WWtauleptonic=0.0246772*2
    sumfracPh_WW=fracPh_WWmain+fracPh_WWtauleptonic
    CrossX_WZ_aMCNLO=24.20
    fracPh_WZmain=0.230996
    fracPh_WZtauleptonic=0.0114347
    sumfracPh_WZ=fracPh_WZmain+fracPh_WZtauleptonic
### Take the relative fractions of WW to WZ from the same MC (aMCNLO), but the ovarall cross section from theory
    CrossXTheory_diboson=59.8+22.9
    CrossX_WW=CrossXTheory_diboson*CrossX_WW_aMCNLO/(CrossX_WW_aMCNLO+CrossX_WZ_aMCNLO)
    CrossX_WZ=CrossXTheory_diboson*CrossX_WZ_aMCNLO/(CrossX_WW_aMCNLO+CrossX_WZ_aMCNLO)
    
    pars.dibosonFiles = [
        (WWfileName, 799899/sumfracPh_WW, CrossX_WW*fracPh_WWmain/sumfracPh_WW),
        (pars.MCDirectory + 'RD_%s_WW_lvtauvtau_CMSSW532.root' % (flavorString),(24898+24898)/sumfracPh_WW, CrossX_WW*fracPh_WWtauleptonic/sumfracPh_WW),
        (pars.MCDirectory + 'RD_%s_WZtoLNuQQ_amcnlo_100k_CMSSW532.root' % (flavorString), 99791/sumfracPh_WZ, CrossX_WZ*fracPh_WZmain/sumfracPh_WZ),
        (pars.MCDirectory + 'RD_%s_WZ_lvtautau_CMSSW532.root' % (flavorString),(24897+24899)/sumfracPh_WZ, CrossX_WZ*fracPh_WZtauleptonic/sumfracPh_WZ),
        ]

##     pars.dibosonFiles = [
##         (pars.MCDirectory + 'RD_%s_WWtoLNuQQ_amcnlo_800k_CMSSW532.root' % (flavorString), 799899, 57.52*0.444408),
##         (pars.MCDirectory + 'RD_%s_WZtoLNuQQ_amcnlo_100k_CMSSW532.root' % (flavorString), 99791, 24.20*0.230996),
##         ]



    

## ### Test with PowHeg Using aMC@NLO Cross-Section (WZ is an inclusive sample)
##     pars.dibosonFiles = [
##         (pars.MCDirectory + 'RD_%s_WWtoLnuQQ_powheg2_CMSSW532.root' % (flavorString), 100000, 57.52*0.444408),
##         (pars.MCDirectory + 'RD_%s_WZtoAnything_powheg2_50k_CMSSW532.root' % (flavorString), 50000, 24.20),
##         ]


## ### Test with aMC@NLO  
##     pars.dibosonFiles = [
##         (pars.MCDirectory + 'RD_%s_WW_aMCNLO_nocut_CMSSW532.root' % (flavorString),
##          49497, 57.52*0.444408),
##         (pars.MCDirectory + 'RD_%s_WZ_CMSSW532.root' % (flavorString),
##          10000267, 22.88),
##         ]

## ### Test with MADGRAPH at LO  
##     pars.dibosonFiles = [
##         (pars.MCDirectory + 'RD_%s_WW_MGLO_nocut99K_CMSSW532.root' % (flavorString),
##          99000, 50.59*0.444408),
##         (pars.MCDirectory + 'RD_%s_WZ_CMSSW532.root' % (flavorString),
##          10000267, 22.88),
##         ]


## ### Test with aMC@NLO (the cross-sections are abs valued with the negative event weights propagated via effwt
##     pars.dibosonFiles = [
##         (pars.MCDirectory + 'RD_%s_WW_minPt20_amcnloQiang_CMSSW532.root' % (flavorString),
##          250000, 2*0.22204*47.99),
##         (pars.MCDirectory + 'RD_%s_WZ_minPt20_amcnloQiang_CMSSW532.root' % (flavorString),
##          100000, 0.230996*20.09),
##         ]
    
    pars.dibosonFracOfData = -1
    pars.dibosonModels = [22]
    #pars.dibosonModels = [-1]
    pars.dibosonModelsAlt = pars.dibosonModels
    pars.dibosonConvModels = pars.GlobalConvModels
    pars.dibosonConvModelsAlt = pars.dibosonConvModels

### WZ separately: ###
##     pars.WZFiles = [(pars.MCDirectory + 'RD_%s_WZ_CMSSW532.root' % (flavorString),
##                      10000267, 22.88),
##                     ]
    pars.WZFiles = [
        (pars.MCDirectory + 'RD_%s_WZtoLNuQQ_amcnlo_100k_CMSSW532.root' % (flavorString), 99791/sumfracPh_WZ, CrossX_WZ*fracPh_WZmain/sumfracPh_WZ),
        (pars.MCDirectory + 'RD_%s_WZ_lvtautau_CMSSW532.root' % (flavorString),(24897+24899)/sumfracPh_WZ, CrossX_WZ*fracPh_WZtauleptonic/sumfracPh_WZ),
        ]
    pars.WZFracOfData = -1
    pars.WZModels = [13]
    pars.WZModelsAlt = pars.WZModels
    pars.WZConvModels = pars.GlobalConvModels
    pars.WZConvModelsAlt = pars.WZConvModels


##     ## sideband region studies
##     pars.dibosonFracOfData = -1
##     pars.dibosonModels = [-1]


    
#####################  WpJ: #######################################
    wpj_kfactor = 1.16
    if pars.btagSelection:
        wpj_kfactor = 2.05 # from arXiv:1011.6647 [hep-ph]
    pars.WpJFiles = [
##         (pars.MCDirectory + 'RD_%s_WpJ_CMSSW532.root' % (flavorString),
##          18353019+50768992, 36257.2),
        (pars.MCDirectory + 'RD_%s_W2Jets_CMSSW532.root' % (flavorString),
         33004921, 1750.0*wpj_kfactor),
        (pars.MCDirectory + 'RD_%s_W3Jets_CMSSW532.root' % (flavorString),
         15059503, 519.0*wpj_kfactor),
        (pars.MCDirectory + 'RD_%s_W4Jets_CMSSW532.root' % (flavorString),
         12842803, 214.0*wpj_kfactor),
        (pars.MCDirectory + 'RD_%s_W1Jets_CMSSW532.root' % (flavorString),
         19871598, 5400.0*wpj_kfactor),
        (pars.MCDirectory + 'RD_%s_ZpJ_CMSSW532.root' % (flavorString),
         30209426, 3503.71),
        ]

## #### Check WGamma* with V+Jets Parameterization
##     pars.WpJFiles = [(pars.MCDirectory + 'RD_%s_qcdwa0-3jets-12fghl_CMSSW532.root' % (flavorString), 1355914, 17.54), ]

    #### Cross check using ScaleUp, ScaleDown and Sherpa samples as default ####
    ##set the scale factors such that the expected yield matches the deafult
##     crosscheck_wpjSU_kfactor=117078.9/93123.7;
##     pars.WpJFiles = [(pars.MCDirectory + 'RD_%s_WpJscaleup_CMSSW532.root' % (flavorString), 20784694, crosscheck_wpjSU_kfactor*36257.2), ]
##     crosscheck_wpjSD_kfactor=117078.9/199303.4;
##     pars.WpJFiles = [(pars.MCDirectory + 'RD_%s_WpJscaledown_CMSSW532.root' % (flavorString), 20760830, crosscheck_wpjSD_kfactor*36257.2), ] 
##     crosscheck_wpjSherpa_kfactor=1.0; ## the expected event count should match the default
##     if pars.btagSelection:
##         if isElectron:
##             crosscheck_wpjSherpa_kfactor=6247.1/2032.6
##         else:
##             crosscheck_wpjSherpa_kfactor=6884.5/2404.1
##     else:
##         if isElectron:
##             crosscheck_wpjSherpa_kfactor=102229.7/28765.4
##         else:
##             crosscheck_wpjSherpa_kfactor=118256.8/33448.1
        
##     pars.WpJFiles = [(pars.MCDirectory + 'RD_%s_WJets_sherpa_CMSSW532.root' % (flavorString), 93158078, crosscheck_wpjSherpa_kfactor*36257.2), ] 
    

    # To implement Template Morphing set pars.WpJModels=[-2]
    # and be sure to edit the WpJ*InputParameters.txt file so that the
    # naming scheme corresponds to the correct components/subcomponents.
    # E.g. the parameters from the shape fit to WpJ default MC should
    # contain the suffix Nom, while the overall yield shouldn't, and the
    # lines
    # fMU_WpJ = 0.0 +/- 100.0 L(-1 - 1)
    # fSU_WpJ = 0.0 +/- 100.0 L(-1 - 1)
    # should be added to the .txt file

    if pars.btagSelection:
        pars.WpJModels = [14]
        pars.WpJModelsAlt = [10]
        #pars.WpJAuxModelsAlt = [3]
        if isElectron:
            #pars.WpJFracOfData = 0.2098 #medium
            #pars.WpJFracOfData = 0.51145 #loose
            #pars.WpJFracOfData = 0.4831 #loose with QCD
            #pars.WpJFracOfData = ; #(6247.1+170)/(6247.1+170+124.7+6180.2+20.1) - Since QCD is absorbed into WJets, reset the expected yield to be #VJets+#QCD events
            pars.WpJFracOfData = -1
        else:
            #pars.WpJFracOfData = 0.1935 #medium
            #pars.WpJFracOfData = 0.48379 #loose
            pars.WpJFracOfData = -1
    else:
        pars.WpJModels = [10]
        pars.WpJModelsAlt = [24]
        pars.WpJAuxModelsAlt = [3]
        #pars.WpJModels = [-2]
        #pars.WpJFracOfData = 0.78 
        pars.WpJFracOfData = -1

    #pars.WpJModels = [-2]
    pars.WpJ_NomFiles = pars.WpJFiles
#    pars.WpJ_NomModels = [23]
    pars.WpJ_NomModels = [-1]
    pars.WpJ_NomAuxModels = [5]
    pars.WpJ_MUFiles = [ (pars.MCDirectory + 'RD_%s_WpJmatchingup_CMSSW532.root' % (flavorString), 20976007, 36257.2), (pars.SecondaryDirectory + 'RD_%s_WpJmatchingup_CMSSW532_UserGen.root' % (flavorString), 20976007, 36257.2), ]##the centrally produced and user generated samples are added with equal weights
    pars.WpJ_MUModels = [-1]
    pars.WpJ_MDFiles = [ (pars.MCDirectory + 'RD_%s_WpJmatchingdown_CMSSW532.root' % (flavorString), 21364575, 36257.2), (pars.SecondaryDirectory + 'RD_%s_WpJmatchingdown_CMSSW532_UserGen.root' % (flavorString), 21364575, 36257.2), ]
    pars.WpJ_MDModels = [-1]
    pars.WpJ_SUFiles = [ (pars.MCDirectory + 'RD_%s_WpJscaleup_CMSSW532.root' % (flavorString), 20784694, 36257.2), (pars.SecondaryDirectory + 'RD_%s_WpJscaleup_CMSSW532_UserGen.root' % (flavorString), 20784694, 36257.2), ]
    pars.WpJ_SUModels = [-1]
    pars.WpJ_SDFiles = [ (pars.MCDirectory + 'RD_%s_WpJscaledown_CMSSW532.root' % (flavorString), 20760830, 36257.2), (pars.SecondaryDirectory + 'RD_%s_WpJscaledown_CMSSW532_UserGen.root' % (flavorString), 20760830, 36257.2), ]
    pars.WpJ_SDModels = [-1]

    #Fitting the alternate samples as a cross check (& to potentially inflate the fit errors)
    pars.WpJ_MUModels = [10]
    pars.WpJ_MDModels = [10]
    pars.WpJ_SUModels = [10]
    pars.WpJ_SDModels = [10]
    

    pars.WpJ_NomModelsAlt = pars.WpJ_NomModels
    pars.WpJ_NomAuxModelsAlt = pars.WpJ_NomAuxModels
    pars.WpJ_MUModelsAlt = pars.WpJ_MUModels
    pars.WpJ_MDModelsAlt = pars.WpJ_MDModels
    pars.WpJ_SUModelsAlt = pars.WpJ_SUModels
    pars.WpJ_SDModelsAlt = pars.WpJ_SDModels
    pars.WpJConvModels = pars.GlobalConvModels
    pars.WpJConvModelsAlt = pars.WpJConvModels
    
    if pars.useTopSideband:
        pars.WpJModels = [-1]

#####################  ZpJ: #######################################
    pars.ZpJFiles = [
        (pars.MCDirectory + 'RD_%s_ZpJ_CMSSW532.root' % (flavorString),
         30209426, 3503.71),
        ]

    pars.ZpJFracOfData = -1
    #pars.ZpJFracOfData = 0.045
    pars.ZpJModels = [14]
    #pars.ZpJModels = [-1]
    #pars.ZpJAuxModels = [4]
    if pars.btagSelection:
        pars.ZpJModels = [0]
        pars.ZpJAuxModels = [3]
        pars.ZpJFracOfData = -1

    pars.ZpJModelsAlt = pars.ZpJModels
    pars.ZpJConvModels = pars.GlobalConvModels
    pars.ZpJConvModelsAlt = pars.ZpJConvModels

#####################  top: #######################################    
    pars.topFiles = [
        (pars.MCDirectory + 'RD_%s_TTbar_CMSSW532.root' % (flavorString),
         6893735, 246.73),
##     pars.topFiles = [
##         (pars.MCDirectory + 'RD_%s_TTbar_CMSSW532.root' % (flavorString),
##          6893735, 225.197),
##         (pars.MCDirectory + 'RD_%s_TT_mcatnlo_CMSSW532.root' % (flavorString),
##          32852589, 225.197),
        (pars.MCDirectory + 'RD_%s_STopS_Tbar_CMSSW532.root' % (flavorString),
         139974, 1.75776),
        (pars.MCDirectory + 'RD_%s_STopS_T_CMSSW532.root' % (flavorString),
         259960, 3.89394),
        (pars.MCDirectory + 'RD_%s_STopT_Tbar_CMSSW532.root' % (flavorString),
         1935066, 30.0042),
        (pars.MCDirectory + 'RD_%s_STopT_T_CMSSW532.root' % (flavorString),
         3758221, 55.531),
        (pars.MCDirectory + 'RD_%s_STopTW_Tbar_CMSSW532.root' % (flavorString),
         493458, 11.1773),
        (pars.MCDirectory + 'RD_%s_STopTW_T_CMSSW532.root' % (flavorString),
         497657, 11.1773),
        ]

##     pars.topFiles = [
##         (pars.MCDirectory + 'RD_WmunuJets_DataAll_GoldenJSON_19p3invfb.root' , 1.0, 1.0 ),
##         ]


    if pars.btagSelection:
        if isElectron:
            #pars.topFracOfData = 0.7722 #medium
            #pars.topFracOfData = 0.47349 #loose
            #pars.topFracOfData = 0.4472 #loose with QCD
            pars.topFracOfData = -1
        else:
            #pars.topFracOfData = 0.7914 #medium
            #pars.topFracOfData = 0.50073 #loose
            pars.topFracOfData = -1
        pars.topModels = [23]
        pars.topAuxModels = [3]
    else:
        pars.topFracOfData = -1
        #pars.topFracOfData = 0.046
        pars.topModels = [13] #anti-btag selection
        pars.topAuxModels = [3]
        
    #pars.topModels = [-1]
    pars.topModelsAlt = pars.topModels
    pars.topAuxModelsAlt = pars.topAuxModels
    pars.topConvModels = pars.GlobalConvModels
    pars.topConvModelsAlt = pars.topConvModels

    if pars.useTopSideband:
        pars.topModels = [13] ## Compare the top shapes (between ttbar default vs mcatlo) in the control region
##        pars.topModels = [-1]

#####################  WHbb: #######################################
    pars.WHbbFiles = [
        (pars.MCDirectory + 'RD_%s_WH_WToLNu_HToBB_M-125_CMSSW532.root' % (flavorString), 999998, 0.6966*0.577*(0.1075+0.1057+0.1125)),
        ]

    pars.WHbbFracOfData = -1
    pars.WHbbModels = [5]

    pars.WHbbModelsAlt = pars.WHbbModels
    pars.WHbbConvModels = pars.GlobalConvModels
    pars.WHbbConvModelsAlt = pars.WHbbConvModels
    
#####################  HWW: #######################################
    pars.HWWFiles = [
        (pars.MCDirectory + 'RD_%s_ggh125_CMSSW532.root' % (flavorString), 1966144, 19.27*0.04735),
        (pars.MCDirectory + 'RD_%s_qqh125_CMSSW532.root' % (flavorString), 1140336, 1.578*0.04735),
        ]

    pars.HWWFracOfData = -1
    pars.HWWModels = [5]
    ##pars.HWWModels = [-1]

    pars.HWWModelsAlt = pars.HWWModels
    pars.HWWConvModels = pars.GlobalConvModels
    pars.HWWConvModelsAlt = pars.HWWConvModels

#####################  WWTau backgrounds: #######################################
    #### W->e/mu v and the other W->tau vtau - include in the signal
    pars.WWTauFiles = [
        (pars.MCDirectory + 'RD_%s_WW_lvtauvtau_CMSSW532.root' % (flavorString), 24898+24898, 0.0246772*2*55.9),
        ]

    pars.WWTauFracOfData = -1
    pars.WWTauModels = [5]

    pars.WWTauModelsAlt = pars.WWTauModels
    pars.WWTauConvModels = pars.GlobalConvModels
    pars.WWTauConvModelsAlt = pars.WWTauConvModels

#####################  WZTau backgrounds : #######################################
    #### W->lv and Z->tautau - include in the signal
    pars.WZTauFiles = [
        (pars.MCDirectory + 'RD_%s_WZ_lvtautau_CMSSW532.root' % (flavorString), 24897+24899, 0.0114347*23.6),
        ]

    pars.WZTauFracOfData = -1
    pars.WZTauModels = [5]

    pars.WZTauModelsAlt = pars.WZTauModels
    pars.WZTauConvModels = pars.GlobalConvModels
    pars.WZTauConvModelsAlt = pars.WZTauConvModels



###################################################################
    pars.WZPlotting = {'color' : kGreen+3, 'title' : 'WZ'}
    pars.dibosonPlotting = {'color' : kAzure+8, 'title' : 'WW+WZ'}
    pars.WpJPlotting = { 'color' : kRed, 'title' : 'V+jets'}
    pars.ZpJPlotting = { 'color' : kYellow, 'title' : 'Z+jets'}
    pars.topPlotting = {'color' : kGreen+2, 'title' : 'top'}
    pars.QCDPlotting = {'color' : kGray, 'title' : 'multijet'}
    pars.WHbbPlotting = {'color' : kBlue, 'title' : 'WHbb'}
    pars.HWWPlotting = {'color' : kBlue, 'title' : 'HWW'}
    
    pars.var = ['Mass2j_PFCor']
    pars.varRanges = {'Mass2j_PFCor': (14, 48., 160., []),}
    pars.sigRegionMin = 70.0
    pars.sigRegionMax = 100.0

    if pars.btagSelection:
        pars.varRanges = {'Mass2j_PFCor': (15, 40., 160., []),
                          }
    pars.varTitles = {'Mass2j_PFCor': 'm_{jj}',
                      }
    pars.varNames = {'Mass2j_PFCor': 'Mass2j_PFCor' }


    if pars.useTopSideband:
        pars.var = ['TopWm']
        pars.varRanges = {'TopWm': (14, 48., 160., []),}
        pars.varTitles = {'TopWm': 'Top m_{jj}',}
        pars.varNames = {'TopWm': 'TopWm' }
    
    
    pars.exclude = {}
    pars.blind = False
    # pars.v1binEdges = [50, 55.,60.,65.,70.,75.,80.,85.,95.,
    #                    105.,115.,125.,135.,150.,165.,180.,200.]
    # pars.v1nbins = len(pars.v1binEdges)-1

    pars.integratedLumi = 19300.

    # pars.binData = False
    pars.binData = False

    #Standard vs QCD cuts:
    pars.QCDcuts = pars.cuts
    pars.cuts += '&&(event_met_pfmet>25)'
##    pars.QCDcuts += '&&(event_met_pfmet>20)'
    pars.QCDcuts += '&&(event_met_pfmet>20)&&(W_electron_pfIsoEA>0.3)'


    pars.includeSignal = includeSignal
    pars.signals = []

    return customizeElectrons(pars) if isElectron else \
        customizeMuons(pars)

    
        
def customizeElectrons(pars):
    
    pars.DataFile = pars.MCDirectory + 'RD_WenuJets_DataAllSingleElectronTrigger_GoldenJSON_19p2invfb.root'

    if ((not pars.useTopSideband) and (not pars.btagSelection)):
        pars.QCDFiles = [(pars.SecondaryDirectory + 'RDQCD_WenuJets_Isog0p3NoElMVA_19p2invfb.root',
                          1,1), #The events come from the data sideband
                         ]
        pars.backgrounds.append('QCD')
        
##         if pars.btagSelection:
##             pars.QCDFracOfData = 0.054
##             pars.QCDModels = [0]
##             pars.yieldConstraints['QCD'] = 0.5
        pars.QCDFracOfData = 0.07
        pars.QCDModels = [17]
##        pars.QCDModels = [0]
        pars.yieldConstraints['QCD'] = 0.5
##         ### Cross-Check
##         pars.yieldConstraints['QCD'] = 10.0
            
        pars.QCDModelsAlt = pars.QCDModels
        pars.QCDConvModels = pars.GlobalConvModels
        pars.QCDConvModelsAlt = pars.QCDConvModels


    pars.doEffCorrections = True
    pars.effToDo = ['lepton']
    pars.leptonEffFiles = {
        'id': ["EffTable2012/scaleFactor-Run2012ABC-GsfElectronToId.txt"],
        'reco': ["EffTable2012/scaleFactor-Run2012ABC-SCToElectron.txt"],
        'HLT': ["EffTable2012/efficiency-Run2012ABC-WP80ToHLTEle.txt"]
        }
    pars.lumiPerEpoch = [pars.integratedLumi]

    pars.integratedLumi = 19200.
    pars.QCDcuts += '&&(W_electron_pt>30)'
    pars.cuts += '&&(W_electron_pt>30)'
    return pars

def customizeMuons(pars):

    pars.DataFile = pars.MCDirectory + 'RD_WmunuJets_DataAll_GoldenJSON_19p3invfb.root'
    
    pars.doEffCorrections = True
    pars.effToDo = ['lepton']
    pars.leptonEffFiles = {
        'id': ["EffTable2012/scaleFactor-Run2012ABC-RecoToIso.txt"],
        'HLT': ["EffTable2012/efficiency-Run2012ABC-IsoToIsoMuHLT.txt"]
        }
    pars.lumiPerEpoch = [pars.integratedLumi]

    pars.QCDcuts += '&&(abs(W_muon_eta)<2.1)&&(W_muon_pt>25.)'
    pars.cuts += '&&(abs(W_muon_eta)<2.1)&&(W_muon_pt>25.)'
    
    return pars
