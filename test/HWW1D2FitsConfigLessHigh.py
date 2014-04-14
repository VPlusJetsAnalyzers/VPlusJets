from RooWjj2DFitterPars import Wjj2DFitterPars
from ROOT import kRed, kAzure, kGreen, kBlue, kCyan, kViolet, kGray
import HWWSignalShapes

# dictionaries of tuples keyed on the Higgs mass.  The tuple structure is
# (0: mavVariableName, 1: mvaCut, 2: min4BodyMass, 3: max4BodyMass, 
#  4: n4BodyBins, 5: dictionary with tuples specifying the model for 
#  each of the two components of the fit, 5: dictionary defining sideband 
#  regions)

mu2Pars = {
    170: ( "mva2j170mu", 0.500, 165.0, 245.0, 10,
           {'diboson':(22, 8),'top':(13, 8),'WpJ':(10, 14, None, None),
           # {'diboson':(22, 8),'top':(13, 8),'WpJ':(10, 23, None, 4),
            'ggH':(13, 5),'qqH':(13, 5)}, {"high":(98, 154), "low":(55, 66)} ),
    180: ( "mva2j180mu", 0.500, 165.0, 245.0, 10,
           {'diboson': (22, 12),'top':(13, 8),'WpJ':(10, 8, None, None),
           # {'diboson': (22, 12),'top':(13, 8),'WpJ':(10, 23, None, 4),
            'ggH':(13, 5),'qqH':(13, 5)}, {"high":(98, 154), "low":(55, 66)} ),
    190: ( "mva2j190mu", 0.500, 165.0, 245.0, 10,
           {'diboson': (22, 12), 'top': (13, 8), 'WpJ': (10, 8, None, None), 
           # {'diboson': (22, 12), 'top': (13, 8), 'WpJ': (10, 23, None, 4), 
            'ggH':(13, 5),'qqH':(13, 5)}, {"high":(98, 154), "low":(55, 66)} ),
    200: ( "mva2j200mu", 0.500, 165.0, 245.0, 10,
           {'diboson': (22, 12), 'top': (13, 8), 'WpJ': (10, 14, None, None), 
           # {'diboson': (22, 12), 'top': (13, 8), 'WpJ': (10, 23, None, 4), 
            'ggH':(13, 5),'qqH':(13, 5)}, {"high":(98, 154), "low":(55, 66)} ),
    250: ( "mva2j250mu", 0.500, 200.0, 400.0, 10,
           # {'diboson': (22, 12), 'top': (5, 10), 'WpJ': (10, 18, None, None), 
           {'diboson': (22, 12), 'top': (5, 10), 'WpJ': (10, 14, None, None), 
            'ggH':(13, 5),'qqH':(7, 5)}, {"high":(98, 154), "low":(55, 66)} ),
    300: ( "mva2j300mu", 0.500, 200.0, 380.0, 9,
           {'diboson': (22, 12), 'top': (5, 12), 'WpJ': (10, 14, None, None), 
            'ggH':(13, 5),'qqH':(7, 5)}, {"high":(98, 154), "low":(55, 66)} ),
    350: ( "mva2j350mu", 0.500, 250., 475., 9,
           {'diboson': (22, 12),'top':(13, 12),'WpJ':(10, 37, None, 2),
           # {'diboson': (22, 12),'top':(13, 12),'WpJ':(10, 23, None, 4),
            'ggH':(13, 5),'qqH':(7, 5)}, {"high":(98, 154), "low":(55, 66)} ),
    400: ( "mva2j400mu", 0.500, 300.0, 600, 10,
           {'diboson':(31, 12),'top':(13, 12),'WpJ':(10, 37, None, 2),
            'ggH':(13, 5),'qqH':(7, 5)}, {"high":(98, 154), "low":(55, 66)} ),
    450: ( "mva2j450mu", 0.500, 340.0, 900.0, 14,
           {'diboson':(31, 12),'top':(5, 12),'WpJ':(10, 35, None, None),
           #{'diboson':(31, 12),'top':(5, 12),'WpJ':(10, 23, None, 5),
            'ggH':(13, 5),'qqH':(7, 5)}, {"high":(98, 154), "low":(55, 66)} ),
    500: ( "mva2j500mu", 0.500, 340.0, 900.0, 14,
           {'diboson': (22, 12), 'top': (5, 12), 'WpJ': (10, 35, None, None),
           # {'diboson': (22, 12), 'top': (5, 12), 'WpJ': (10, 23, None, 5),
            'ggH':(13, 5),'qqH':(7, 5)}, {"high":(98, 154), "low":(55, 66)} ),
    550: ( "mva2j550mu", 0.600, 340.0, 900.0, 14,
           {'diboson': (22, 12), 'top': (5, 12), 'WpJ': (10, 35, None, None),
           # {'diboson': (22, 12), 'top': (5, 12), 'WpJ': (10, 23, None, 5),
            'ggH':(13, 5),'qqH':(7, 5)}, {"high":(98, 154), "low":(55, 66)} ),
    600: ( "mva2j600mu", 0.600, 340.0, 860.0, 13,
           {'diboson': (22, 12), 'top': (5, 12), 'WpJ': (10, 33, None, None),
           # {'diboson': (22, 12), 'top': (5, 12), 'WpJ': (10, 23, None, 5),
            'ggH':(13, 5),'qqH':(7, 5)}, {"high":(98, 154), "low":(55, 66)} ),
    }

el2Pars = dict(mu2Pars)

import HWW1D2FitsConfigNominal

def theConfig(**kwargs):
    return HWW1D2FitsConfigNominal.theConfig(mu2Pars, el2Pars, **kwargs)
    
