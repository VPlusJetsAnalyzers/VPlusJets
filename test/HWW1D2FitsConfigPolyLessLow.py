from RooWjj2DFitterPars import Wjj2DFitterPars
from ROOT import kRed, kAzure, kGreen, kBlue, kCyan, kViolet, kGray
import HWWSignalShapes

# dictionaries of tuples keyed on the Higgs mass.  The tuple structure is
# (0: mavVariableName, 1: mvaCut, 2: min4BodyMass, 3: max4BodyMass, 
#  4: n4BodyBins, 5: dictionary with tuples specifying the model for 
#  each of the two components of the fit, 5: dictionary defining sideband 
#  regions)


import HWW1D2FitsConfigLessLow
import HWW1D2FitsConfigPoly

mu2Pars = dict(HWW1D2FitsConfigLessLow.mu2Pars)
for mass in mu2Pars:
    mu2Pars[mass][5]['WpJ'] = HWW1D2FitsConfigPoly.mu2Pars[mass][5]['WpJ']
el2Pars = dict(mu2Pars)

import HWW1D2FitsConfigNominal

def theConfig(**kwargs):
    return HWW1D2FitsConfigNominal.theConfig(mu2Pars, el2Pars, **kwargs)
