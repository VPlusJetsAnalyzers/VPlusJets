from RooWjj2DFitterPars import Wjj2DFitterPars
from ROOT import kRed, kAzure, kGreen, kBlue, kCyan, kViolet, kGray
import HWWSignalShapes

import HWW1D2FitsConfig

import copy

mu2Pars = copy.deepcopy(HWW1D2FitsConfig.mu2Pars)

polyOrders = {
    170 : 6,
    180 : 6,
    190 : 6,
    200 : 6,
    250 : 6,
    300 : 6,
    350 : 6,
    400 : 6,
    450 : 7,
    500 : 7,
    550 : 7,
    600 : 7,
}

polyModels = {
    170 : 23,
    180 : 23,
    190 : 23,
    200 : 23,
    250 : 23,
    300 : 23,
    350 : 23,
    400 : 38,
    450 : 38,
    500 : 38,
    550 : 38,
    600 : 38,
}

for mass in mu2Pars:
    mu2Pars[mass][5]['WpJ'] = (mu2Pars[mass][5]['WpJ'][0], polyModels[mass], 
                               mu2Pars[mass][5]['WpJ'][2], polyOrders[mass])
el2Pars = dict(mu2Pars)

import HWW1D2FitsConfigNominal

def theConfig(**kwargs):
    return HWW1D2FitsConfigNominal.theConfig(mu2Pars, el2Pars, **kwargs)
