from RooWjj2DFitterPars import Wjj2DFitterPars
from ROOT import kRed, kAzure, kGreen, kBlue, kCyan, kViolet, kGray
import HWWSignalShapes

import HWW1D2FitsConfig

mu2Pars = dict(HWW1D2FitsConfig.mu2Pars)

polyOrders = {
    170 : 6,
    180 : 6,
    190 : 6,
    200 : 6,
    250 : 6,
    300 : 6,
    350 : 6,
    400 : 5,
    450 : 5,
    500 : 5,
    550 : 6,
    600 : 6,
}

for mass in mu2Pars:
    mu2Pars[mass][5]['WpJ'] = (mu2Pars[mass][5]['WpJ'][0], 23, 
                               mu2Pars[mass][5]['WpJ'][2], polyOrders[mass])
el2Pars = dict(mu2Pars)

import HWW1D2FitsConfigNominal

def theConfig(**kwargs):
    return HWW1D2FitsConfigNominal.theConfig(mu2Pars, el2Pars, **kwargs)
