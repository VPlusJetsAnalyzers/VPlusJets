#! /usr/bin/env python

import subprocess

def runCommand(cmd):
    print ' '.join(cmd)
    subprocess.call(cmd)

from optparse import OptionParser

parser = OptionParser()
parser.add_option('-H', '--mH', dest='mH', default=300, type='int',
                  help='Higgs Mass Point')
parser.add_option('-m', '--mode', default="HWW1D2FitsConfig", 
                  dest='modeConfig',
                  help='which config to select look at HWW2DConfig.py for ' +\
                      'an example.  Use the file name minus the .py extension.'
                  )
parser.add_option('--electrons', dest='isElectron', action='store_true',
                  default=False, help='do electrons instead of muons')
parser.add_option('--sig', dest='sig', action='store_true',
                  default=False, help='include signal hypothesis in fit')
parser.add_option('--mva', dest='mvaCut', type='float',
                  help='override cut value for mva')
parser.add_option('--xrootd', dest='xrootd', action='store_true',
                  help='use xrootd file opening.')
parser.add_option('--toy', dest='toy', action='store_false', default=True,
                  help='use pseudo-data instead of data file')
parser.add_option('--injectS', type='float', dest='sigInject',
                  help='amount of signal to inject')
parser.add_option('--xc', dest='xc', action='store_true',
                  help='use cross-check background to generate')

(opts, args) = parser.parse_args()

import os
import fnmatch

commonCmd = [ 'python', 'genHWW1D2FitsToy.py', '-b', '-j', str(2), 
              '--mode', opts.modeConfig, '--mH', str(opts.mH)]

if opts.mvaCut:
    commonCmd.extend(['--mva', str(opts.mvaCut)])

if opts.xrootd:
    commonCmd.append('--xrootd')

if opts.sigInject != None:
    commonCmd.extend(['--injectS', str(opts.sigInject)])

if opts.toy:
    commonCmd.extend(['--seed', '0'])

searchString = '*HWW%iParameters' % opts.mH
if opts.isElectron:
    searchString += '_el'
    commonCmd += [ '--electrons' ]
searchString += '_m??.txt'
if opts.sig:
    pass
else:
    commonCmd += [ '--nosig' ]
inputPars = [ '1D2FitsParameters/%s' % n for n in fnmatch.filter(
        os.listdir('1D2FitsParameters/'), searchString ) 
              ]
if opts.xc:
    for (i, fname) in enumerate(inputPars):
        print fname,'->',
        if fname.find('WpJ') >= 0:
            new_fname = fname.replace('mWW', 'mWW_syst_p')
            print new_fname,'->',
            inputPars[i] = new_fname
        print inputPars[i]

#print inputPars
muCmd = list(commonCmd) + inputPars

if opts.mH >= 400:
    muCmd += [ '1D2FitsParameters/' + \
                   searchString.replace('*', 'ggH').replace('m??', 'mWW_up'),
               '1D2FitsParameters/' + \
                   searchString.replace('*', 'ggH').replace('m??', 'mWW_down')
               ]

runCommand(muCmd)
