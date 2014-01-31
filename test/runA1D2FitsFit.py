#! /usr/bin/env python

import subprocess

def runCommand(cmd):
    print ' '.join(cmd)
    subprocess.call(cmd)

from optparse import OptionParser

parser = OptionParser()
parser.add_option('-H', '--mH', dest='mH', default=350, type='int',
                  help='Higgs Mass Point')
parser.add_option('--electrons', dest='isElectron', action='store_true',
                  default=False, help='do electrons instead of muons')
parser.add_option('--reuse', dest='reuse', help='use previous fit to get data')
parser.add_option('--sig', dest='sig', action='store_true',
                  default=False, help='include signal hypothesis in fit')
parser.add_option('--expLimit', type='int', default=0, dest='limit',
                  help='number of toys for expected limit, 0 is off')
parser.add_option('--mva', dest='mvaCut', type='float',
                  help='override cut value for mva')
parser.add_option('--sideband', dest='sb', type='int',
                  help='use sideband dataset and model instead')
parser.add_option('--xrootd', dest='xrootd', action='store_true',
                  help='use xrootd file opening.')
parser.add_option('--obsLimit', dest='obsLimit', action='store_true',
                  default=False,
                  help='calculate observed limit too')
parser.add_option('--toy', dest='toy', action='store_true',
                  help='use pseudo-data instead of data file')
parser.add_option('--toyOut', dest='toyOut', help='filename for toy output')
parser.add_option('--injectS', type='float', dest='sigInject',
                  help='amount of signal to inject')
parser.add_option('--genConfig',
                  dest='genConfig',
                  help='which config to select look at HWW2DConfig.py for ' +\
                      'an example.  Use the file name minus the .py extension.'
                  )
parser.add_option('--xc', dest='xc', action='store_true',
                  help='use cross-check background to generate')

(opts, args) = parser.parse_args()

import os
import fnmatch

commonCmd = [ 'python', 'runHWW1D2FitsFitter.py', '-b', '-j', str(2), 
              '--mH', str(opts.mH)] #, '--obsLimit']

if opts.limit:
    commonCmd += ['--limit', str(opts.limit)]

if opts.mvaCut:
    commonCmd.extend(['--mva', str(opts.mvaCut)])

if opts.sb:
    commonCmd.extend(['--sideband', str(opts.sb)])
    commonCmd[1] = 'runHWW1D2FitsSideband.py'

if opts.xrootd:
    commonCmd.append('--xrootd')

if opts.obsLimit:
    commonCmd.append('--obsLimit')

if opts.sigInject != None:
    commonCmd.extend(['--injectS', str(opts.sigInject)])

if opts.toy:
    commonCmd.append('--toy')
    commonCmd.extend(['--seed', '0'])

if opts.toyOut:
    print "toyOut:",opts.toyOut
    commonCmd.extend(['--toyOut', opts.toyOut])

if opts.reuse != None:
    flavor = 'electron' if opts.isElectron else 'muon'
    wsname = 'HWW%ilnujj_%s_2jets_1D2Fit_output.root' % (opts.mH, flavor)
    if len(opts.reuse) > 0:
        wsname = opts.reuse
    commonCmd += [ '--ws', wsname ]

if opts.genConfig:
    commonCmd.extend(['--genConfig', opts.genConfig])

if opts.xc:
    commonCmd.append('--xc')

searchString = '*HWW%iParameters' % opts.mH
if opts.sb:
    searchString += '_sideband%i' % opts.sb
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
#print inputPars
muCmd = list(commonCmd) + inputPars

if opts.mH >= 400:
    muCmd += [ '1D2FitsParameters/' + searchString.replace('*', 'ggH').replace('m??', 'mWW_up'),
               '1D2FitsParameters/' + searchString.replace('*', 'ggH').replace('m??', 'mWW_down')
               ]

runCommand(muCmd)
