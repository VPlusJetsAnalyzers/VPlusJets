#! /usr/bin/env python

import subprocess
import sys

def runCommand(cmd, bn, interf = 0, reuse = False):
    if interf == 2:
        bn += '_up'
        cmd += [ '--interference', str(interf) ]
    elif interf == 3:
        bn += '_down'
        cmd += [ '--interference', str(interf) ]
    elif interf == 1:
        cmd += [ '--interference', str(interf) ]

    cmd += ['--bn', bn, bn+'.txt']
    if reuse:
        cmd += ['--ws', bn+'.root']
    print ' '.join(cmd)
    sys.stdout.flush()
    subprocess.call(cmd)

from optparse import OptionParser

parser = OptionParser()
parser.add_option('-H', '--mH', dest='mH', default=350, type='int',
                  help='Higgs Mass Point')
parser.add_option('--electrons', dest='isElectron', action='store_true',
                  default=False, help='do electrons instead of muons')
parser.add_option('--reuse', dest='reuse', action='store_true',
                  default=False, help='use previous fit to get data')
parser.add_option('--mva', dest='mvaCut', type='float',
                  help='override cut value for mva')
parser.add_option('--sideband', dest='sb', type='int',
                  help='use sideband dataset and model instead')
parser.add_option('--xrootd', dest='xrootd', action='store_true',
                  help='use xrootd file opening.')
parser.add_option('--altConfig',
                  dest='altConfig',
                  help='which config to select look at HWW2DConfig.py for ' +\
                      'an example.  Use the file name minus the .py extension.'
                  )

(opts, args) = parser.parse_args()

components = [ 'ggH', 'qqH', 'diboson', 'top', 'WpJ' ]
if len(args) > 0:
    components = args

mjjConfig = 'HWW1D2FitsConfig'
if opts.altConfig:
    mjjConfig = opts.altConfig

cmdBase = [ 'python', 'fitComponentShapePdf.py', '-b', '-j', '2', 
            '--nativeBins',
            '-m', mjjConfig, '--mH', str(opts.mH) ]

if opts.isElectron:
    cmdBase.append('--electrons')

if opts.mvaCut:
    cmdBase.extend(['--mva', str(opts.mvaCut)])

if opts.sb:
    cmdBase.extend(['--sideband', str(opts.sb)])

if opts.xrootd:
    cmdBase.append('--xrootd')

for component in components:
    cmd = list(cmdBase) + [ '--comp', component, ]
    if component == 'WpJ':
        cmd.append('--makeFree')
    cmd2 = list(cmd)
    cmd2[cmd2.index(mjjConfig)] = 'HWW1D2FitsConfig_mWW'
    cmd2.extend(['--mjj', mjjConfig])
    bn = '1D2FitsParameters/%sHWW%iParameters' % (component, opts.mH)
    if opts.sb:
        bn += '_sideband%i' % opts.sb
    if opts.isElectron:
        bn += '_el'
    if (component == 'ggH') and (opts.mH >= 400):
        for i in range(1, 4):
            runCommand(list(cmd2), bn + '_mWW', interf = i, reuse = opts.reuse)
    else:
        runCommand(cmd2, bn + '_mWW', reuse = opts.reuse)
    runCommand(cmd, bn + '_mjj', reuse = opts.reuse)
