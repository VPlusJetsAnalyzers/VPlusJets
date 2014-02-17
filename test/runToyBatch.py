#! /usr/bin/env python

import subprocess

from optparse import OptionParser

parser = OptionParser()
parser.add_option('-H', '--mH', dest='mH', default=300, type='int',
                  help='Higgs Mass Point')
parser.add_option('-s', '--start', dest='start', default=0, type='int',
                  help='starting sequence')
parser.add_option('-e', '--end', dest='end', default=50, type='int',
                  help='ending sequence')
parser.add_option('--injectS', type='float', dest='sigInject', default=0.,
                  help='amount of signal to inject')
parser.add_option('--xc', dest='xc',
                  help='use cross-check background to generate')
parser.add_option('--genConfig', dest='genConfig', help='use this config to'+\
                      'to generate with')
parser.add_option('-m', '--mode', 
                  dest='modeConfig',
                  help='which config to select look at HWW2DConfig.py for ' +\
                      'an example.  Use the file name minus the .py extension.'
                  )

(opts, args) = parser.parse_args()

cmd = ['python', 'runA1D2FitsFit.py', '--mH', str(opts.mH), #'--expLimit', '20',
       '--toy', '--injectS', str(opts.sigInject)]

logfilename = 'HWW%i_toyFit%.1fS_%i_%i.log' % (opts.mH, opts.sigInject, 
                                               opts.start, opts.end)
logfile = open(logfilename, 'w')

print "logfilename:",logfilename

if opts.modeConfig:
    cmd.extend(['-m', opts.modeConfig])

if opts.genConfig:
    cmd.extend(['--genConfig', opts.genConfig])

if opts.xc:
    cmd.extend(['--xc', opts.xc])

for t in range(opts.start, opts.end):
    # if opts.gen_config:
    #     xccmd = ['python', 'genA1D2FitsToy.py',  
    #              '-m', opts.gen_config, '--injectS', str(opts.sigInject),
    #              '--mH', str(opts.mH)]
    #     if opts.xc:
    #         xccmd.append('--xc')
    #     print ' '.join(xccmd)
    #     subprocess.call(xccmd, stdout=logfile, stderr=subprocess.STDOUT)
    #     cmd.extend([ '--reuse', 'HWW%ilnujj_muon_2jets_1D2Fit_gen.root' %\
    #                      (opts.mH) ])
    tcmd = list(cmd) + ['--toyOut', 'toyFitResults_HWW%i_Sig%.1f_%i.asc' % \
                            (opts.mH,opts.sigInject,t)]
    print ' '.join(tcmd)
    subprocess.call(tcmd, stdout=logfile, stderr=subprocess.STDOUT)
