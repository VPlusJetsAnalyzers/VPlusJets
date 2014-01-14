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
parser.add_option('--xc', dest='xc', action='store_true',
                  help='use cross-check background to generate')

(opts, args) = parser.parse_args()

cmd = ['python', 'runA1D2FitsFit.py', '--mH', str(opts.mH), '--expLimit', '20',
       '--obsLimit', '--toy', '--injectS', str(opts.sigInject)]

logfilename = 'HWW%i_toyFit%.1fS_%i_%i.log' % (opts.mH, opts.sigInject, 
                                               opts.start, opts.end)
logfile = open(logfilename, 'w')

print "logfilename:",logfilename

for t in range(opts.start, opts.end):
    if opts.xc:
        xccmd = ['python', 'genA1D2FitsToy.py', '--xc', 
                 '-m', 'HWW1D2FitsConfigPoly', '--injectS', str(opts.sigInject),
                 '--mH', str(opts.mH)]
        print ' '.join(xccmd)
        subprocess.call(xccmd, stdout=logfile, stderr=subprocess.STDOUT)
        cmd.extend([ '--reuse', 'HWW%ilnujj_muon_2jets_1D2Fit_gen.root' %\
                         (opts.mH) ])
    tcmd = list(cmd) + ['--toyOut', 'toyFitResults_HWW%i_Sig%.1f_%i.asc' % \
                            (opts.mH,opts.sigInject,t)]
    print ' '.join(tcmd)
    subprocess.call(tcmd, stdout=logfile, stderr=subprocess.STDOUT)
