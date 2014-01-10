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

(opts, args) = parser.parse_args()

cmd = ['python', 'runA1D2FitsFit.py', '--mH', str(opts.mH), '--expLimit', '20',
       '--obsLimit', '--toy', '--injectS', str(opts.sigInject)]

for t in range(opts.start, opts.end):
    tcmd = list(cmd) + ['--toyOut', 'toyFitResults_HWW%i_Sig%.1f_%i.asc' % \
                            (opts.mH,opts.sigInject,t)]
    print ' '.join(tcmd)
    subprocess.call(tcmd)
