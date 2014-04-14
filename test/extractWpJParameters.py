from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False,
                  help='no X11 windows')
parser.add_option('--degrade', dest='degrade', type='float', default=1.0,
                  help='factor to multipy errors.')
parser.add_option('--fixToZero', dest='fixToZero', action='store_true',
                  help='fixLowParametersToZero')
(opts, args) = parser.parse_args()

import pyroot_logon
import ROOT as r

for fname in args:
    f = r.TFile.Open(fname)
    w_mWW = f.Get('w_mWW')
    out_fname = fname.replace('.root', '.txt')
    pars = w_mWW.pdf('WpJ_extended').getParameters(w_mWW.set('obsSet'))
    pIter = pars.createIterator()
    p = pIter.Next()
    lines = []
    while p:
        if p.GetName().find('_nrm') < 0:
            if not p.isConstant() and (abs(p.getVal()*10) < p.getError()):
                print p.GetName(), 'is not significantly different from zero.',
                print 'Consider fixing it to zero.'
                if opts.fixToZero:
                    print 'fixing',p.GetName(),'to zero'
                    p.setVal(0.0)
                    p.setConstant()
            lines.append("%s = %f +/- %f %s" % \
                             (p.GetName()[:-3], p.getVal(), 
                              p.getError()*opts.degrade, 
                              'C' if p.isConstant() else '')
                         )

        p = pIter.Next()
    print '\n'.join(lines)
    out_f = open(out_fname, 'w')
    lines.append('')
    out_f.write('\n'.join(lines))
    print 'written to',out_fname
    out_f.close()
    pars.IsA().Destructor(pars)
    f.Close()
