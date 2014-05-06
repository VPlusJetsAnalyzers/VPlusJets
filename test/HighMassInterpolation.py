import subprocess
import sys

def runCommand(cmd, bn, interf = 0, paramsName = None):
    interfName = ''
    if interf == 2:
        interfName = '_up'
        cmd += [ '--interference', str(interf) ]
    elif interf == 3:
        interfName = '_down'
        cmd += [ '--interference', str(interf) ]
    elif interf == 1:
        cmd += [ '--interference', str(interf) ]
    bn += interfName
    cmd += ['--bn', bn]       
    if paramsName:
        cmd += ['%s%s.txt' % (paramsName,interfName)]
    print ' '.join(cmd)
    sys.stdout.flush()
    if opts.noexec:
        return None
    subprocess.call(cmd)
    return bn + '.txt'

def rewriteParamFile(pfn, suffix):
    if not pfn:
        return pfn
    newLines = []
    with open(pfn, 'r') as pf:
        for l in pf:
            line = l.split()
            line[0] = line[0] + suffix
            newLines.append(' '.join(line))
    with open(pfn, 'w') as pf:
        pf.write('\n'.join(newLines))
    return pfn

masses = range(170,300,5)
masses.extend(range(300,400,10))
masses.extend(range(400,600,20))
mcmasses = [170, 180, 190, 200, 250, 300,
            350, 400, 450, 500, 550, 600]

cprimes = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
new_BRs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
# cprimes = [1.0]
# new_BRs = [0.0]

if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-H', '--mH', dest='mH', default=350, type='int',
                      help='Higgs Mass Point')
    parser.add_option('--electrons', dest='isElectron', action='store_true',
                      default=False, help='do electrons instead of muons')
    parser.add_option('--altConfig',
                      dest='altConfig', default='HWW1D2FitsConfigPoly',
                      help='which config to select look at HWW2DConfig.py ' +\
                          'for an example.  Use the file name minus the ' +\
                          '.py extension.'
                      )
    parser.add_option('--xrootd', dest='xrootd', action='store_true',
                      help='use xrootd file opening.')
    parser.add_option('--noexec', dest='noexec', action='store_true',
                      help='just test don\'t actually run')
    parser.add_option('--inputDir',
                      dest='inputDir', default='HWWPolyFits24Apr2014_mva0.0',
                      help='directory with the nomiminal fit results'
                      )
    parser.add_option('--outputDir',
                      dest='outputDir', default='HighMassFittingFiles',
                      help='directory with the nomiminal fit results'
                      )
    (opts, args) = parser.parse_args()

    import os
    import pyroot_logon
    import ROOT as r

    try:
        os.makedirs(opts.outputDir)
    except OSError:
        print opts.outputDir, 'could not be created.'

    cmd = ['python', 'fitComponentShapePdf.py', '-b', 
           '--mode', 'HWW1D2FitsConfig_mWW', '--mjj', opts.altConfig,
           '--mH', str(opts.mH)]
    flavorTag = ''
    paramTag = '_mu'
    fileTag = 'muon'
    if opts.isElectron:
        cmd.append('--electrons')
        flavorTag = '_el'
        paramTag = flavorTag
        fileTag = 'electron'
    if opts.xrootd:
        cmd.append('--xrootd')

    filename = 'HWW%dlnujj_%s_2jets_1D2Fit_output.root' % (opts.mH, fileTag)
    for cprime in cprimes:
        for BR_new in new_BRs:
           paramFiles = []
           for component in ['ggH', 'qqH']:
               bn = 'HighMassSignalFits/%sHWW%dParameters%s_mWW_Cprime_%.1f_BRnew_%.1f'\
                   % (component, opts.mH, flavorTag, cprime, BR_new)
               pn = '1D2FitsParameters/%sHWW%dParameters%s_mWW' % (component,
                                                                   opts.mH, 
                                                                   flavorTag)
               #pn = None
               subCmd = list(cmd) + ['--comp', component, 
                                     '--Cprime', str(cprime), 
                                     '--BRnew', str(BR_new)]
               if (opts.mH >= 400) and (component == 'ggH'):
                   for i in range(1,4):
                       paramFilename = runCommand(list(subCmd), bn, i, pn)
                       paramFiles.append(paramFilename)
               else:
                   paramFilename = runCommand(subCmd, bn, paramsName = pn)
                   paramFiles.append(paramFilename)

           oldFile = r.TFile.Open('%s/%s' % (opts.inputDir, filename))
           old_ws = oldFile.Get('w_mWW')
           for paramFilename in paramFiles:
               rewriteParamFile(paramFilename, paramTag)
               old_ws.allVars().readFromFile(paramFilename)
           old_ws.writeToFile('%s/HWW%dlnujj_%s_1D2Fit_output_Cprime_%.1f_BRnew_%.1f.root' % (opts.outputDir, opts.mH, fileTag, cprime, BR_new))
