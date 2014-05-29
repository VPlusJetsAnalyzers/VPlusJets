import os
cmssw_base = os.environ['CMSSW_BASE']

from ROOT import gROOT, gSystem, TClass
import subprocess

def RooFitInclude():
    print "adding RooFit ...",
    scramCmd = ['scram','tool','info','roofitcore']
    grepCmd = ['grep', 'INCLUDE']
    pscram = subprocess.Popen(scramCmd, stdout = subprocess.PIPE)
    pgrep = subprocess.Popen(grepCmd, stdin=pscram.stdout,
                             stdout=subprocess.PIPE)
    pscram.stdout.close()
    output = pgrep.communicate()[0]
    if (pgrep.returncode == 0):
        print "done"
        roofitinc = output.split("=")[1].rstrip()
        # print roofitinc
        gROOT.GetInterpreter().AddIncludePath(roofitinc)
        roofitinc = '-I"' + roofitinc + '"'
        # print roofitinc
        gSystem.AddIncludePath(roofitinc)
        return True
    else:
        print "failed"
        print 'scram returned:',pscram.returncode,'grep:',pgrep.returncode
        return False

import pyroot_plain
from pyroot_plain import cmsLabel

if (gSystem.DynamicPathName("libFWCoreFWLite.so",True)):
    gSystem.Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libMMozerpowhegweight.so")
    if os.access(os.path.expandvars("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so"), os.R_OK):
        gSystem.Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so")
    gROOT.GetInterpreter().AddIncludePath(cmssw_base + '/src')
    gSystem.AddIncludePath('-I"' + cmssw_base + '/src"')
    if not RooFitInclude():
        workingdir = os.getcwd()
        print 'changing to', cmssw_base, 'directory'
        os.chdir(cmssw_base)
        RooFitInclude()
        print 'returning to working directory', workingdir
        os.chdir(workingdir)
    gROOT.ProcessLine('.L EffTableReader.cc+')
    gROOT.ProcessLine('.L EffTableLoader.cc+')
    gROOT.ProcessLine('.L CPWeighter.cc+')
    # print 'importing RooFit PDF classes'
    if not TClass.GetClass('RooPowerLaw'):
        gROOT.ProcessLine('.L RooPowerLaw.cc+')
    if not TClass.GetClass('RooPowerFunction'):
        gROOT.ProcessLine('.L RooPowerFunction.cxx+')
    if not TClass.GetClass('RooPowerExpPdf'):
        gROOT.ProcessLine('.L RooPowerExpPdf.cxx+')
    if not TClass.GetClass('RooErfExpPdf'):
        gROOT.ProcessLine('.L RooErfExpPdf.cxx+')
    if not TClass.GetClass('RooErfPdf'):
        gROOT.ProcessLine('.L RooErfPdf.cxx+')
    if not TClass.GetClass('RooTH1DPdf'):
        gROOT.ProcessLine('.L RooTH1DPdf.cxx+')
    if not TClass.GetClass('RooChebyshevPdf'):
        gROOT.ProcessLine('.L RooChebyshevPDF.cc+')
    if not TClass.GetClass('alphaFunction'):
        gROOT.ProcessLine('.L alphaFunction.cxx+')
    if not TClass.GetClass('RooATGCFunction'):
        gROOT.ProcessLine('.L RooATGCFunction.cxx+')
    if not TClass.GetClass('RooExpPoly'):
        gROOT.ProcessLine('.L RooExpPoly.cxx+')
    if not TClass.GetClass('RooRelBWRunningWidth'):
        gROOT.ProcessLine('.L RooRelBWRunningWidth.cxx+')

print 'end of pyroot_logon'

if __name__ == '__main__':
    import ROOT as r
