import pyroot_logon
import ROOT as r

def morphSignalShape(basis_m, morph_m, target_m, Cprime = None,
                     BRnew = None, inputDirectory = '.', 
                     outputDirectory = '.', showPlots = False):
    cprime_tag = ''
    if Cprime != None:
        cprime_tag = '_Cprime_%.1f' % Cprime
    brnew_tag = ''
    if BRnew != None:
        brnew_tag = '_BRnew_%.1f' % BRnew

    for flavor in ['muon', 'electron']:
        basisFilename = '%s/HWW%dlnujj_%s_1D2Fit_output%s%s.root' % \
            (inputDirectory, basis_m, flavor, cprime_tag, brnew_tag)
        morphFilename = '%s/HWW%dlnujj_%s_1D2Fit_output%s%s.root' % \
            (inputDirectory, morph_m, flavor, cprime_tag, brnew_tag)
        outFilename = '%s/HWW%dlnujj_%s_1D2Fit_output%s%s.root' % \
            (outputDirectory, target_m, flavor, cprime_tag, brnew_tag)

        f_basis = r.TFile.Open(basisFilename)
        f_morph = r.TFile.Open(morphFilename)

        alpha = (target_m - basis_m)/float(morph_m - basis_m)

        print 'for basis:',basis_m,'morphed with:',morph_m,'to target:',target_m,
        print 'alpha',alpha

        new_ws = morphShapes(f_basis.Get('w_mWW'), f_morph.Get('w_mWW'), alpha)

        # data = new_ws.data("")
        # while data:
        #     # print data
        #     # print "data:",data.GetName(),data.IsA().GetName()
        #     if data.IsA().GetName() == 'RooDataSet':
        #         data.SetName('data_unbinned')
        #     else:
        #         data.SetName('data_obs')
        #     data = new_ws.data("")

        dataList = new_ws.allData()
        for d in dataList:
            if d.IsA().GetName() == 'RooDataSet':
                d.SetName('data_unbinned')
            else:
                d.SetName('data_obs')
            
        # new_ws.Print()
        new_ws.writeToFile(outFilename)
        print outFilename, 'written.'
        f_basis.Close()
        f_morph.Close()
        if showPlots:
            overlayPlots(basisFilename, morphFilename, outFilename)

pdfs = ['ggH_extended', 'qqH_extended', 'ggH_extended_interf_ggHUp',
        'ggH_extended_interf_ggHDown']
def morphShapes(ws_basis, ws_morph, alpha):
    # print ws_basis
    # ws_basis.Print()

    morphSpace = r.RooWorkspace(ws_basis)

    for signal in pdfs:
        if ws_basis.pdf(signal) and ws_morph.pdf(signal):
            parameters = ws_basis.pdf(signal).getParameters(ws_basis.set('obsSet'))
            pIter = parameters.createIterator()
            parameter = pIter.Next()
            while parameter:
                pb = parameter.getVal()
                # print parameter.GetName(), pb,
                try:
                    pm = ws_morph.var(parameter.GetName()).getVal()
                    pnew = pb + (pm-pb)*alpha
                    # print 'and',pm,'->',pnew,
                    morphSpace.var(parameter.GetName()).setVal(pnew)
                except TypeError:
                    morphSpace.var(parameter.GetName()).setVal(pb)
                finally:
                    parameter = pIter.Next()
                    # print

    return morphSpace

def overlayPlots(basisFilename, morphFilename, outFilename):
    f_basis = r.TFile.Open(basisFilename)
    f_morph = r.TFile.Open(morphFilename)
    f_out = r.TFile.Open(outFilename)

    ws_basis = f_basis.Get('w_mWW')
    ws_morph = f_morph.Get('w_mWW')
    ws_out = f_out.Get('w_mWW')

    frame = ws_basis.var('fit_mlvjj').frame()
    colors = [r.kBlue, r.kRed, r.kRed+1, r.kRed-1]
    for (i,pdf) in enumerate(pdfs):
        if ws_out.pdf(pdf):
            ws_out.pdf(pdf).plotOn(frame, r.RooFit.LineColor(colors[i]))
            ws_basis.pdf(pdf).plotOn(frame, r.RooFit.LineColor(colors[i]),
                                     r.RooFit.LineStyle(2))
            ws_morph.pdf(pdf).plotOn(frame, r.RooFit.LineColor(colors[i]),
                                     r.RooFit.LineStyle(4))
    frame.Draw()
    r.gPad.Update()
    r.gPad.WaitPrimitive()
    return frame

        

if __name__ == '__main__':
    from HighMassInterpolation import masses, mcmasses, cprimes, new_BRs
    import bisect
    import sys
    lowpoint = 0

    r.gInterpreter.GenerateDictionary("list< RooAbsData* >",
                                      'list;RooAbsData.h')
    directory = 'HighMassFittingFiles'

    # morphSignalShape(350, 400, 390, Cprime = 1.0,
    #                  BRnew = 0.0, inputDirectory = directory, 
    #                  outputDirectory = '.', showPlots = True)

    # sys.exit(0)


    for mass in masses:
        if not mass in mcmasses:
            for cprime in cprimes:
                for BR_new in new_BRs:
                    low_idx = bisect.bisect(mcmasses, mass, lo=lowpoint)
                    print 'target:',mass,'low_idx',low_idx,
                    low_mass = mcmasses[low_idx-1]
                    high_mass = mcmasses[low_idx]
                    lowpoint = low_idx-1
                    basis_m = high_mass
                    morph_m = low_mass
                    if (mass - low_mass) < (high_mass - mass):
                        basis_m = low_mass
                        morph_m = high_mass
                    if basis_m == 170:
                        basis_m = high_mass
                        morph_m = low_mass
                    if low_mass == 350:
                        basis_m = low_mass
                        morph_m = 400
                    print 'basis:',basis_m,'other:',morph_m
                
                    morphSignalShape(basis_m, morph_m, mass, cprime,
                                     BR_new, directory, 
                                     directory)
