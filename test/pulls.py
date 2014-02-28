import ROOT as r

def createPull(theData, curve, curveUp = None, curveDown = None):
    return createResid(theData, curve, curveUp, curveDown, True)
        
def createResid(theData, curve, curveUp = None, curveDown = None,
                normalize = False):
    # from ROOT import TGraphErrors, TMath

    pullGraph = r.TGraphErrors(theData.GetN())

    Npull = 0
    # print 'normalize',normalize

    for datapt in range(0, theData.GetN()):
        binmin = theData.GetX()[datapt] - theData.GetEXlow()[datapt]
        binmax = theData.GetX()[datapt] + theData.GetEXhigh()[datapt]

        binN = theData.GetY()[datapt]*(binmax-binmin)

        curveN = curve.average(binmin, binmax)*(binmax-binmin)

        diffUp = 0.
        diffDown = 0.
        if curveUp and curveDown:
            diffUp = curveUp.average(binmin,binmax)*(binmax-binmin) - curveN
            diffDown = curveN - curveDown.average(binmin,binmax)*(binmax-binmin)

            if diffUp < 0:
                #print diffUp
                hold = diffDown
                diffDown = -1.*diffUp
                diffUp = -1.*diffDown
                
        
        if (binN > curveN):
            errN = r.TMath.Sqrt((theData.GetEYlow()[datapt]*(binmax-binmin))**2 \
                                  + diffUp**2)
        else:
            errN = r.TMath.Sqrt((theData.GetEYhigh()[datapt]*(binmax-binmin))**2 \
                                  + diffDown**2)

        pull = 0.
        flagit=False
        if errN <= 0:
            continue
        if abs(errN - binN) < 0.01:
            # print 'suspicious point at',theData.GetX()[datapt],
            # print 'binN',binN,
            # print 'errN',errN,
            # print 'curveN',curveN
            flagit=True

        if normalize:
            if (errN > 0):
                pull = (binN - curveN)/errN
                if (pull < -4):
                    print 'suspicious point at',theData.GetX()[datapt],
                    print 'binN',binN,
                    print 'errN',errN,
                    print 'curveN',curveN
                    if flagit:
                        continue
                errN = 1.
                
        elif not normalize:
            if (errN > 0):
                pull = (binN - curveN)/(binmax-binmin)
            errN /= (binmax-binmin)

        # print 'bin: (', binmin, ',', binmax, ') N:', binN, '+/-', errN,
        # print 'curve N:', curveN
        # print 'pull:', pull, 'errN:',errN

        pullGraph.SetPoint(Npull, theData.GetX()[datapt], pull)
        pullGraph.SetPointError(Npull, 0., errN)
        Npull += 1

    pullGraph.Set(Npull)
    return pullGraph

def computeChi2(theData, curve):
    pullGraph = createPull(theData, curve)
    chi2 = 0.
    for pnt in range(0, pullGraph.GetN()):
        chi2 += pullGraph.GetY()[pnt]**2

    return (chi2, pullGraph.GetN())

def splitErrCurve(curve):
    # from ROOT import RooCurve

    upperCurve = r.RooCurve()
    upperCurve.SetName('%s_upper' % curve.GetName())

    lowerCurve = r.RooCurve()
    lowerCurve.SetName('%s_lower' % curve.GetName())

    lastx = curve.GetX()[0]
    lasty = curve.GetY()[0]
    inCurve = False
    currentCurve = upperCurve
    for pt in range(0, curve.GetN()):
        x = curve.GetX()[pt]
        y = curve.GetY()[pt]

        if (lastx > x):
            currentCurve = lowerCurve

        if (x == lastx) and inCurve:
            if (lasty != y):
                inCurve = False
        elif not inCurve and (lastx != x):
            inCurve = True
            currentCurve.addPoint(lastx, lasty if (lasty>0.) else 0.)

        if (x != lastx) and inCurve:
            currentCurve.addPoint(x, y if (y>0.) else 0.)

        lastx = x
        lasty = y

    upperCurve.Sort()
    lowerCurve.Sort()
    if lowerCurve.GetY()[0] > upperCurve.GetY()[0]:
        tmp = lowerCurve
        lowerCurve = upperCurve
        upperCurve = tmp
        
    return (upperCurve, lowerCurve)

def curveToHist(curve, hist, debug = False):
    from math import sqrt
    # from ROOT import gPad
    #hist.Sumw2()
    for binx in range(1, hist.GetNbinsX()+1):
        low = hist.GetBinLowEdge(binx)
        high = low + hist.GetBinWidth(binx)
        hist.SetBinContent(binx, curve.average(low,high))
        #hist.SetBinError(binx, sqrt(hist.GetBinContent(binx)))
    if debug:
        print 'drawing background histogram...'
        curve.Draw('al')
        hist.Draw('same')
        r.gPad.Update()
        r.gPad.WaitPrimitive()

    return hist

from array import array

def subtractCurves(curve, mcurve, loX = None, hiX = None, debug = False):

    subCurve = r.RooCurve(curve)
    idxs = array('i', [0]*subCurve.GetN())
    r.TMath.BubbleLow(len(idxs), subCurve.GetX(), idxs)
    if loX == None:
        loX = subCurve.GetX()[idxs[0]]
    if hiX == None:
        hiX = subCurve.GetX()[idxs[-1]]
    print 'range for x:',loX, hiX
    for n in range(0, subCurve.GetN()):
        x = subCurve.GetX()[n]
        y = subCurve.GetY()[n]
        if debug:
            print '(',x,',',y,')', mcurve.interpolate(x),
        y -= mcurve.interpolate(x)
        if (abs(x -loX) <= 1e-6) or (abs(x-hiX) <= 1e-6):
            y = 0.
        if debug:
            print 'new y:',y
        subCurve.SetPoint(n,x,y)

    return subCurve

def clipCurve(curve):
    if (curve.GetX()[0] == curve.GetX()[1]) and (curve.GetY()[0] < 1e-6):
        curve.RemovePoint(0)
    if (curve.GetX()[curve.GetN()-1] == curve.GetX()[curve.GetN()-2]) and \
            (curve.GetY()[curve.GetN()-1] < 1e-6):
        curve.RemovePoint(curve.GetN()-1)

    return curve
