import pandas as pd
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit
import ROOT
import array
import sys
import os
import subprocess
from ROOT import TFile, gROOT, TGaxis, gStyle


def SetGlobalStyle(**kwargs):
    '''
    Method to set global style.

    Parameters
    ----------

    - padrightmargin (float), default = 0.035
    - padleftmargin (float), default = 0.12
    - padtopmargin (float), default = 0.035
    - padbottommargin (float), default = 0.1

    - titlesize (float), default = 0.050
    - titlesizex (float), default = 0.050
    - titlesizey (float), default = 0.050
    - titlesizez (float), default = 0.050

    - labelsize (float), default = 0.045
    - labelsizex (float), default = 0.045
    - labelsizey (float), default = 0.045
    - labelsizez (float), default = 0.045

    - titleoffset (float), default = 1.2
    - titleoffsetx (float), default = 1.2
    - titleoffsey (float), default = 1.2
    - titleoffsetz (float), default = 1.2

    - opttitle (int), default = 0
    - optstat (int), default = 0

    - padtickx (int), default = 1
    - padticky (int), default = 1

    - maxdigits (int), default no max value

    - palette (int), default kBird
    '''

    # pad margins
    if 'padrightmargin' in kwargs:
        gStyle.SetPadRightMargin(kwargs['padrightmargin'])
    else:
        gStyle.SetPadRightMargin(0.035)

    if 'padleftmargin' in kwargs:
        gStyle.SetPadLeftMargin(kwargs['padleftmargin'])
    else:
        gStyle.SetPadLeftMargin(0.12)

    if 'padtopmargin' in kwargs:
        gStyle.SetPadTopMargin(kwargs['padtopmargin'])
    else:
        gStyle.SetPadTopMargin(0.035)

    if 'padbottommargin' in kwargs:
        gStyle.SetPadBottomMargin(kwargs['padbottommargin'])
    else:
        gStyle.SetPadBottomMargin(0.1)

    # title sizes
    if 'titlesize' in kwargs:
        gStyle.SetTitleSize(kwargs['titlesize'], 'xyz')
    else:
        gStyle.SetTitleSize(0.050, 'xyz')

    if 'titlesizex' in kwargs:
        gStyle.SetTitleSize(kwargs['titlesizex'], 'x')
    if 'titlesizey' in kwargs:
        gStyle.SetTitleSize(kwargs['titlesizex'], 'y')
    if 'titlesizez' in kwargs:
        gStyle.SetTitleSize(kwargs['titlesizex'], 'z')

    # label sizes
    if 'labelsize' in kwargs:
        gStyle.SetLabelSize(kwargs['labelsize'], 'xyz')
    else:
        gStyle.SetLabelSize(0.045, 'xyz')

    if 'labelsizex' in kwargs:
        gStyle.SetLabelSize(kwargs['labelsizex'], 'x')
    if 'labelsizey' in kwargs:
        gStyle.SetLabelSize(kwargs['labelsizey'], 'y')
    if 'labelsizez' in kwargs:
        gStyle.SetLabelSize(kwargs['labelsizez'], 'z')

    # title offsets
    if 'titleoffset' in kwargs:
        gStyle.SetTitleOffset(kwargs['titleoffset'], 'xyz')
    else:
        gStyle.SetTitleOffset(1.2, 'xyz')

    if 'titleoffsetx' in kwargs:
        gStyle.SetTitleOffset(kwargs['titleoffsetx'], 'x')
    if 'titleoffsety' in kwargs:
        gStyle.SetTitleOffset(kwargs['titleoffsety'], 'y')
    if 'titleoffsetz' in kwargs:
        gStyle.SetTitleOffset(kwargs['titleoffsetz'], 'z')

    # other options
    if 'opttitle' in kwargs:
        gStyle.SetOptTitle(kwargs['opttitle'])
    else:
        gStyle.SetOptTitle(0)

    if 'optstat' in kwargs:
        gStyle.SetOptStat(kwargs['optstat'])
    else:
        gStyle.SetOptStat(0)

    if 'padtickx' in kwargs:
        gStyle.SetPadTickX(kwargs['padtickx'])
    else:
        gStyle.SetPadTickX(1)

    if 'padticky' in kwargs:
        gStyle.SetPadTickY(kwargs['padticky'])
    else:
        gStyle.SetPadTickY(1)

    gStyle.SetLegendBorderSize(0)

    if 'maxdigits' in kwargs:
        TGaxis.SetMaxDigits(kwargs['maxdigits'])

    if 'palette' in kwargs:
        gStyle.SetPalette(kwargs['palette'])

    gROOT.ForceStyle()


def SetObjectStyle(obj, **kwargs):
    '''
    Method to set root object style.

    Parameters
    ----------

    - obj: object to set style

    - linecolor (int) default 1 (black)
    - linealpha (float) default 1
    - linewidth (int) default 2
    - linestyle (int) default 1

    - markercolor (int) default 1 (black)
    - markeralpha (float) default 1
    - markerstyle (int) default 20 (full circle)
    - markersize (int) default 20 (full circle)

    - fillcolor (int) default no filling
    - fillalpha (float) default 1
    - fillstyle (int) default 0 (no style)

    - color (int) sets same color for line, marker and fill
    - alpha (float) sets same alpha for line, marker and fill
    '''

    # alpha parameters
    lalpha = kwargs.get('linealpha', 1)
    malpha = kwargs.get('markeralpha', 1)
    falpha = kwargs.get('fillalpha', 1)
    if 'alpha' in kwargs:
        lalpha = kwargs['alpha']
        malpha = kwargs['alpha']
        falpha = kwargs['alpha']
    if 'linealpha' in kwargs:
        lalpha = kwargs['linealpha']
    if 'markeralpha' in kwargs:
        malpha = kwargs['markeralpha']
    if 'fillalpha' in kwargs:
        falpha = kwargs['fillalpha']

    # line styles
    if 'linecolor' in kwargs:
        if lalpha < 1:
            obj.SetLineColorAlpha(kwargs['linecolor'], lalpha)
        else:
            obj.SetLineColor(kwargs['linecolor'])
    else:
        if lalpha < 1:
            obj.SetLineColorAlpha(1, lalpha)
        else:
            obj.SetLineColor(1)

    if 'linewidth' in kwargs:
        obj.SetLineWidth(kwargs['linewidth'])
    else:
        obj.SetLineWidth(2)

    if 'linestyle' in kwargs:
        obj.SetLineStyle(kwargs['linestyle'])
    else:
        obj.SetLineStyle(1)

    # marker styles
    if 'markercolor' in kwargs:
        if malpha < 1:
            obj.SetMarkerColorAlpha(kwargs['markercolor'], malpha)
        else:
            obj.SetMarkerColor(kwargs['markercolor'])
    else:
        if malpha < 1:
            obj.SetMarkerColorAlpha(1, malpha)
        else:
            obj.SetMarkerColor(1)

    if 'markersize' in kwargs:
        obj.SetMarkerSize(kwargs['markersize'])
    else:
        obj.SetMarkerSize(1)

    if 'markerstyle' in kwargs:
        obj.SetMarkerStyle(kwargs['markerstyle'])
    else:
        obj.SetMarkerStyle(20)

    # fill styles
    if 'fillcolor' in kwargs:
        if falpha < 1:
            obj.SetFillColorAlpha(kwargs['fillcolor'], falpha)
        else:
            obj.SetFillColor(kwargs['fillcolor'])

    if 'fillstyle' in kwargs:
        obj.SetFillStyle(kwargs['fillstyle'])

    #global color
    if 'color' in kwargs:
        if lalpha < 1:
            obj.SetLineColorAlpha(kwargs['color'], lalpha)
        else:
            obj.SetLineColor(kwargs['color'])
        if malpha < 1:
            obj.SetMarkerColorAlpha(kwargs['color'], malpha)
        else:
            obj.SetMarkerColor(kwargs['color'])
        if falpha < 1:
            obj.SetFillColorAlpha(kwargs['color'], falpha)
        else:
            obj.SetFillColor(kwargs['color'])


def GetInvMassHistAndFit(infile, ptmin, ptmax, nbin, hasfeflections=False):
    inFile = ROOT.TFile.Open(infile)
    
    cMassVsV2 = inFile.Get(f'cSimFit_Pt{ptmin}_{ptmax}')
    hInvMass = cMassVsV2.GetPad(1).GetListOfPrimitives().FindObject(f'MassForFit{nbin}')
    # hInvMass = cMassVsV2.GetPad(1).GetListOfPrimitives().FindObject(f'MassForFit')
    # hInvMass = cMassVsV2.GetPad(1).GetListOfPrimitives().FindObject(f'hist_mass')
    fMassTot = cMassVsV2.GetPad(1).GetListOfPrimitives().FindObject('fMassTotFunc')
    fMassBkg = cMassVsV2.GetPad(1).GetListOfPrimitives().FindObject('fMassBkgFunc')
    hV2VsMass = cMassVsV2.GetPad(2).GetListOfPrimitives().FindObject('hDummy')
    fV2Tot = cMassVsV2.GetPad(2).GetListOfPrimitives().FindObject('fVnTotFunc')
    fV2Bkg = cMassVsV2.GetPad(2).GetListOfPrimitives().FindObject('fVnBkgFunc')
    if hasfeflections:
        fMassRefl = cMassVsV2.GetPad(1).GetListOfPrimitives().FindObject('fMassRflFunc')    
        SetObjectStyle(fMassRefl, linewidth=1, color=ROOT.kGreen+2, fillalpha=0.3, fillstyle=1001)
    
    SetObjectStyle(hInvMass, color=ROOT.kBlack)
    SetObjectStyle(fMassTot, linewidth=3, linecolor=ROOT.kAzure+4)
    SetObjectStyle(fMassBkg, linewidth=3, linecolor=ROOT.kRed-4, linestyle=2)
    SetObjectStyle(hV2VsMass, color=ROOT.kBlack)
    SetObjectStyle(fV2Tot, linewidth=3, linecolor=ROOT.kAzure+4)
    SetObjectStyle(fV2Bkg, linewidth=3, linecolor=ROOT.kRed-4, linestyle=2)
    
    fMassTot.SetNpx(1000)
    fMassBkg.SetNpx(1000)
    fV2Tot.SetNpx(1000)
    fV2Bkg.SetNpx(1000)

    if hasfeflections:
        return hInvMass, fMassTot, fMassBkg, hV2VsMass, fV2Tot, fV2Bkg, fMassRefl
    else:
        return hInvMass, fMassTot, fMassBkg, hV2VsMass, fV2Tot, fV2Bkg


def GetV2HistAndFit(infile, dir, ptmin, ptmax, nbin, hasfeflections=False):
    inFile = ROOT.TFile.Open(infile)
    if dir:
        cMassVsV2 = inFile.Get(f'{dir}/cFrac_{ptmin}_{ptmax}') 
    else: 
        cMassVsV2 = inFile.Get(f'cFrac_{ptmin}_{ptmax}')
    gv2 = cMassVsV2.GetListOfPrimitives().FindObject('Graph')
    hV2VsFrac = cMassVsV2.GetListOfPrimitives().FindObject(f'hV2VsFrac_{nbin}')
    tf1 = gv2.GetFunction("linear")
    
    SetObjectStyle(gv2, color=ROOT.kBlack, markersize=1.2, linewidth=1)
    SetObjectStyle(hV2VsFrac, linewidth=1, linecolor=ROOT.kAzure+4, markersize=0, alpha=0.2)
    SetObjectStyle(tf1, linewidth=1, linecolor=ROOT.kRed-4, markersize=0, linestyle=9)
    
    return gv2, hV2VsFrac, tf1


def GetCanvas4sub(name, xmins, xmaxs, ymins_mass, ymaxs_mass, ymins_v2, ymaxs_v2, axisnametop, axisnamebottom):
    """
    Creates a canvas with 4 adjacent subpads (one for mass on left top, one for v2 on left bottom
    one for cut_variation on right top, one for v2 vs FD fraction on right bottom),
    sharing the x-axis while maintaining different y-axis ranges.

    Args:
        name (str): Name of the canvas.
        xmins (float): Minimum x-axis value (common for both pads).
        xmaxs (float): Maximum x-axis value (common for both pads).
        ymins_mass (float): Minimum y-axis value for the mass panel.
        ymaxs_mass (float): Maximum y-axis value for the mass panel.
        ymins_v2 (float): Minimum y-axis value for the v2 panel.
        ymaxs_v2 (float): Maximum y-axis value for the v2 panel.
        axisnametop (str): Y-axis title for the mass plot.
        axisnamebottom (str): Y-axis title for the v2 plot.

    Returns:
        tuple: (canvas, frames)
    """
    # Create canvas with 2 rows (top = mass, bottom = v2)
    canvas = ROOT.TCanvas(name, name, 1200, 1100)  
    canvas.Divide(2, 2)  # Small spacing between top and bottom

    frames = []
    for i in range(4):
        canvas.cd(i + 1)
        pad = ROOT.gPad
        # Set the correct y-axis range for each pad
        if i == 0:  # Mass plot (top)
            frame = pad.DrawFrame(xmins, ymins_mass, xmaxs, ymaxs_mass, axisnametop)
        elif i == 1:  
            frame = pad.DrawFrame(0.5, 0, 20.5, 20000, ';Minimum BDT score for prompt #Lambda_{c}^{+}; raw yield')
        elif i == 2:  # v2 plot (bottom)
            frame = pad.DrawFrame(xmins, ymins_v2, xmaxs, ymaxs_v2, axisnamebottom)
        elif i == 3:  # v2 vs Fraction (bottom)
            frame = pad.DrawFrame(0, 0, 1.05, 0.3, ';Non-prompt fraction; #it{v}_{2}^{obs.}{SP, |#Delta#it{#eta}| > 1.3}')
        frame.SetTitle("")

        frames.append(frame)

    return canvas, frames


def GetCanvas2sub(name, xmins, xmaxs, ymins_mass, ymaxs_mass, ymins_v2, ymaxs_v2, axisnameleft, axisnameright):
    """
    Creates a canvas with 2 horizontal subpads (one for mass on left, one for v2 on right),
    sharing the y-axis while maintaining different x-axis ranges.
    
    Args:
        name (str): Name of the canvas.
        xmins (float): Minimum x-axis value for the mass panel.
        xmaxs (float): Maximum x-axis value for the mass panel.
        ymins_mass (float): Minimum y-axis value (common for both pads).
        ymaxs_mass (float): Maximum y-axis value (common for both pads).
        ymins_v2 (float): Minimum x-axis value for the v2 panel.
        ymaxs_v2 (float): Maximum x-axis value for the v2 panel.
        axisnameleft (str): X-axis title for the mass plot.
        axisnameright (str): X-axis title for the v2 plot.

    Returns:
        tuple: (canvas, frames)

    # Create canvas with 2 columns (left = mass, right = v2)
    """
    canvas = ROOT.TCanvas(name, name, 1000, 600)  
    canvas.Divide(2, 1, 0.0, 0.0)
    frames = []
    for i in range(2):
        canvas.cd(i + 1)
        pad = ROOT.gPad
        # Set the correct x-axis range for each pad
        if i == 0:  # Mass plot (left)
            frame = pad.DrawFrame(xmins, ymins_mass, xmaxs, ymaxs_mass, axisnameleft)
            frame.SetTitle("")
            frame.GetYaxis().SetTickLength(0.040)
            frame.GetXaxis().SetTickLength(0.025)
        elif i == 1:  # v2 plot (right)
            frame = pad.DrawFrame(xmins, ymins_v2, xmaxs, ymaxs_v2, axisnameright)
            # set thick y-axis on the right side
            frame.GetYaxis().SetTitleSize(0.8)
            frame.GetYaxis().SetTitleOffset(1.6)
        frame.SetTitle("")

        frames.append(frame)

    return canvas, frames



def GetLegend(xmin=0.19, ymin=0.62, xmax=0.75, ymax=0.77, textsize=0.04, ncolumns=2, header=' ', fillstyle=0):
    """
    Creates a formatted legend.
    
    Args:
        xmin, ymin, xmax, ymax (float): Legend position.
        textsize (float): Text size in legend.
        ncolumns (int): Number of columns in legend.
        header (str): Header text for legend.
        fillstyle (int): Fill style.
    
    Returns:
        ROOT.TLegend: Configured legend.
    """
    leg = ROOT.TLegend(xmin, ymin, xmax, ymax)
    leg.SetTextSize(textsize)
    leg.SetNColumns(ncolumns)
    leg.SetFillStyle(fillstyle)
    leg.SetHeader(header)
    return leg


def DrawLineAt0(min, max, title=False):
    line = ROOT.TLine(min, 0, max, 0)
    line.SetLineStyle(9)
    line.SetLineWidth(2)
    line.SetLineColor(ROOT.kGray+2) 
    return line

def SystX(graph, syst, percentage=0.2):
    for i in range(graph.GetN()):
        stat = 2*graph.GetErrorXlow(i)
        syst.SetPointEXlow(i, stat*percentage)
        syst.SetPointEXhigh(i, stat*percentage)
    return syst

def GetPrediction(file, xmin, xmax, scale, isUnc, model='tamu'):
    if model == 'tamu':
        splinev2, _, _, _ = ReadTAMUv2(file)
        gpred = ROOT.TGraphAsymmErrors(1)
        for ipt in range(xmin, xmax):
            index = ipt
            ipt *= scale

            gpred.AddPoint(ipt, splinev2['yCent'](ipt))
            if isUnc:
                gpred.SetPointError(index, 0., 0.,
                                    float(splinev2['yCent'](ipt) - splinev2['yMin'](ipt)),
                                    float(splinev2['yMax'](ipt) - splinev2['yCent'](ipt)))
    elif model == 'lido':
        splinev2, _, _, _ = ReadLIDOV2(file)
        gpred = ROOT.TGraphAsymmErrors(1)
        for ipt in range(xmin, xmax):
            ipt *= scale

            gpred.AddPoint(ipt, splinev2['yCent'](ipt))
            if isUnc:
                gpred.SetPointError(ipt, 0., 0.,
                                    float(splinev2['yCent'](ipt) - splinev2['yMin'](ipt)),
                                    float(splinev2['yMax'](ipt) - splinev2['yCent'](ipt)))
    elif model == 'lgr':
        splinev2, _, _, _ = ReadLGRV2(file)
        gpred = ROOT.TGraphAsymmErrors(0)
        index=0
        for ipt in range(xmin, xmax):
            ipt *= scale

            gpred.AddPoint(ipt, splinev2['yCent'](ipt))
            if isUnc:
                index+=1
                errlow=splinev2['yMax'](ipt) - splinev2['yCent'](ipt)
                print(errlow)
                gpred.SetPointError(index, 0., 0.,
                                    float(splinev2['yCent'](ipt) - splinev2['yMin'](ipt)),
                                    float(splinev2['yMax'](ipt) - splinev2['yCent'](ipt)))
    elif model == 'phsd':
        splinePHSD3040, splinePHSD4050, splinePHSD6080, _, _, _ = ReadPHSDV2(file)
        gpred3040 = ROOT.TGraphAsymmErrors(0)
        index=0
        for ipt in range(xmin, xmax):
            ipt *= scale

            gpred3040.AddPoint(ipt, splinePHSD3040['yCent'](ipt))
            if isUnc:
                index+=1
                errlow=splinePHSD3040['yMax'](ipt) - splinePHSD3040['yCent'](ipt)
                print(errlow)
                gpred3040.SetPointError(index, 0., 0.,
                                    float(splinePHSD3040['yCent'](ipt) - splinePHSD3040['yMin'](ipt)),
                                    float(splinePHSD3040['yMax'](ipt) - splinePHSD3040['yCent'](ipt)))
        gpred3040.RemovePoint(0)
        gpred3040.RemovePoint(0)

        gpred4050 = ROOT.TGraphAsymmErrors(0)
        index=0
        for ipt in range(xmin, xmax):
            ipt *= scale

            gpred4050.AddPoint(ipt, splinePHSD4050['yCent'](ipt))
            if isUnc:
                index+=1
                errlow=splinePHSD4050['yMax'](ipt) - splinePHSD4050['yCent'](ipt)
                print(errlow)
                gpred4050.SetPointError(index, 0., 0.,
                                    float(splinePHSD4050['yCent'](ipt) - splinePHSD4050['yMin'](ipt)),
                                    float(splinePHSD4050['yMax'](ipt) - splinePHSD4050['yCent'](ipt)))
        gpred4050.RemovePoint(0)
        gpred4050.RemovePoint(0)

        gpred6080 = ROOT.TGraphAsymmErrors(0)
        index=0
        for ipt in range(xmin, xmax):
            ipt *= scale

            gpred6080.AddPoint(ipt, splinePHSD6080['yCent'](ipt))
            if isUnc:
                index+=1
                errlow=splinePHSD6080['yMax'](ipt) - splinePHSD6080['yCent'](ipt)
                print(errlow)
                gpred6080.SetPointError(index, 0., 0.,
                                    float(splinePHSD6080['yCent'](ipt) - splinePHSD6080['yMin'](ipt)),
                                    float(splinePHSD6080['yMax'](ipt) - splinePHSD6080['yCent'](ipt)))
        gpred6080.RemovePoint(0)
        gpred6080.RemovePoint(0)
        
        return gpred3040, gpred4050, gpred6080
    return gpred

def compute_central_graph(graph_upper, graph_lower):
    """
    Given two TGraphAsymmErrors (upper and lower limits), computes the central values.
    Ensures all unique x-values from both graphs are included.

    Parameters:
        graph_upper (ROOT.TGraphAsymmErrors): Graph representing the upper limit.
        graph_lower (ROOT.TGraphAsymmErrors): Graph representing the lower limit.

    Returns:
        ROOT.TGraphAsymmErrors: Graph with central values and computed uncertainties.
    """
    # Extract x and y values from both graphs
    x_upper = np.array([graph_upper.GetPointX(i) for i in range(graph_upper.GetN())])
    print(x_upper)
    y_upper = np.array([graph_upper.GetPointY(i) for i in range(graph_upper.GetN())])

    x_lower = np.array([graph_lower.GetPointX(i) for i in range(graph_lower.GetN())])
    print(x_lower)
    y_lower = np.array([graph_lower.GetPointY(i) for i in range(graph_lower.GetN())])

    # Get all unique x-values from both graphs
    x_common = np.sort(np.unique(np.concatenate((x_upper, x_lower))))

    # Interpolate/extrapolate the y-values for both graphs at these x-values
    y_upper_interp = np.interp(x_common, x_upper, y_upper, left=y_upper[0], right=y_upper[-1])
    y_lower_interp = np.interp(x_common, x_lower, y_lower, left=y_lower[0], right=y_lower[-1])

    # Compute central values
    y_central = 0.5 * (y_upper_interp + y_lower_interp)

    # Compute uncertainties as differences from central values
    y_err_low = np.abs(y_central - y_lower_interp)
    y_err_high = np.abs(y_upper_interp - y_central)
    print(x_common)

    # Create a new TGraphAsymmErrors
    graph_central = ROOT.TGraphAsymmErrors(len(x_common))
    for i, x in enumerate(x_common):
        graph_central.SetPoint(i, x, y_central[i])
        graph_central.SetPointError(i, 
                                    0, 0,  # No x-errors assumed
                                    y_err_low[i], y_err_high[i])

    return graph_central


def ReadPHSDV2(fileName):
    '''
    Method to read PHSD v2 files

    Inputs
    ----------
    - fileName: file name

    Returns
    ----------
    splinePHSD: dictionary of splines with PHSD predictions {yCent, yMin, yMax}
    dfPHSD: pandas dataframe with original values
    ptMin: minimum pt for which the model is valid
    ptMax: maximum pt for which the model is valid
    '''

    dfPHSD = pd.read_csv(fileName, sep=' ', comment='#')
    dfPHSD['3040_min'] = dfPHSD['3040'] - dfPHSD['3040_unc']
    dfPHSD['3040_max'] = dfPHSD['3040'] + dfPHSD['3040_unc']
    dfPHSD['4050_min'] = dfPHSD['4050'] - dfPHSD['4050_unc']
    dfPHSD['4050_max'] = dfPHSD['4050'] + dfPHSD['4050_unc']
    dfPHSD['6080_min'] = dfPHSD['6080'] - dfPHSD['6080_unc']
    dfPHSD['6080_max'] = dfPHSD['6080'] + dfPHSD['6080_unc']

    splinePHSD3040, ptMin, ptMax = InterpolateModel(dfPHSD['pT'], dfPHSD['3040'], dfPHSD['3040_max'], dfPHSD['3040_min'])
    splinePHSD4050, _, _ = InterpolateModel(dfPHSD['pT'], dfPHSD['4050'], dfPHSD['4050_max'], dfPHSD['4050_min'])
    splinePHSD6080, _, _ = InterpolateModel(dfPHSD['pT'], dfPHSD['6080'], dfPHSD['6080_max'], dfPHSD['6080_min'])

    return splinePHSD3040, splinePHSD4050, splinePHSD6080, dfPHSD, ptMin, ptMax


def ReadLIDOV2(fileName):
    '''
    Method to read LIDO Raa files

    Inputs
    ----------
    - fileName: file name

    Returns
    ----------
    splineLIDO: dictionary of splines with LIDO predictions {yCent, yMin, yMax}
    dfLIDO: pandas dataframe with original values
    ptMin: minimum pt for which the model is valid
    ptMax: maximum pt for which the model is valid
    '''

    dfLIDO = pd.read_csv(fileName, sep=' ', comment='#')
    dfLIDO['v2_min'] = dfLIDO['v2'] - dfLIDO['v2-error']
    dfLIDO['v2_max'] = dfLIDO['v2'] + dfLIDO['v2-error']

    splineLIDO, ptMin, ptMax = InterpolateModel(dfLIDO['pT'], dfLIDO['v2'], dfLIDO['v2_max'], dfLIDO['v2_min'])

    return splineLIDO, dfLIDO, ptMin, ptMax

def ReadLGRV2(fileName):
    '''
    Method to read LGR Raa files

    Inputs
    ----------
    - fileName: file name

    Returns
    ----------
    splineLGR: dictionary of splines with LGR predictions {yCent, yMin, yMax}
    dfLGR: pandas dataframe with original values
    ptMin: minimum pt for which the model is valid
    ptMax: maximum pt for which the model is valid
    '''

    dfLGR = pd.read_csv(fileName, sep=' ', comment='#')
    splineLGR, ptMin, ptMax = InterpolateModel(dfLGR['pT'], dfLGR['v2'], dfLGR['v2_min'], dfLGR['v2_max'])

    return splineLGR, dfLGR, ptMin, ptMax

def InterpolateModel(ptCent, yCent, yMin=None, yMax=None):
    '''
    Helper function to interpolate model predictions.
    The returned splines will raise an error if applied out of the data boundary.

    Parameters
    -----------
    ptCent: list of pT centres to interpolate
    yCent: list of central values to interpolate
    yMin: list of min values to interpolate
    yMax: list of max values to interpolate

    Returns:
    -----------
    splinesAll: dictionary with splines {yCent, yMin, yMax}
    ptMin: minimum pt for which the interpolation is valid
    ptMax: maximum pt for which the interpolation is valid
    '''

    splinesAll = {}
    splinesAll['yCent'] = InterpolatedUnivariateSpline(ptCent, yCent, ext='raise', check_finite=True)

    if yMin is not None and yMin.any():
        splinesAll['yMin'] = InterpolatedUnivariateSpline(ptCent, yMin, ext='raise', check_finite=True)
    if yMax is not None and yMax.any():
        splinesAll['yMax'] = InterpolatedUnivariateSpline(ptCent, yMax, ext='raise', check_finite=True)

    return splinesAll, min(ptCent), max(ptCent)


def ReadTAMUv2(fileNameTAMUv2):
    '''
    Helper function to read TAMU v2 txt files

    Parameters
    -----------
    fileNameTAMU: TAMU file name

    Returns:
    -----------
    splineTAMU: dictionary with splines {yCent, yMin, yMax}
    dfTAMU: pandas dataframe with original values
    ptMin: minimum pt for which the model is valid
    ptMax: maximum pt for which the model is valid
    '''
    dfTAMU = pd.read_csv(fileNameTAMUv2, sep=' ', comment='#')
    if 'v2min' in dfTAMU and 'v2max' in dfTAMU:
        dfTAMU['v2'] = (dfTAMU['v2min'] + dfTAMU['v2max']) / 2 #central value taken as average of min and max
        splineTAMU, ptMin, ptMax = InterpolateModel(dfTAMU['pT'], dfTAMU['v2'],
                                                    dfTAMU['v2min'], dfTAMU['v2max'])
    else:
        splineTAMU, ptMin, ptMax = InterpolateModel(dfTAMU['pT'], dfTAMU['v2'])

    return splineTAMU, dfTAMU, ptMin, ptMax 


def get_particle_info(particleName):
    '''
    Get particle information

    Input:
        - particleName: 
            the name of the particle

    Output:
        - particleTit: 
            the title of the particle
        - massAxisTit: 
            the title of the mass axis
        - decay: 
            the decay of the particle
        - massForFit: 
            float, the mass of the particle
    '''

    if particleName == 'Dplus':
        particleTit = 'D^{+}'
        massAxisTit = '#it{M}(K#pi#pi) (GeV/#it{c}^{2})'
        massForFit = ROOT.TDatabasePDG.Instance().GetParticle(411).Mass()
        decay = 'D^{+} #rightarrow K^{#minus}#pi^{+}#pi^{+}'
    elif particleName == 'Ds':
        particleTit = 'D_{s}^{+}'
        massAxisTit = '#it{M}(KK#pi) (GeV/#it{c}^{2})'
        decay = 'D_{s}^{+} #rightarrow #phi#pi^{+} #rightarrow K^{+}K^{#minus}#pi^{+}'
        massForFit = ROOT.TDatabasePDG.Instance().GetParticle(431).Mass()
    elif particleName == 'LctopKpi':
        particleTit = '#Lambda_{c}^{+}'
        massAxisTit = '#it{M}(pK#pi) (GeV/#it{c}^{2})'
        decay = '#Lambda_{c}^{+} #rightarrow pK^{#minus}#pi^{+}'
        massForFit = ROOT.TDatabasePDG.Instance().GetParticle(4122).Mass()
    elif particleName == 'LctopK0s':
        massAxisTit = '#it{M}(pK^{0}_{s}) (GeV/#it{c}^{2})'
        decay = '#Lambda_{c}^{+} #rightarrow pK^{0}_{s}'
        massForFit = 2.25 # please carefully check the mass of Lc->pK0s, it is constant
        # massForFit = ROOT.TDatabasePDG.Instance().GetParticle(4122).Mass()
    elif particleName == 'Dstar':
        particleTit = 'D^{*+}'
        massAxisTit = '#it{M}(K#pi#pi) - #it{M}(K#pi) (GeV/#it{c}^{2})'
        decay = 'D^{*+} #rightarrow D^{0}#pi^{+} #rightarrow K^{#minus}#pi^{+}#pi^{+}'
        massForFit = ROOT.TDatabasePDG.Instance().GetParticle(413).Mass() - ROOT.TDatabasePDG.Instance().GetParticle(421).Mass()
    elif particleName == 'Dzero':
        particleTit = 'D^{0}'
        massAxisTit = '#it{M}(K#pi) (GeV/#it{c}^{2})'
        decay = 'D^{0} #rightarrow K^{#minus}#pi^{+}'
        massForFit = ROOT.TDatabasePDG.Instance().GetParticle(421).Mass()
    else:
        print(f'ERROR: the particle "{particleName}" is not supported! Choose between Dzero, Dplus, Ds, Dstar, and Lc. Exit!')
        sys.exit()

    return particleTit, massAxisTit, decay, massForFit


def read_txt(file_path, sep=",", header=None, nrows=None):
    return pd.read_csv(file_path, sep=sep, header=header, nrows=nrows, engine='python').astype('float64')


def kEt(m, pt):
    '''kEt = sqrt(m^2 + pt^2) - m'''
    if isinstance(pt, list):
        return [np.sqrt(m**2+ipt**2)-m for ipt in pt]
    else:
        return np.sqrt(m**2+pt**2)-m


def nq_scaling(x, nq):
    if isinstance(x, list):
        return [ix/nq for ix in x]
    else:
        return x/nq


def fit(ptCent, yCent, med, getParams=False):
    '''fitting to extend the model to high pt;be used in pt weight; can't trust the model at high pt'''
    # ------------------------------------------------------------
    # 1. Extract descending portion data (starting from peak point)
    # ------------------------------------------------------------
    peak_idx = np.argmax(yCent)  # Find index of maximum value (peak_idx=11)
    x_down = ptCent[peak_idx:]   # Extract x data for descending portion
    y_down = yCent[peak_idx:]    # Extract y data for descending portion

    # ------------------------------------------------------------
    # 2. Define fitting function (exponential decay model)
    # ------------------------------------------------------------
    def power_law(x, A, alpha):
        return A * x**(-alpha)

    def exp_decay(x, A, alpha):
        """
        Exponential decay function: y = A * exp(-alpha * x)
        """
        return A * np.exp(-alpha * x)
    # ------------------------------------------------------------
    # 3. Perform fitting
    # ------------------------------------------------------------
    # Initial parameter guesses (A = maximum value, alpha = decay rate)
    p0 = [y_down[0], 0.001]  # A=0.00361, alpha=0.1

    # Perform nonlinear least squares fitting
    params, cov = curve_fit(exp_decay, x_down, y_down, p0=p0)
    A_fit, alpha_fit = params
    if getParams:
        return A_fit, alpha_fit
    # # Calculate goodness of fit R²
    # y_fit = exp_decay(x_down, A_fit, alpha_fit)
    # residuals = y_down - y_fit
    # ss_res = np.sum(residuals**2)
    # ss_tot = np.sum((y_down - np.mean(y_down))**2)
    # r_squared = 1 - (ss_res / ss_tot)

    # # ------------------------------------------------------------
    # # 4. Generate extrapolated data (extended to pt=25)
    # # ------------------------------------------------------------
    x_extended = np.linspace(med, 24, 100)  # Extrapolate starting from peak point
    y_extended = exp_decay(x_extended, A_fit, alpha_fit)
    return x_extended, y_extended


def preprocess(file_path, do_interp=False, do_fit_extend=False, catania=False, sep=" ", do_ncq=False, header=None, nrows=None):
    '''preprocess the model data, including interpolation (only be used to caculate chi^2) and fitting (can't trust) to extend to high pt'''
    df = read_txt(file_path, sep, header=header, nrows=nrows)
    if do_ncq:
        # df = df[df.iloc[:, 0] <= 24]
        return df
    med = 11
    if not do_interp and not do_fit_extend:
        print(f'use source data in {file_path}')
        if catania:
            return df.iloc[:, 0], df.iloc[:, 1], df.iloc[:, 2]
        return df.iloc[:, 0], df.iloc[:, 1]
    if do_interp:
        # pchip = InterpolatedUnivariateSpline(df.iloc[:, 0], df.iloc[:, 1]) PchipInterpolator
        pchip = PchipInterpolator(df.iloc[:, 0], df.iloc[:, 1])
        if catania:
            pchip2 = PchipInterpolator(df.iloc[:, 0], df.iloc[:, 2])
        # x_interp = np.linspace(1, med, 100)  # max(df.iloc[:, 0]
        # x_interp = np.linspace(0, 24, 100)  # max(df.iloc[:, 0]
        x_interp = np.linspace(1, max(df.iloc[:, 0]), 300)  # 
        y_pchip = pchip(x_interp)
    if not do_fit_extend:
        # return x_interp, y_pchip
        if catania:
            return max(df.iloc[:, 0]), pchip, pchip2
        return max(df.iloc[:, 0]), pchip
    x_extended, y_extended = fit(x_interp, y_pchip, med, getParams=False)
    return np.concatenate((x_interp, x_extended)), np.concatenate((y_pchip, y_extended))


def preprocess_data(data_file_path, get_source_data=False, compine_syst_stat=True, header=9):
    '''preprocess the light flavour data, use weighted average to combine the 3040 and 4050 to get the 3050 result for comparation with HF'''
    v2_index = 1
    data_columns = ['PT', 'v2', 'Stat +', 'Stat -', 'Syst +', 'Syst -']
    if header==12:
        data_columns = ['PT', 'v2', 'stat +', 'stat -', 'sys +', 'sys -']
        v2_index = 3
    if not isinstance(data_file_path, list):
        data_file_path = list(data_file_path)
    if get_source_data:
        df = read_txt(data_file_path[0], header=header)
        # print(df)
        return df
    df1 = read_txt(data_file_path[0], header=header)
    df2 = read_txt(data_file_path[1], header=header)
    # print(df1.columns.tolist())
    # print(df2.columns.tolist())
    weight1 = 1 / (df1[data_columns[2]] ** 2)
    weight2 = 1 / (df2[data_columns[2]] ** 2)
    weighted_avg = (df1.iloc[:, v2_index] * weight1 + df2.iloc[:, v2_index] * weight2) / (weight1+ weight2)
    df_new = df1.copy()
    columns = df_new.columns.tolist()
    columns[v2_index] = 'v2'
    columns[0] = 'PT [GeV/c]'
    df_new.columns = columns
    df_new['v2'] = weighted_avg.values
    if compine_syst_stat:
        df_new['Total Error'] = np.sqrt(df_new[data_columns[2]]**2 + df_new[data_columns[4]]**2)
    return df_new


def preprocess_ncq(data, do_ket_nq=False, do_pt_nq=True, ismodel=False):
    '''data: dict of {name: df, ...}
        HF is not supprorted here because just one error parameter in df (total)
        return same structure as input
    '''
    pdg_db = ROOT.TDatabasePDG.Instance()
    # if do_ket_nq and do_pt_nq:
    #     print('\033[93mWARNING: must choose one from ket and nq, exit!.\033[0m')
    #     return
    if ismodel:
        print('processing model')
        for key in data.keys():
            for subkey in data[key].keys():
                if subkey == 'lc':
                    print('processing lc model')
                    _, _, _, mass = get_particle_info('LctopKpi')
                    nq = 3
                elif subkey == 'd0':
                    print('processing d0 model')
                    _, _, _, mass = get_particle_info('Dzero')
                    nq = 2
                for i in range(len(data[key][subkey])):
                    if do_ket_nq:
                        print(f'doing ket/nq for {key}{subkey}')
                        data[key][subkey][i].iloc[:, 0] = data[key][subkey][i].iloc[:, 0].apply(lambda x: kEt(mass, x))
                        data[key][subkey][i].iloc[:, 0] = data[key][subkey][i].iloc[:, 0].apply(lambda x: nq_scaling(x, nq))
                    else:
                        print(f'doing pt/nq for {key}{subkey}')
                        data[key][subkey][i].iloc[:, 0] = data[key][subkey][i].iloc[:, 0].apply(lambda x: nq_scaling(x, nq))
                    data[key][subkey][i].iloc[:, 1] = data[key][subkey][i].iloc[:, 1].apply(lambda x: nq_scaling(x, nq))
                    if key == "catania":
                        data[key][subkey][i].iloc[:, 2] = data[key][subkey][i].iloc[:, 2].apply(lambda x: nq_scaling(x, nq))
            # return data
    else:
        print('processing data')
        for key in data.keys():
            if key == "lambda":
                mass = pdg_db.GetParticle(3122).Mass()
                nq = 3
            elif key == "ks":
                mass = pdg_db.GetParticle(310).Mass()
                nq = 2
            if do_ket_nq:
                print(f'doing ket/nq for {key}')
                data[key]['kEt'] = data[key]["PT [GeV/c]"].apply(lambda x: kEt(mass, x))
                data[key]['kEt/nq'] = data[key]["kEt"].apply(lambda x: nq_scaling(x, nq))
                print(type(data[key]["PT [GeV/c]"]))
            else:
                print(f'doing pt/nq for {key}')
                data[key]['pt/nq'] = data[key]["PT [GeV/c]"].apply(lambda x: nq_scaling(x, nq))
            data[key]['v2/nq'] = data[key]["v2"].apply(lambda x: nq_scaling(x, nq))
            data[key]['Total Error/nq'] = data[key]["Total Error"].apply(lambda x: nq_scaling(x, nq))
    return data


def read_hists(file_path, markerstyle, markersize=1, colors=[], gname=['gV2PromptStat', 'gSystTotPrompt']):
    '''
    read the histograms/graph with stat and syst errors from root file
    return [histograms/graph with stat, histograms/graph with syst]
    '''
    # print(gname)
    if not isinstance(colors, list):
        colors = [colors]
    file = TFile.Open(file_path)
    if not file:
        print('error: failed to open file')
        return
    h_prompt_cent = file.Get(gname[0])
    SetObjectStyle(h_prompt_cent, color=colors[0], markerstyle=markerstyle, markersize=markersize, linewidth=2, fillalpha=0.2)
    h_prompt_systtot = file.Get(gname[1])
    SetObjectStyle(h_prompt_systtot, color=colors[0], linewidth=2)
    h_prompt_systtot.SetFillStyle(0)
    if not h_prompt_cent or not h_prompt_systtot:
        print('failed to get hist')
        return
    # h_fd_cent = file.Get('gV2FDStat')
    # SetObjectStyle(h_fd_cent, color=colors[1], markerstyle=markerstyle, linewidth=1, fillalpha=0.2)
    # h_fd_systtot = file.Get('gSystTotFD')
    # SetObjectStyle(h_fd_systtot, color=colors[1])
    # h_fd_systtot.SetFillStyle(0)
    # return [h_prompt_cent, h_fd_cent, h_prompt_systtot, h_fd_systtot]
    return [h_prompt_cent, h_prompt_systtot]


def get_band(low_x, high_x, low_y, high_y, color, doxlim=False):
    '''create and fill a band between two graphs'''
    graph_high = ROOT.TGraph(len(high_x), array.array('d', high_x), array.array('d', high_y))
    graph_high.SetLineColor(color)
    graph_high.SetLineWidth(0)
    graph_low = ROOT.TGraph(len(low_x), array.array('d', low_x), array.array('d', low_y))
    graph_low.SetLineColor(color)
    graph_low.SetLineWidth(0)

    # Find intersection of x-axis ranges for both lines
    x_min = max(min(low_x), min(high_x))
    x_max = min(max(low_x), max(high_x))
    x_min=0.5
    if doxlim:
        x_max=3.3
    # Create polygon to fill area between lines
    polyline_x = []
    polyline_y = []
    # Add points from first line
    for x, y in zip(low_x, low_y):
        if x >= x_min and x <= x_max:
            polyline_x.append(x)
            polyline_y.append(y)
    # Add points from second line (reversed)
    for x, y in zip(reversed(high_x), reversed(high_y)):
        if x >= x_min and x <= x_max:
            polyline_x.append(x)
            polyline_y.append(y)
    # Create polygon
    polyline = ROOT.TPolyLine(len(polyline_x), array.array('d', polyline_x), array.array('d', polyline_y))
    # polyline.SetLineColor(color)
    polyline.SetFillColor(color)  # Set fill color
    polyline.SetFillStyle(1001)  # Set fill style to solid
    polyline.SetLineWidth(0)

    # return [graph_low, graph_high]
    return polyline


def get_latex():
    lat_large = ROOT.TLatex()
    lat_large.SetNDC()
    lat_large.SetTextFont(42)
    # lat_large.SetTextColor(ROOT.kBlack)
    lat_large.SetTextSize(0.05)
    lat_mid = ROOT.TLatex()
    lat_mid.SetNDC()
    lat_mid.SetTextFont(42)
    lat_mid.SetTextSize(0.045)
    latLabel = ROOT.TLatex()
    latLabel.SetNDC()
    latLabel.SetTextFont(42)
    latLabel.SetTextSize(0.04)
    return lat_large, lat_mid, latLabel


def get_edges(df):
    '''Infer edges from centers data'''
    centers = np.array(df["PT [GeV/c]"])
    d = np.diff(centers) / 2
    edges = np.concatenate([[centers[0] - d[0]], centers[:-1] + d, [centers[-1] + d[-1]]])
    # Adjust edges to multiples of 0.125
    edges = np.round(edges / 0.125) * 0.125
    
    return edges


def color():
    red_palette = [
        "#8B0000",   # DarkRed
        "#CD5C5C",   # IndianRed (good for mattia, but need more testing)
        "#DC143C",   # Crimson
        "#FF6347",   # Tomato
        "#FF7F50",   # Coral
        "#E9967A",   # DarkSalmon
        "#FA8072",   # Salmon
        "#FFA07A"    # LightSalmon
    ]
    red_palette = [
        "#2F0A28",   # Wine - deep dark purple-red
        "#5E1914",   # Terracotta
        "#7C0A02",   # Bloodstone
        "#8B0000",   # DarkRed
        "#9A1B1B",   # Rusty Red
        "#A52A2A",   # Brown Red
        "#B22222",   # Firebrick
        "#C04000",   # Copper Red
        "#CD5C5C",   # IndianRed
        "#D2691E",   # Chocolate
        "#DC143C",   # Crimson
        "#E25822",   # Flame Orange
        "#E97451",   # Burnt Sienna
        "#FF4500",   # OrangeRed
        "#FF6F61",   # Coral Pink
        "#FFA07A"    # Light Salmon
    ]
    blue_palette = [
        "#191970",   # MidnightBlue
        "#4169E1",   # RoyalBlue
        "#4682B4",   # SteelBlue
        "#5F9EA0",   # CadetBlue
        "#6495ED",   # CornflowerBlue
        "#87CEEB",   # SkyBlue
        "#87CEFA",   # LightSkyBlue
        "#B0E0E6"    # PowderBlue
    ]


def fill_hist(data, hist='', columns=["PT [GeV/c]", "v2", "Stat +"]):
    # Fill histogram
    if not hist:
        hist = ROOT.TH1F("hist", "Histogram", 100, 0, 10)
    for i in range(len(data[columns[0]])):
        x = data[columns[0]][i]
        y = data[columns[1]][i]
        yerr = data[columns[2]][i]
        hist.SetBinContent(i + 1, y)
        hist.SetBinError(i + 1, yerr)
    hist.SetStats(ROOT.kFALSE)


def fill_graph(data, columns=["PT [GeV/c]", "v2", "Total Error"], compine_syst_stat=True):
    # print(columns)
    n_points = len(data[columns[0]])
    # graph = ROOT.TGraphErrors(n_points)
    graph = ROOT.TGraphAsymmErrors(n_points)
    
    # Fill graph
    for i in range(n_points):
        x = data[columns[0]][i]
        y = data[columns[1]][i]
        # yerr = data["Stat +"][i]
        if compine_syst_stat:
            yerr = data[columns[2]][i]
        else:
            yerr=0
        graph.SetPoint(i, x, y)  # Set point (x, y)
        dx=0
        # dx=x/10
        # if dx < 0.4:
        #     dx = 0.25
        if 'syst' in columns[2]:
            # print('box for syst')
            dx = x/10
            dx = dx*0.8
            if dx < 0.4:
                dx = 0.25
                dx = dx*0.8
        # graph.SetPointError(i, dx, yerr)  # Set error (dx, dy), dx=0 indicates no x-axis error
        graph.SetPointError(i, dx, dx, yerr, yerr)  # Set error (dx, dy), dx=0 indicates no x-axis error

    # Set graph style
    graph.SetMarkerStyle(20)  # Set marker style
    graph.SetMarkerSize(1.2)  # Set marker size
    graph.SetLineColor(ROOT.kBlue)  # Set line color
    graph.SetLineWidth(2)  # Set line width

    return graph


def graph_to_hist_with_errors(graph, hist_name, pt_bins, title="", use_syst_errors=False, graph_hist=''):
    """
    Convert TGraph (supports TGraphErrors) to TH1 histogram with errors, considering errors from original data
    
    Parameters:
        graph: Input TGraph or TGraphErrors object
        hist_name: Name for the output histogram
        pt_bins: List containing pt bin boundaries, e.g., [0, 1, 2, 3]
        title: Histogram title
        use_y_errors: Whether to use y-direction errors from TGraph as weights
    
    Returns:
        Converted TH1F histogram object (includes error information)
    """
    # Validate input
    if not isinstance(graph, (ROOT.TGraph, ROOT.TGraphErrors)):
        raise TypeError("Input must be a TGraph or TGraphErrors object")
    
    if len(pt_bins) < 2:
        raise ValueError("pt_bins must contain at least two boundary values")
    if use_syst_errors and not graph_hist:
        raise ValueError("If using systematic errors, must provide graph histogram")

    # Create histogram and enable error calculation
    nbins = len(pt_bins) - 1
    hist = ROOT.TH1F(hist_name, title, nbins, np.asarray(pt_bins, dtype='d'))
    hist.Sumw2()  # Critical: Enable sum of squares of errors calculation
    
    # Get graph data
    n_points = graph.GetN()
    x_vals = graph.GetX()
    y_vals = graph.GetY()
    # Check if graph has errors
    # has_errors = isinstance(graph, ROOT.TGraphErrors)
    has_errors = isinstance(graph, ROOT.TGraphAsymmErrors)
    # x_errs = graph.GetEX() if has_errors else None
    y_errs = graph.GetEYhigh() if has_errors else graph.GetEY()
    if use_syst_errors and not graph_hist:
        y_errs = graph.GetEYhigh() if has_errors else graph.GetEY()
    # Iterate through all points and fill histogram
    for ibin in range(n_points):
        x = x_vals[ibin]
        y = y_vals[ibin]
        
        hist.SetBinContent(ibin+1, y)
        hist.SetBinError(ibin+1, graph.GetEYhigh()[ibin])
    return hist


def get_interp_hist(hists_target, x_max, input_df=[], name='', cent=True): 
    '''
    Interpolated models hist and get the model value at experiment data bin center to caculate chi^2 between model and data
    input: hist 
    cent: use bin center
    return hist with value at experiment data bin center
    '''
    target_bins = get_edges_from_hist(hists_target)
    interpolate_bins = interpolate_pt_bins(target_bins)
    new_hist = create_hist_safely(name, name, interpolate_bins)
    # print(target_bins)
    # print(interpolate_bins)
    # print(new_hist.GetNbinsX()+1)
    if cent:
        new_hist = hists_target.Clone(name)
    for iPt in range(1, new_hist.GetNbinsX()+1):
        ptCent = new_hist.GetBinCenter(iPt)
        ptmax = new_hist.GetBinLowEdge(iPt+10) 
        if iPt%10==0:
            # print(ptCent)
            ptmax = new_hist.GetBinLowEdge(iPt+10)
        if cent:
            ptmax = ptCent
        if ptmax < x_max:
            if len(input_df)==1:
                new_hist.SetBinContent(iPt, input_df[0](ptCent))
            elif len(input_df)==2:
                new_hist.SetBinContent(iPt, np.mean([input_df[0](ptCent), input_df[1](ptCent)]))
        else:
            new_hist.SetBinContent(iPt, 1e-10)
        new_hist.SetBinError(iPt, 0)
    # new_hist =  rebin_safely(new_hist, name, target_bins, fixed_rebin=10)
    return new_hist


def get_edges_from_hist(hist):
    n_bins = hist.GetNbinsX() 
    bin_edges = [hist.GetBinLowEdge(i) for i in range(1, n_bins + 2)]
    return np.array(bin_edges, 'd') 


def create_hist_safely(name, title, bin_edges):
    """
    Safely create a variable-width TH1F histogram (with automatic input validation)
    
    Parameters:
        name (str): Histogram name (must be unique)
        title (str): Histogram title (can include ROOT formatting)
        bin_edges (list): Array of bin boundaries (must be monotonically increasing with length ≥ 2)
        
    Returns:
        ROOT.TH1F: Initialized histogram object
        
    Exceptions:
        ValueError: Thrown when input does not meet requirements
    """
    #-------------------------------------------
    # 1. Input validation
    #-------------------------------------------
    # Check if bin_edges is an iterable object
    if not hasattr(bin_edges, '__iter__'):
        raise ValueError("bin_edges must be an iterable object (e.g., list, np.array)")
    
    # Convert to numpy array and validate data type
    try:
        bin_edges_array = np.asarray(bin_edges, dtype='d')
    except:
        raise ValueError("Failed to convert bin_edges to a float array")
    
    # Check if bin boundary count is valid (at least 2 elements)
    if len(bin_edges_array) < 2:
        raise ValueError("bin_edges must contain at least 2 elements (minimum 1 bin)")
    
    # Check if bin boundaries are monotonically increasing
    if not np.all(np.diff(bin_edges_array) > 0):
        raise ValueError("bin_edges must be strictly monotonically increasing")
    
    #-------------------------------------------
    # 2. Calculate key parameters
    #-------------------------------------------
    n_bins = len(bin_edges_array) - 1  # Number of bins = number of boundaries - 1
    
    #-------------------------------------------
    # 3. Create histogram
    #-------------------------------------------
    # Note: ROOT constructor requires number of bins and pointer to boundary array
    hist = ROOT.TH1F(name, title, n_bins, bin_edges_array)
    
    # Automatically enable error calculation (to avoid crashes in subsequent operations)
    # hist.Sumw2()
    
    return hist


def rebin_safely(hist, new_name, new_bin_edges, is_density_hist=False, fixed_rebin=False):
    """
    Safely rebin a histogram, supporting both density and non-density histograms.
    
    Parameters:
        hist (TH1): Input histogram
        new_name (str): Name for the new histogram
        new_bin_edges (list): Array of new bin boundaries (e.g., [0, 1, 2] represents two bins)
        is_density_hist (bool): Whether the input histogram is a density histogram (content = events/bin width)
        
    Returns:
        TH1: Adjusted histogram
    """
    # Clone original histogram to avoid modifying raw data
    hist_clone = hist.Clone(f"{hist.GetName()}_clone")
    if not new_name:
        new_name = f"{hist.GetName()}_rebin"
    # If input is a density histogram, first convert to event count histogram
    if is_density_hist:
        for ibin in range(1, hist_clone.GetNbinsX() + 1):
            old_width = hist_clone.GetBinWidth(ibin)
            content = hist_clone.GetBinContent(ibin)
            error = hist_clone.GetBinError(ibin)
            hist_clone.SetBinContent(ibin, content * old_width)  # Convert to event count
            hist_clone.SetBinError(ibin, error * old_width)      # Synchronously adjust error
    
    # Perform rebinning operation
    n_new_bins = len(new_bin_edges) - 1
    hist_rebin = hist_clone.Rebin(n_new_bins, new_name, np.array(new_bin_edges, 'd'))
    
    if is_density_hist:
        # Convert rebinned histogram back to density histogram (events / new bin width)
        for ibin in range(1, hist_rebin.GetNbinsX() + 1):
            new_width = hist_rebin.GetBinWidth(ibin)
            content = hist_rebin.GetBinContent(ibin)
            error = hist_rebin.GetBinError(ibin)
            hist_rebin.SetBinContent(ibin, content / new_width)  # Adjust content to density
            hist_rebin.SetBinError(ibin, error / new_width)      # Adjust error
    elif fixed_rebin:
        # Convert rebinned histogram to density histogram (events / new bin width)
        for ibin in range(1, hist_rebin.GetNbinsX() + 1):
            new_width = fixed_rebin
            content = hist_rebin.GetBinContent(ibin)
            error = hist_rebin.GetBinError(ibin)
            hist_rebin.SetBinContent(ibin, content / new_width)  # Adjust content to density
            hist_rebin.SetBinError(ibin, error / new_width)      # Adjust error    
    return hist_rebin


def interpolate_pt_bins(pt_bins, points_per_interval=9):
    """
    Insert a specified number of points between each pair of adjacent values in pt_bins
    Parameters:
        pt_bins: Original pt interval array (e.g., [2,3,4,...24])
        points_per_interval: Number of points to insert between each pair of values (default 9)
    Returns:
        interpolated: Complete array after interpolation
    """
    interpolated = []
    # Iterate through all adjacent intervals (from 0th element to second-to-last element)
    for i in range(len(pt_bins) - 1):
        a = pt_bins[i]       # Current interval left endpoint
        b = pt_bins[i + 1]   # Current interval right endpoint
        # 1. Generate (points_per_interval + 2) points within [a, b] (including a and b)
        # 2. [1:-1] removes the duplicate right endpoint b (to avoid overlap with next interval's a)
        interval_points = np.linspace(a, b, points_per_interval + 2)[1:]
        # 3. First add left endpoint a (only the first interval needs to add a manually, subsequent intervals' a is previous interval's b)
        if i == 0:
            interpolated.append(a)
        elif i == len(pt_bins) - 2:
            interval_points = interval_points[:-1]
        # 4. Add 10 inserted points for current interval
        interpolated.extend(interval_points.tolist())
    # 5. Finally add the rightmost endpoint 24 (supplemented after all intervals are traversed)
    interpolated.append(pt_bins[-1])
    
    # Optional: Keep 2 decimal places (to avoid floating-point precision issues, e.g., 2.00, 2.10, etc.)
    interpolated = [round(x, 2) for x in interpolated]
    return interpolated


def merge_asymmetric_errors(graph_stat, graph_syst):
    """
    Final compatible version: copy old graph + modify errors with SetPointError
    Solves the problem that SetErrorYlow does not exist
    """
    # 1. Basic check: consistent number of data points
    n_stat = graph_stat.GetN()
    n_syst = graph_syst.GetN()
    if n_stat != n_syst:
        raise ValueError(f"Mismatched number of data points: {n_stat} in statistical graph vs {n_syst} in systematic graph")

    # 2. Copy statistical error graph (deep copy, preserving coordinates, x errors, and style)
    merged_graph = graph_stat.Clone(f"merged_{graph_stat.GetName()}")
    merged_graph.SetTitle(f"Merged Errors: {graph_stat.GetTitle()}")

    # 3. Calculate merged errors point by point + modify new graph errors with SetPointError
    for i in range(n_stat):
        # Check x coordinates match (to avoid data misalignment)
        x_stat = graph_stat.GetX()[i]
        x_syst = graph_syst.GetX()[i]
        if abs(x_stat - x_syst) > 1e-6:
            raise ValueError(f"X coordinate mismatch at point {i}: statistical graph {x_stat:.4f} vs systematic graph {x_syst:.4f}")

        # Step 1: Get x errors from statistical error graph (retained, not modified)
        exl_stat = graph_stat.GetErrorXlow(i)  # x lower error (statistical)
        exh_stat = graph_stat.GetErrorXhigh(i) # x upper error (statistical)

        # Step 2: Get statistical and systematic errors (y direction), merge errors
        eyl_stat = graph_stat.GetErrorYlow(i)   # y lower error (statistical)
        eyh_stat = graph_stat.GetErrorYhigh(i)  # y upper error (statistical)
        eyl_syst = graph_syst.GetErrorYlow(i)   # y lower error (systematic)
        eyh_syst = graph_syst.GetErrorYhigh(i)  # y upper error (systematic)
        
        total_eyl = np.sqrt(eyl_stat**2 + eyl_syst** 2)  # Combined y lower error
        total_eyh = np.sqrt(eyh_stat**2 + eyh_syst** 2)  # Combined y upper error

        # Step 3: Key! Use SetPointError to set errors for new graph (x errors retained, y errors replaced)
        merged_graph.SetPointError(
            i,              # ith point
            exl_stat,       # x lower error (uses statistical error)
            exh_stat,       # x upper error (uses statistical error)
            total_eyl,      # y lower error (after merging)
            total_eyh       # y upper error (after merging)
        )

    return merged_graph


def model_chi2(data_asymm, h_model, ndf=0):
    ''' caculate chi^2 between model and data'''
    print(h_model.GetName())
    chi2 = 0.0
    if  not ndf:
        for ibin in range(1,h_model.GetNbinsX()+1):
            if h_model.GetBinContent(ibin) <= 1e-9:
                h_model.SetBinContent(ibin, 0)
                ndf = ibin-1
                print(h_model.GetBinCenter(ibin-1))
                break
    residuals = []
    if isinstance(data_asymm, ROOT.TH1F):
        for ibin in range(1, ndf+1):
            residual = data_asymm.GetBinContent(ibin) - h_model.GetBinContent(ibin)  # Residual
            residuals.append(residual)
            chi2 += (residual ** 2) / (data_asymm.GetBinError(ibin) ** 2)
    else:
        x_vals = data_asymm.GetX()
        x_vals = list(x_vals)
        y_vals = data_asymm.GetY()
        y_vals = list(y_vals)
        print(x_vals, y_vals)
        model_y_vals = []
        y_errs = []
        chi2_list = []
        for i in range(ndf):
            # Extract data points (including asymmetric errors)
            # data_asymm.GetPoint(i, x_data, y_data)
            x_data = x_vals[i]
            y_data = y_vals[i]
            y_err_low = data_asymm.GetErrorYlow(i)   # Lower error (data_i - Δy_- )
            y_err_high = data_asymm.GetErrorYhigh(i) # Upper error (data_i + Δy_+ )
            # if y_err_low != y_err_high:
            #     print(f"Warning: Asymmetric errors at point {i} {x_data}, y_err_low={y_err_low}, y_err_high={y_err_high}")
            # Extract model points
            x_model = h_model.GetBinCenter(i+1)
            y_model = h_model.GetBinContent(i+1)
            model_y_vals.append(y_model)
            tolerance = 1e-6
            # Check pt matching
            if abs(x_data - x_model) > tolerance:
                raise ValueError(f"pt mismatch at point {i}: data pt={x_data}, model pt={x_model}")
            # Calculate residual
            residual = y_data - y_model
            residuals.append(residual)
            # Select corresponding error based on residual sign (avoid division by zero)
            if residual > 0:
                # Data higher than model, use upper error
                if y_err_high < 1e-10:
                    print(f"warning: {i} y_err_high=0")
                    ndf -= 1
                    continue
                sigma = y_err_high
            else:
                # Data lower than model (residual ≤ 0)
                if y_err_low < 1e-10:
                    print(f"warning: {i}: y_err_low=0")
                    ndf -= 1
                    continue
                sigma = y_err_low
            y_errs.append(sigma)
            # Accumulate chi-squared
            chi2 += (residual ** 2) / (sigma **2)
            chi2_list.append((residual ** 2) / (sigma **2))
    # print(x_data)
    print(f'chi2:{chi2:.2f}; ndf: {ndf}; chi2/ndf: {chi2/ndf:.2f}')
    # print('x_vals', list(x_vals)[:ndf])
    print('y_vals', list(y_vals)[:ndf])
    print('sigma', list(y_errs)[:ndf])
    print('model_y_vals', model_y_vals)
    print('residuals', residuals)
    print('chi2_list', chi2_list)
    if len(residuals) != ndf:
        print(residuals)
        print('len(residuals):', len(residuals))
    return chi2, ndf, chi2/ndf


def compute_ratio_graph(graph_num, graph_den):
    """
    Calculate the ratio of two TGraphAsymmErrors (graph_num / graph_den) and return a TGraphAsymmErrors of the ratio.
    
    Parameters:
        graph_num: Numerator TGraphAsymmErrors (contains y-direction errors)
        graph_den: Denominator TGraphAsymmErrors (assumed to have no errors, only its y-values are used)
    
    Returns:
        ratio_graph: TGraphAsymmErrors of the ratio, including propagated errors
    """
    # Check input types and convert TH1F to TGraphAsymmErrors if necessary
    if isinstance(graph_num, ROOT.TH1F):
        graph_num = ROOT.TGraphAsymmErrors(graph_num)
    if isinstance(graph_den, ROOT.TH1F):
        graph_den = ROOT.TGraphAsymmErrors(graph_den)
    if not (isinstance(graph_num, ROOT.TGraphAsymmErrors) and 
            isinstance(graph_den, ROOT.TGraphAsymmErrors)):
        raise TypeError("Inputs must be of type TGraphAsymmErrors")
    
    n_points = graph_num.GetN()
    if n_points != graph_den.GetN():
        raise ValueError(f"Mismatch in number of points: {n_points} vs {graph_den.GetN()}")
    
    # Extract x-axis data (access via pointer instead of passing array)
    x = []
    x_err_low = []  # x left error
    x_err_high = [] # x right error
    for i in range(n_points):
        x.append(graph_num.GetX()[i])  # Directly get i-th x value
        x_err_low.append(graph_num.GetErrorXlow(i))
        x_err_high.append(graph_num.GetErrorXhigh(i))
    
    # Extract y-values and errors (numerator and denominator)
    y_num = [graph_num.GetY()[i] for i in range(n_points)]  # Numerator y-values
    y_den = [graph_den.GetY()[i] for i in range(n_points)]  # Denominator y-values
    
    # Y-direction errors of the numerator (lower and upper errors)
    y_num_err_low = [graph_num.GetErrorYlow(i) for i in range(n_points)]
    y_num_err_high = [graph_num.GetErrorYhigh(i) for i in range(n_points)]
    
    # Calculate ratio and error propagation
    ratio_y = []
    ratio_err_low = []
    ratio_err_high = []
    for i in range(n_points):
        if y_den[i] == 0:
            raise ZeroDivisionError(f"Denominator value is zero at point {i}, cannot compute ratio")
        
        # Ratio = numerator y-value / denominator y-value
        ratio = y_num[i] / y_den[i]
        ratio_y.append(ratio)
        
        # Error propagation (denominator has no errors, numerator errors scaled proportionally)
        ratio_err_low.append(y_num_err_low[i] / y_den[i])  # Lower error
        ratio_err_high.append(y_num_err_high[i] / y_den[i])  # Upper error
    
    # Convert to numpy arrays (required for TGraphAsymmErrors constructor)
    x_np = np.array(x, dtype=float)
    ratio_y_np = np.array(ratio_y, dtype=float)
    x_err_low_np = np.array(x_err_low, dtype=float)
    x_err_high_np = np.array(x_err_high, dtype=float)
    ratio_err_low_np = np.array(ratio_err_low, dtype=float)
    ratio_err_high_np = np.array(ratio_err_high, dtype=float)
    
    # Create TGraphAsymmErrors for the ratio
    ratio_graph = ROOT.TGraphAsymmErrors(
        n_points,
        x_np, ratio_y_np,
        x_err_low_np, x_err_high_np,  # x-direction errors
        ratio_err_low_np, ratio_err_high_np  # y-direction errors
    )
    
    return ratio_graph



def scale_x_errors(graph, scale_factor=0.8, target_graph='', target_bins=[]):
    """
    Scale the x errors (both low and high errors) of a TGraphAsymmErrors object by a specified factor.
    
    Parameters:
        graph: TGraphAsymmErrors object whose x errors will be scaled
        scale_factor: Scaling factor (default is 0.8)
        target_graph: Optional TGraphAsymmErrors object. If provided, its x errors will be used 
                      as the base for scaling instead of the original errors of 'graph'
    """
    n_points = graph.GetN()  # Get the number of data points
    
    # Process each data point in a loop
    for i in range(n_points):
        # Get original x and y values (these remain unchanged)
        x = graph.GetX()[i]
        y = graph.GetY()[i]
        
        # Get original x errors (low and high)
        x_err_low = graph.GetErrorXlow(i)  # x low error (x - x_err_low)
        x_err_high = graph.GetErrorXhigh(i)  # x high error (x + x_err_high)
        
        # If target_graph is provided, use its x errors as the base
        if target_graph:
            x_err_low = target_graph.GetErrorXlow(i)  # Get x low error from target graph
            x_err_high = target_graph.GetErrorXhigh(i)  # Get x high error from target graph
            # print('Scaling x errors according to target graph', x_err_high)
        elif target_bins:
            x_err_low = (target_bins[i+1] - target_bins[i])/2
            x_err_high = x_err_low
            # print(x_err_low, x_err_high)

        # Get original y errors (these remain unchanged)
        y_err_low = graph.GetErrorYlow(i)
        y_err_high = graph.GetErrorYhigh(i)

        # Scale the x errors by the specified factor
        new_x_err_low = x_err_low * scale_factor
        new_x_err_high = x_err_high * scale_factor
        
        # Update the data point (keep x and y values unchanged, only modify x errors)
        graph.SetPoint(i, x, y)
        graph.SetPointError(i, new_x_err_low, new_x_err_high, y_err_low, y_err_high)
    
    return graph


def pdf2eps_imagemagick(pdf_paths, target_format='png'):
    """use ImageMagick convert PDF to PNG / EPS"""
    for pdf in pdf_paths:
        target_path = os.path.splitext(pdf)[0] + f'.{target_format}'
        try:
            subprocess.run([
                'convert', '-trim', 
                '-density', '300', 
                pdf, 
                target_path
            ], check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            # raise Exception(f"{pdf} ImageMagick convert failed")
            print(f"{pdf} ImageMagick convert failed")

