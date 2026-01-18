'''
Script with miscellanea utils methods for the analysis
'''

import ctypes
import sys
import numpy as np
import pandas as pd
from ROOT import TFile, TCanvas, TH1F, TF1, TList, TGraph, TGraphErrors, TGraphAsymmErrors # pylint: disable=import-error,no-name-in-module

def ComputeRatioDiffBins(hNum, hDen, uncOpt=''):
    '''
    Method to compute ratio between histograms with different bins (but compatible)

    Parameters
    ----------
    - hNum: histogram for numerator
    - hDen: histogram for denominator
    - uncOpt: uncertainty option as in ROOT.TH1.Divide

    Returns
    ----------
    - hRatio: ratio histogram
    '''

    ptMinNum = hNum.GetBinLowEdge(1)
    ptMaxNum = hNum.GetXaxis().GetBinUpEdge(hNum.GetNbinsX())
    ptMinDen = hDen.GetBinLowEdge(1)
    ptMaxDen = hDen.GetXaxis().GetBinUpEdge(hDen.GetNbinsX())
    if ptMinNum < ptMinDen:
        ptMin = ptMinDen
    else:
        ptMin = ptMinNum
    if ptMaxNum > ptMaxDen:
        ptMax = ptMaxDen
    else:
        ptMax = ptMaxNum

    if hNum.GetNbinsX() < hDen.GetNbinsX():
        if np.array(hNum.GetXaxis().GetXbins(), 'd').any(): # variable binning
            ptLimsRatio = np.array(hNum.GetXaxis().GetXbins(), 'd')
        else: # constant binning
            binWidth = hNum.GetBinWidth(1)
            ptLimsRatio = np.array([ptMinDen + iBin * binWidth for iBin in range(hNum.GetNbinsX()+1)], 'd')
    else:
        if np.array(hDen.GetXaxis().GetXbins(), 'd').any(): # variable binning
            ptLimsRatio = np.array(hDen.GetXaxis().GetXbins(), 'd')
        else: # constant binning
            binWidth = hDen.GetBinWidth(1)
            ptLimsRatio = np.array([ptMinDen + iBin * binWidth for iBin in range(hDen.GetNbinsX()+1)], 'd')
    ptLimsRatio = ptLimsRatio[(ptLimsRatio >= ptMin) & (ptLimsRatio <= ptMax)]
    nPtBins = len(ptLimsRatio)-1

    hRatio = TH1F('hRatio', f';{hNum.GetXaxis().GetTitle()};ratio', nPtBins, ptLimsRatio)
    hNumReb = TH1F('hNumReb', '', nPtBins, ptLimsRatio)
    hDenReb = TH1F('hDenReb', '', nPtBins, ptLimsRatio)

    for iPtRatio in range(1, hRatio.GetNbinsX()+1):
        deltaPt = ptLimsRatio[iPtRatio]-ptLimsRatio[iPtRatio-1]
        num, numUnc, den, denUnc = (0 for _ in range(4))
        for iPtNum in range(1, hNum.GetNbinsX()+1):
            if hNum.GetBinLowEdge(iPtNum) >= ptLimsRatio[iPtRatio-1] and \
                hNum.GetXaxis().GetBinUpEdge(iPtNum) <= ptLimsRatio[iPtRatio]:
                num += hNum.GetBinContent(iPtNum) * hNum.GetBinWidth(iPtNum)
                numUnc += hNum.GetBinError(iPtNum)**2 * hNum.GetBinWidth(iPtNum)**2 # considered uncorrelated
        hNumReb.SetBinContent(iPtRatio, num/deltaPt)
        hNumReb.SetBinError(iPtRatio, np.sqrt(numUnc)/deltaPt)
        for iPtDen in range(1, hDen.GetNbinsX()+1):
            if hDen.GetBinLowEdge(iPtDen) >= ptLimsRatio[iPtRatio-1] and \
                hDen.GetXaxis().GetBinUpEdge(iPtDen) <= ptLimsRatio[iPtRatio]:
                den += hDen.GetBinContent(iPtDen) * hDen.GetBinWidth(iPtDen)
                denUnc += hDen.GetBinError(iPtDen)**2 * hDen.GetBinWidth(iPtDen)**2 # considered uncorrelated
        hDenReb.SetBinContent(iPtRatio, den/deltaPt)
        hDenReb.SetBinError(iPtRatio, np.sqrt(denUnc)/deltaPt)

    hRatio.Divide(hNumReb, hDenReb, 1., 1., uncOpt)

    return hRatio


def ScaleGraph(graph, scaleFactor):
    '''
    Helper method to scale a TGraph

    Parameters
    ----------
    - graph: graph to scale
    - scaleFactor: scale factor
    '''
    for iPt in range(graph.GetN()):
        x, y = ctypes.c_double(), ctypes.c_double()
        graph.GetPoint(iPt, x, y)
        graph.SetPoint(iPt, x.value, y.value * scaleFactor)
        if isinstance(graph, TGraphAsymmErrors):
            yUncLow = graph.GetErrorYlow(iPt)
            yUncHigh = graph.GetErrorYhigh(iPt)
            graph.SetPointEYlow(iPt, yUncLow * scaleFactor)
            graph.SetPointEYhigh(iPt, yUncHigh * scaleFactor)
        elif isinstance(graph, TGraphErrors):
            yUncLow = graph.GetErrorYlow(iPt)
            yUncHigh = graph.GetErrorYhigh(iPt)
            graph.SetPointError(iPt, graph.GetErrorX(iPt), graph.GetErrorY(iPt) * scaleFactor)
        elif isinstance(graph, TGraph):
            continue


def ComputeRatioGraph(gNum, gDen, useDenUnc=True):
    '''
    Helper method to divide two TGraph (assuming same binning)

    Parameters
    ----------
    - gNum: graph to divide (numerator)
    - gDen: graph to divide (denominator)

    Returns
    ----------
    - gRatio: resulting graph
    '''
    if gNum.GetN() != gDen.GetN():
        print('ERROR: only graphs with same number of bins can be divided!')
        return None

    gRatio = TGraphAsymmErrors(1)
    for iPt in range(gNum.GetN()):
        x, num = ctypes.c_double(), ctypes.c_double()
        xd, den = ctypes.c_double(), ctypes.c_double()
        gNum.GetPoint(iPt, x, num)
        xUncLow = gNum.GetErrorXlow(iPt)
        xUncHigh = gNum.GetErrorXhigh(iPt)
        numUncLow = gNum.GetErrorYlow(iPt)
        numUncHigh = gNum.GetErrorYhigh(iPt)
        gDen.GetPoint(iPt, xd, den)
        denUncLow = gDen.GetErrorYlow(iPt)
        denUncHigh = gDen.GetErrorYhigh(iPt)

        ratio, ratioUncLow, ratioUncHigh = 0., 0., 0.
        if num.value != 0. and den.value != 0.:
            ratio = num.value/den.value
            if useDenUnc:
                ratioUncLow = np.sqrt((numUncLow/num.value)**2 + (denUncLow/den.value)**2) * ratio
                ratioUncHigh = np.sqrt((numUncHigh/num.value)**2 + (denUncHigh/den.value)**2) * ratio
            else:
                ratioUncLow = numUncLow / num.value * ratio
                ratioUncHigh = numUncHigh / num.value * ratio

        gRatio.SetPoint(iPt, x.value, ratio)
        gRatio.SetPointError(iPt, xUncLow, xUncHigh, ratioUncLow, ratioUncHigh)

    return gRatio


def DivideGraphByHisto(gNum, hDen, useHistoUnc=True):
    '''
    Helper method to divide a TGraph by a TH1 (assuming same binning)

    Parameters
    ----------
    - gNum: graph to divide (numerator)
    - hDen: histogram (denominator)

    Returns
    ----------
    - gRatio: resulting graph
    '''
    if gNum.GetN() != hDen.GetNbinsX():
        print('ERROR: only graphs and histos with same number of bins can be divided!')
        return None

    gRatio = TGraphAsymmErrors(0)
    for iPt in range(gNum.GetN()):
        x, num = ctypes.c_double(), ctypes.c_double()
        gNum.GetPoint(iPt, x, num)
        xUncLow = gNum.GetErrorXlow(iPt)
        xUncHigh = gNum.GetErrorXhigh(iPt)
        numUncLow = gNum.GetErrorYlow(iPt)
        numUncHigh = gNum.GetErrorYhigh(iPt)
        ptBinHisto = hDen.GetXaxis().FindBin(x.value)
        den = hDen.GetBinContent(ptBinHisto)
        if useHistoUnc:
            ratioUncLow = np.sqrt((numUncLow/num.value)**2 + (hDen.GetBinError(ptBinHisto)/den)**2) * num.value/den
            ratioUncHigh = np.sqrt((numUncHigh/num.value)**2 + (hDen.GetBinError(ptBinHisto)/den)**2) * num.value/den
        else:
            ratioUncLow = numUncLow/num.value * num.value/den
            ratioUncHigh = numUncHigh/num.value * num.value/den
        gRatio.SetPoint(iPt, x.value, num.value/den)
        gRatio.SetPointError(iPt, xUncLow, xUncHigh, ratioUncLow, ratioUncHigh)

    return gRatio


def ComputeWeightedAverage(values, weights, uncValues, uncWeights=None):
    '''
    Helper method to compute a weighted average

    Parameters
    ----------
    - values: values to be averaged
    - weights: weights for the average
    - uncValues: uncertainties on the values to be averaged
    - uncWeights: uncertainties on the weights (optional)

    Returns
    ----------
    - average: weighted average
    - unc: uncertainty on weighted average

    '''

    if len(values) != len(weights):
        print('ERROR: number of values and weights for weighted average different! Returning None')
        return None, None

    if len(values) != len(uncValues):
        print('ERROR: number of values and uncertainties for weighted average different! Returning None')
        return None, None

    if uncWeights:
        if len(weights) != len(uncWeights):
            print('ERROR: number of weights and uncertainties for weighted average different! Returning None')
            return None, None
    else:
        uncWeights = [0. for _ in weights]

    num, sumOfWeights, sumOfDerToVal, sumOfDerToWeight = (0. for _ in range(4))
    for val, wi in zip(values, weights):
        sumOfWeights += wi
        num += wi * val

    for val, wi, uncV, uncW in zip(values, weights, uncValues, uncWeights):
        sumOfDerToVal += wi**2 * uncV**2 / sumOfWeights**2
        sumOfDerToWeight += (val * sumOfWeights - num)**2 * uncW**2 / sumOfWeights**4

    average = num / sumOfWeights
    unc = np.sqrt(sumOfDerToVal + sumOfDerToWeight)

    return average, unc

