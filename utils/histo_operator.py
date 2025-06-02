from ROOT import TH1D, TMath

def rebin_histo(hOriginal, rebin: int, firstBin=0):
    """
    Rebin histogram by a given rebin factor, from bin firstBin to lastBin.
    Use all bins if firstBin=0.
    If ngroup is not an exact divider of the number of bins,
    the bin width is kept as rebin * original width
    and the range of rebinned histogram is adapted.
    Parameters:
    -----------
    hOriginal: TH1D
        The original histogram to be rebinned.
    rebin: int
        The number of bins to rebin the histogram into.
    firstBin: int, optional
        The first bin to use for rebinned histogram. Default is 0 (use all bins).
    Returns:
    -----------
    TH1D
        The rebinned histogram.
    """
    nBinsOrig = hOriginal.GetNbinsX()
    firstBinUsed = 1
    lastBinUsed = nBinsOrig
    nBinsOrigUsed = nBinsOrig
    nBinsFinal = nBinsOrig / rebin

    # get the first bin to use    
    if firstBin >= 1:
        firstBinUsed = firstBin
        nBinsFinal = (nBinsOrig - firstBin + 1) / rebin
        nBinsOrigUsed = nBinsFinal * rebin
        lastBinUsed = firstBinUsed + nBinsOrigUsed - 1
    # get the last bin to use
    else:
        remainder = nBinsOrigUsed % rebin
        if remainder != 0:
            nBinsOrigUsed -= remainder
            lastBinUsed = firstBinUsed + nBinsOrigUsed - 1

    nBinsFinal = round(nBinsFinal)
    print(f"Rebin {nBinsOrig} bins to {nBinsFinal} bins:")
    print(f"\tit will use {nBinsOrigUsed} bins in range {firstBinUsed}-{lastBinUsed} and will have {nBinsFinal} bins")
    
    lowLim = hOriginal.GetXaxis().GetBinLowEdge(firstBinUsed)
    hiLim = hOriginal.GetXaxis().GetBinUpEdge(lastBinUsed)
    hRebin = TH1D(f"{hOriginal.GetName()}-rebin", hOriginal.GetTitle(), nBinsFinal, lowLim, hiLim)
    lastSummed = firstBinUsed - 1

    for iBin in range(1, nBinsFinal + 1):
        sum = 0.0
        sumE2 = 0.0
        for iOrigBin in range(rebin):
            sum += hOriginal.GetBinContent(lastSummed + 1)
            sumE2 += (hOriginal.GetBinError(lastSummed + 1) * hOriginal.GetBinError(lastSummed + 1))
            lastSummed += 1
        hRebin.SetBinContent(iBin, sum)
        hRebin.SetBinError(iBin, TMath.Sqrt(sumE2))

    return hRebin