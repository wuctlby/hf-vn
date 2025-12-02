'''
Module with function definitions and fit utils
'''

from ROOT import TMath, TF1, kBlue, kGreen, TDatabasePDG, TH1D # pylint: disable=import-error,no-name-in-module

def RebinHisto(h_orig, reb, show_print, first_use = 0):
    '''
    Rebin histogram, from bin firstUse to lastUse
    Use all bins if firstUse=-1
    If ngroup is not an exact divider of the number of bins,
    the bin width is kept as reb*original width
    and the range of rebinned histogram is adapted
    '''
    
    n_bin_orig = h_orig.GetNbinsX()
    first_bin_orig = 1
    last_bin_orig = n_bin_orig
    n_bin_orig_used = n_bin_orig
    n_bin_final = n_bin_orig / reb
    
    if first_use >= 1: 
        first_bin_orig = first_use
        n_bin_final = (n_bin_orig-first_use+1) / reb
        n_bin_orig_used = n_bin_final * reb
        last_bin_orig = first_bin_orig + n_bin_orig_used - 1
    else:
        exc = n_bin_orig_used % reb
        if exc != 0: 
            n_bin_orig_used -= exc
            last_bin_orig = first_bin_orig + n_bin_orig_used - 1

    n_bin_final = round(n_bin_final)
    if (show_print):
        print(f"Rebin from {n_bin_orig} bins to {n_bin_final} bins -- Used bins = {n_bin_orig_used} in range {first_bin_orig}-{last_bin_orig}\n")
    low_lim = h_orig.GetXaxis().GetBinLowEdge(first_bin_orig)
    hi_lim = h_orig.GetXaxis().GetBinUpEdge(last_bin_orig)
    hRebin = TH1D(f"{h_orig.GetName()}-rebin", h_orig.GetTitle(), n_bin_final, low_lim, hi_lim)
    last_summed = first_bin_orig-1
    
    for iBin in range(1, n_bin_final+1):
        sum = 0.
        sume2 = 0.
        for _ in range(reb):
            sum += h_orig.GetBinContent(last_summed+1)
            sume2 += (h_orig.GetBinError(last_summed+1) * h_orig.GetBinError(last_summed+1))
            last_summed += 1
            
        hRebin.SetBinContent(iBin, sum)
        hRebin.SetBinError(iBin, TMath.Sqrt(sume2))
    
    return hRebin