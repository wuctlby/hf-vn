import ROOT
import array
import os
import numpy as np
import subprocess
from ROOT import TFile, TLegend, TCanvas, TGraphErrors, TH1F, TF1, TPad, gStyle, TColor, kFullSquare, kFullCrossX, kFullCross, kTeal, kFullDiamond, kOpenCircle, kOpenSquare
import sys
sys.path.append('./')
from plot_utils import GetLegend, GetCanvas2sub, SetGlobalStyle, SetObjectStyle, GetInvMassHistAndFit, GetV2HistAndFit, SystX
from plot_utils import get_edges_from_hist, preprocess_ncq, read_txt, DrawLineAt0, rebin_safely, preprocess, preprocess_data, read_hists, get_band, get_latex, merge_asymmetric_errors, fill_graph, model_chi2, graph_to_hist_with_errors, get_interp_hist
from plot_utils import scale_x_errors, compute_ratio_graph, pdf2eps_imagemagick
SetGlobalStyle(padleftmargin=0.16, padbottommargin=0.14, padtopmargin=0.08,
               opttitle=1, titleoffsety=1.6, labelsize=0.05, titlesize=0.05,
               labeloffset=0.01, titleoffset=1.2, labelfont=42, titlefont=42, palette=ROOT.kRainBow)

colors = [ROOT.TColor.GetColorTransparent(c, 0.6) for c in [ROOT.kRed+1, ROOT.kAzure+4, ROOT.kSpring+2, ROOT.kOrange-3, ROOT.kGray+1
                                                            , ROOT.kOrange+7, ROOT.kMagenta+2, ROOT.kBlack]]
#, ROOT.kViolet+4, ROOT.kCyan+2, ROOT.kMagenta+1, ,

def make_tgraph_from_file(filename, name="gData", title=";x;y"):
    """
    Read a two-column text file and create a ROOT TGraph.

    Parameters
    ----------
    filename : str
        Path to input file (two columns: x y).
    name : str
        Name of the TGraph.
    title : str
        Title of the TGraph (use ROOT format: 'title;x-axis;y-axis').

    Returns
    -------
    ROOT.TGraph
    """

    x_vals = []
    y_vals = []
    yunc_vals = []

    with open(filename, "r") as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            x, y, yunc = map(float, line.split())
            x_vals.append(x)
            y_vals.append(y)
            yunc_vals.append((yunc-y) * np.sqrt(2))


    graph = ROOT.TGraphErrors(len(x_vals))
    graph.SetName(name)
    graph.SetTitle(title)

    for i, (x, y, yunc) in enumerate(zip(x_vals, y_vals, yunc_vals)):
        graph.SetPoint(i, x, y)
        graph.SetPointError(i, 0, yunc)

    return graph

latex = ROOT.TLatex()
latex.SetTextFont(42)
latex.SetTextSize(0.05)
latexmedium = ROOT.TLatex()
latexmedium.SetTextFont(42)
latexmedium.SetTextSize(0.058)
latexdetail =ROOT.TLatex()
latexdetail.SetTextFont(42)
latexdetail.SetTextSize(0.04)
latexsmall =ROOT.TLatex()
latexsmall.SetTextFont(42)
latexsmall.SetTextSize(0.03)
latexdetail2 =ROOT.TLatex()
latexdetail2.SetTextFont(42)
latexdetail2.SetTextSize(0.045)
latexlarge = ROOT.TLatex()
latexlarge.SetTextFont(42)
latexlarge.SetTextSize(0.07)

def main():

    # Read histograms
    lcfile = TFile.Open('/home/spolitan/alice/hcv2-prl-figures/input/v2_Lc_d0.root')
    hlcv2 = lcfile.Get('h1')
    hlcv2.SetDirectory(0)
    d0file = TFile.Open('/home/spolitan/alice/hcv2-prl-figures/input/v2VsFracD0_OO_010.root')
    hd0v2 = d0file.Get('hV2VsPtPrompt')
    hd0v2.SetDirectory(0)
    d0filecms = TFile.Open('/home/spolitan/alice/hcv2-prl-figures/input/v2_OO_pt_cms_pbpb_6080.root')
    gd0v2cmspPb = d0filecms.Get('g0')
    dpfile = TFile.Open('/home/spolitan/alice/hcv2-prl-figures/input/v2VsFrac_010_Dplus_combined.root')
    hdplusv2 = dpfile.Get('hV2VsPtPrompt')
    hdplusv2.SetDirectory(0)
    dsfile = TFile.Open('/home/spolitan/alice/ds_v2_oo/output/181225/v2/v2VsFrac.root')
    hdsv2 = dsfile.Get('hV2VsPtPrompt')
    hdsv2.SetDirectory(0)
    k0sfile = TFile.Open('/home/spolitan/alice/hcv2-prl-figures/input/v2ofK0s_020.root')
    hk0sv2 = k0sfile.Get('v2datapoints')
    hk0sv2.SetDirectory(0)
    lambdafile = TFile.Open('/home/spolitan/alice/hcv2-prl-figures/input/v2ofLambda_020.root')
    hlambdav2 = lambdafile.Get('v2datapoints')
    hlambdav2.SetDirectory(0)
    gd0cmspp = TFile.Open('/home/spolitan/alice/hcv2-prl-figures/input/HEPData-ins1817310-v1-Table_1.root')
    d0ppv2 = gd0cmspp.Get('Table 1/Graph1D_y1')

    d0ppcent = TGraphErrors()
    d0ppcent.SetPoint(0, 3.5, 0.07201955978346933)
    d0ppcent.SetPointError(0, 1.5, 0.008751435712419337)
    SetObjectStyle(d0ppcent, markercolor=colors[-1], linecolor=colors[-1], markersize=1.2, markerstyle=20)
    
    predfile = TFile.Open('/home/spolitan/alice/hcv2-prl-figures/input/oo_v2_prediction.root')
    gpred = predfile.Get('OO_v2_prediction')
    gpred.SetLineColor(ROOT.kGray+1)
    gpred.SetLineWidth(3)
    gpred.SetFillColorAlpha(ROOT.kRed+1, 0.3)
    gpred.SetLineStyle(9)

    v2OO_prev = make_tgraph_from_file('/home/spolitan/alice/hf-vn/figures/models/input_v2d0_oo_run2.txt', name="glc_run2", title=";p_{T} (GeV/c);v_{2}")
    

    SetObjectStyle(hd0v2, markercolor=colors[1], linecolor=colors[1], markersize=1.2, markerstyle=20, linewidth=3)
    SetObjectStyle(hlcv2, markercolor=colors[0], linecolor=colors[0], markersize=1.2, markerstyle=20, linewidth=3)
    SetObjectStyle(hdplusv2, markercolor=colors[2], linecolor=colors[2], markersize=1.2, markerstyle=kFullDiamond, linewidth=3)
    SetObjectStyle(hdsv2, markercolor=colors[3], linecolor=colors[3], markersize=1.2, markerstyle=kFullCrossX, linewidth=3)
    SetObjectStyle(hk0sv2, markercolor=kTeal-5, linecolor=kTeal-5, markersize=1.2, markerstyle=kFullSquare, linewidth=3)
    SetObjectStyle(hlambdav2, markercolor=colors[5], linecolor=colors[5], markersize=1.2, markerstyle=kFullCross, linewidth=3)
    SetObjectStyle(gd0v2cmspPb, markercolor=colors[4], linecolor=colors[4], markersize=1.2, markerstyle=kOpenCircle, linewidth=3)
    SetObjectStyle(d0ppv2, markercolor=ROOT.kBlack, linecolor=ROOT.kBlack, markersize=1.2, markerstyle=kOpenSquare, linewidth=3)
    SetObjectStyle(v2OO_prev, markercolor=colors[0], linecolor=colors[0], markersize=1.2, markerstyle=20, linewidth=3)
    
    xranges = [1, 8.0]
    yranges = [-0.12, 0.44]
    axisnames = ['#it{p}_{T} (GeV/#it{c})', '#it{v}_{2}{SP, |#Delta#it{#eta}|>1.3}']

    cD0v2run3, _ = GetCanvas2sub('cDv2run3', 
                                        xranges[0],
                                        xranges[1],
                                        yranges[0],
                                        yranges[1],
                                        yranges[0],
                                        yranges[1],
                                        axisnames[0],
                                        axisnames[1],
                                     )
    cD0v2run3.cd(1)
    
    hdplusv2.GetYaxis().SetDecimals()
    hdplusv2.GetYaxis().SetNdivisions(505)
    hdplusv2.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
    hdplusv2.GetYaxis().SetTitle('#it{v}_{2}')
    hdplusv2.GetYaxis().SetRangeUser(yranges[0], yranges[1])
    hdplusv2.GetXaxis().SetRangeUser(xranges[0], xranges[1])
    hdplusv2.SetTitle('')

    hd0v2.GetYaxis().SetDecimals()
    hd0v2.GetYaxis().SetNdivisions(505)
    hd0v2.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
    hd0v2.GetYaxis().SetTitle('#it{v}_{2}{SP, |#Delta#it{#eta}|>1.3}')
    hd0v2.GetYaxis().SetRangeUser(yranges[0], yranges[1])
    hd0v2.GetXaxis().SetRangeUser(1, 8)
    hd0v2.SetTitle('')

    legD = TLegend(0.2, 0.6, 0.95, 0.8)
    legD.SetTextFont(42)
    legD.SetTextSize(0.03)
    legD.SetBorderSize(0)
    legD.SetNColumns(2) 
    legD.AddEntry(hd0v2, 'Prompt D^{0} {SP, |#Delta#it{#eta}|>1.3}', 'p')
    legD.AddEntry(gd0v2cmspPb, ' D^{0} pPb (CMS) {2PC}', 'p')
    legD.AddEntry(hdplusv2, 'Prompt D^{+} {SP, |#Delta#it{#eta}|>1.3}', 'p')
    legD.AddEntry(d0ppv2, ' D^{0} pp (CMS) {2PC}', 'p')
    legD.AddEntry(hdsv2, 'Prompt D_{s}^{+} {SP, |#Delta#it{#eta}|>1.3}', 'p')
    #legD.AddEntry(hlcv2, 'Inc. #Lambda_{c}^{+} (0#minus10%)', 'p')

    line = DrawLineAt0(xranges[0]+0.5, xranges[1])
    #hlcv2.Draw('')
    hdplusv2.Draw('')
    #gd0v2cmspPb.Draw('same PZ')
    #d0ppv2.Draw('same PZ')
    #hdsv2.Draw('same')
    hd0v2.Draw('same')
    v2OO_prev.Draw('same PZ')
    #d0ppcent.Draw('same PZ')
    latex.DrawLatexNDC(0.22, 0.92, 'ALICE, ongoing study')
    latexdetail.DrawLatexNDC(0.22, 0.86, 'O#minusO (0#minus10%) #sqrt{#it{s}_{NN}} = 5.36 TeV')
    latexdetail.DrawLatexNDC(0.8, 0.92, '|#it{y}| < 0.8')
    legD.Draw()
    line.Draw()
   
    cD0v2run3.cd(2)
    legD2 = TLegend(0.1, 0.85, 0.3, 0.95)
    legD2.SetTextFont(42)
    legD2.SetTextSize(0.04)
    legD2.SetBorderSize(0) 
    legD2.AddEntry(gpred, 'D^{0} (O#minusO, 0#minus10%), PRC 102, 041901', 'l')
    hd0v2.Draw('')
    gpred.GetXaxis().SetRangeUser(1, 8)
    gpred.Draw('same L')
    line.Draw()
    legD2.Draw()

    # Save plot
    cD0v2run3.SaveAs("./plot_oo_v2.pdf")


    # hf vs lf
    canv = TCanvas("canv_hf_lf_v2", "canv_hf_lf_v2", 800, 800)
    leghflf = TLegend(0.2, 0.6, 0.45, 0.75)
    leghflf.SetTextFont(42)
    leghflf.SetTextSize(0.035)
    leghflf.SetBorderSize(0) 
    leghflf.AddEntry(hk0sv2, 'K_{S}^{0} (0#minus20%) {2PC, |#Delta#it{#eta}|>1.2}', 'p')
    leghflf.AddEntry(hlambdav2, '#Lambda (0#minus20%) {2PC, |#Delta#it{#eta}|>1.2}', 'p')
    leghflf.AddEntry(hd0v2, 'Prompt D^{0} (0#minus10%) {SP, |#Delta#it{#eta}|>1.3}', 'p')
    leghflf.AddEntry(hdplusv2, 'Prompt D^{+} (0#minus10%) {SP, |#Delta#it{#eta}|>1.3}', 'p')
    hk0sv2.GetYaxis().SetDecimals()
    #hk0sv2.GetYaxis().SetNdivisions(505)
    hk0sv2.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
    hk0sv2.GetYaxis().SetTitle('#it{v}_{2}')
    hk0sv2.GetYaxis().CenterTitle(False)
    hk0sv2.GetYaxis().SetRangeUser(-0.1, 0.5)
    hk0sv2.GetXaxis().SetRangeUser(0.0, 8)
    hk0sv2.SetTitle('')

    line2 = DrawLineAt0(0., 8-0.5)
    hk0sv2.Draw('')
    hlambdav2.Draw('same')
    hdplusv2.Draw('same')
    hd0v2.Draw('same')
    latex.DrawLatexNDC(0.22, 0.86, 'ALICE, ongoing study')
    latexdetail.DrawLatexNDC(0.22, 0.80, 'O#minusO #sqrt{#it{s}_{NN}} = 5.36 TeV')
    #latexdetail.DrawLatexNDC(0.8, 0.86, '|#it{y}| < 0.8')
    leghflf.Draw()
    line2.Draw()
    canv.SaveAs("./plot_oo_v2_hf_lf.pdf")


    # hf pp vs OO
    canv = TCanvas("canv_hf_oo_pp_v2", "canv_hf_oo_pp_v2", 800, 800)
    leghflf = TLegend(0.2, 0.65, 0.45, 0.83)
    leghflf.SetTextFont(42)
    leghflf.SetTextSize(0.035)
    leghflf.SetBorderSize(0) 
    leghflf.AddEntry(gd0v2cmspPb, ' D^{0} pPb (CMS)', 'p')
    leghflf.AddEntry(hd0v2, 'Inc. D^{0} OO, 0#minus10%', 'p')
    leghflf.AddEntry(d0ppcent, 'Inc. D^{0} pp, 0#minus20%', 'p')
    leghflf.AddEntry(d0ppv2, ' D^{0} pp (CMS)', 'p')

    hd0v2.GetYaxis().SetDecimals()
    hd0v2.GetYaxis().SetNdivisions(505)
    hd0v2.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
    hd0v2.GetYaxis().SetTitle('#it{v}_{2}')
    hd0v2.GetYaxis().SetRangeUser(-0.1, 0.2)
    hd0v2.GetXaxis().SetRangeUser(0.1, 8)
    hd0v2.SetTitle('')

    line2 = DrawLineAt0(0.1+0.2, 8-0.5)
    hd0v2.Draw('')
    gd0v2cmspPb.Draw('same PZ')
    #d0ppcent.Draw('same PZ')
   # d0ppv2.Draw('same PZ')
    latex.DrawLatexNDC(0.22, 0.86, 'ALICE, ongoing study')
    leghflf.Draw()
    line2.Draw()
    canv.SaveAs("./plot_oo_v2_hf_oo_pp.pdf")

    # preview vs reality
    canv = TCanvas("canv_oo_v2_preview_vs_real", "canv_oo_v2_preview_vs_real", 800, 800)
    leg = TLegend(0.2, 0.68, 0.5, 0.78)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.SetBorderSize(0) 
    leg.AddEntry(hd0v2, 'Prompt D^{0} (0#minus10%, #it{L}_{int} = 5.0 nb^{-1})', 'pe')
    leg.AddEntry(v2OO_prev, 'Prompt D^{0} (0#minus80%, #it{L}_{int} = 1.0 nb^{-1})', 'pez')

    legmodel = TLegend(0.6, 0.15, 0.9, 0.25)
    legmodel.SetTextFont(42)
    legmodel.SetTextSize(0.035)
    legmodel.SetBorderSize(0)
    legmodel.AddEntry(gpred, 'PRC 102, 041901', 'l')

    hd0v2.GetYaxis().SetDecimals()
    hd0v2.GetYaxis().SetNdivisions(505)
    hd0v2.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
    hd0v2.GetYaxis().SetTitle('#it{v}_{2}')
    hd0v2.GetYaxis().SetRangeUser(-0.05, 0.22)
    hd0v2.GetXaxis().SetRangeUser(1, 8)
    hd0v2.SetTitle('') 

    line2 = DrawLineAt0(1+0.2, 8-0.5)
    hd0v2.Draw('')
    v2OO_prev.Draw('same PZ')

    gpred.GetXaxis().SetRangeUser(1, 9)
    gpred.Draw('same L')

    latex.DrawLatexNDC(0.22, 0.86, 'ALICE, performance study')
    latexdetail.DrawLatexNDC(0.22, 0.80, 'O#minusO #sqrt{#it{s}_{NN}} = 5.36 TeV')
    leg.Draw()
    legmodel.Draw()
    line2.Draw()
    canv.SaveAs("./plot_oo_v2_preview_vs_real.pdf")


main()