import ROOT
import sys
from ReadModel import ReadTAMUv2, ReadLGR
from plot_utils import LoadGraphAndSyst, GetCanvas, GetLegend, GetCanvas3sub, DrawLineAt0, SaveCanvas, GetPrediction, DrawStatSystEmpty, SetGlobalStyle, SetObjectStyle, GetInvMassHistAndFit, GetCanvas2sub, GetV2HistAndFit, ReadLIDOV2, HepDataHandeler,GetJPsiGraph,GetDeuGraph, compute_central_graph, GetCanvas4sub, compute_ds_d0_v2_diff

SetGlobalStyle(padleftmargin=0.16, padbottommargin=0.14, padtopmargin=0.08,
               opttitle=1, titleoffsety=1.6, labelsize=0.05, titlesize=0.05,
               labeloffset=0.01, titleoffset=1.2, labelfont=42, titlefont=42, palette=ROOT.kRainBow)

#_________________________________________________________________________________________________________________________________________
# Settings
plot_invmass_fit = False
plot_v2ffd_fit = False
plot_cutvar = False
plot_d0_prompt_np = False
plot_npd0_run2_run3 = False
plot_npd0_vsmodel = False
plot_d_v2_run3_3050 = False
plot_charm_3050 = True
plot_d0_ds_v2_run3_3050 = False
plot_d0_ds_v2_run3_3050_wdiff = False
plot_d_v2_run3_vscent = False
plot_d_v2_run3_vscent_vmodel = False
plot_d0_v2_run3_6080_vsmodel = False
plot_d_v2_run2_run3 = False
plot_ds_v2_run2_run3 = False
plot_ds_v2_vsmodel = False
plot_d_v2_vsmodel = False
plot_v2_compilation = False
plot_v2_d_vsdeuteron = False
load_v2_files = True
outdir = './'
suffix = 'pass4'

if load_v2_files:
    #_________________________________________________________________________________________________________________________________________
    # Load input files
    # Run 2
#    gist_run2_av, gist_run2_av_syst, gist_run2_av_fd  = LoadGraphAndSyst('/Users/spolitan/flowD/run2/Daverage_v2_q2_0-100_PbPb5TeV_3050.root',
#                                                       'gAverageVnPromptStat',
#                                                       'gAverageVnPromptSystData',
#                                                       'gAverageVnPromptSystFeedDown',
#                                                       ROOT.kGray+1,
#                                                       ROOT.kOpenStar)
#    gist_run2_ds, gist_run2_ds_syst, gist_run2_ds_fd  = LoadGraphAndSyst('/Users/spolitan/Desktop/CorrV2_SP_ML_pass3.root',
#                                                        'gVnPromptStat',
#                                                        'gVnPromptSystData',
#                                                        'gVnPromptSystFeedDown',
#                                                        ROOT.kGray,
#                                                        ROOT.kFullDoubleDiamond)
#    gist_run2_dp, gist_run2_dp_syst, gist_run2_dp_fd  = LoadGraphAndSyst('/Users/spolitan/flowD/input/DplusPrompt_3050_v2_EvShapeSP_q2_0-100.root',
#                                                        'gVnPromptStat',
#                                                        'gVnPromptSystTot',
#                                                        'gVnPromptSystFeedDown',
#                                                        ROOT.kSpring+3,
#                                                        ROOT.kFullDoubleDiamond)
#    gist_d0_cms_run2, gist_d0_cms_run2_syst = LoadGraphAndSyst('/Users/spolitan/alice/DmesonAnalysis/run3/flow/comparison/cms_d0_v2_3050.root',
#                                                       'graph',
#                                                       'syst',
#                                                       '',
#                                                       ROOT.kBlack,
#                                                       ROOT.kFullCrossX)
    # h
    
    gh_run2_alice_3040, ghsyst_run2_alice_3040 = HepDataHandeler('/Users/spolitan/cernbox/flowRun3/input/HEPData-ins1672822-v1-Table_5.root',
                                                 'Table 5',
                                                 '1',
                                                  stat='e1',
                                                  syst='e2')
    SetObjectStyle(gh_run2_alice_3040, linecolor=ROOT.kGray, markercolor=ROOT.kGray, markersize=1.2, markestyle=ROOT.kOpenCircle)
    SetObjectStyle(ghsyst_run2_alice_3040, linecolor=ROOT.kGray, markercolor=ROOT.kGray, markersize=1.2, markestyle=ROOT.kOpenCircle)
    gh_run2_alice_4050, ghsyst_run2_alice_4050 = HepDataHandeler('/Users/spolitan/cernbox/flowRun3/input/HEPData-ins1672822-v1-Table_6.root',
                                                 'Table 6',
                                                 '1',
                                                  stat='e1',
                                                  syst='e2')
    SetObjectStyle(gh_run2_alice_4050, linecolor=ROOT.kGray, markercolor=ROOT.kGray, markersize=1.6, markestyle=ROOT.kOpenCircle)
    SetObjectStyle(ghsyst_run2_alice_4050, linecolor=ROOT.kGray, markercolor=ROOT.kGray, markersize=1.6, markestyle=ROOT.kOpenCircle)
    gh_run2_alice_6070, ghsyst_run2_alice_6070 = HepDataHandeler('/Users/spolitan/cernbox/flowRun3/input/HEPData-ins1672822-v1-Table_8.root',
                                                 'Table 8',
                                                 '1',
                                                  stat='e1',
                                                  syst='e2')
    SetObjectStyle(gh_run2_alice_6070, linecolor=ROOT.kGray, markercolor=ROOT.kGray, markersize=1.6, markestyle=ROOT.kOpenCircle)
    SetObjectStyle(ghsyst_run2_alice_6070, linecolor=ROOT.kGray, markercolor=ROOT.kGray, markersize=1.6, markestyle=ROOT.kOpenCircle)
    # deuteron
    gist_deu_pass4_3040, gist_deu_pass4_syst_3040 = GetDeuGraph('/Users/spolitan/cernbox/flowRun3/input/deuterons/HEPData-ins1798556-v1-Table_1e.csv', ROOT.kAzure-4, ROOT.kFullSquare)
    gist_deu_pass4_4050, gist_deu_pass4_syst_4050 = GetDeuGraph('/Users/spolitan/cernbox/flowRun3/input/deuterons/HEPData-ins1798556-v1-Table_1f.csv', ROOT.kAzure-4, ROOT.kFullSquare)
    gist_deu_pass4_6070, gist_deu_pass4_syst_6070 = GetDeuGraph('/Users/spolitan/cernbox/flowRun3/input/deuterons/HEPData-ins1798556-v1-Table_1h.csv', ROOT.kAzure-4, ROOT.kFullSquare)
    gist_deu_pass4_3040_276, gist_deu_pass4_syst_3040_276 = HepDataHandeler('/Users/spolitan/cernbox/flowRun3/input/deuterons/HEPData-ins1611301-v1-Table_3.root',
                                                 'Table 3',
                                                 '5',
                                                  stat='e1',
                                                  syst='e2',
                                                 )
    gist_deu_pass4_4050_276, gist_deu_pass4_syst_4050_276 = HepDataHandeler('/Users/spolitan/cernbox/flowRun3/input/deuterons/HEPData-ins1611301-v1-Table_3.root',
                                                 'Table 3',
                                                 '6',
                                                  stat='e1',
                                                  syst='e2',
                                                 )
    SetObjectStyle(gist_deu_pass4_3040_276, linecolor=ROOT.kBlue, markercolor=ROOT.kBlue, markersize=1.2, markestyle=ROOT.kFullCrossX)
    SetObjectStyle(gist_deu_pass4_syst_3040_276, linecolor=ROOT.kBlue, markercolor=ROOT.kBlue, markersize=1.2, markestyle=ROOT.kFullCrossX)
    SetObjectStyle(gist_deu_pass4_4050_276, linecolor=ROOT.kBlue, markercolor=ROOT.kBlue, markersize=1.2, markestyle=ROOT.kFullCrossX)
    SetObjectStyle(gist_deu_pass4_syst_4050_276, linecolor=ROOT.kBlue, markercolor=ROOT.kBlue, markersize=1.2, markestyle=ROOT.kFullCrossX)
    
    gkzero = ROOT.TFile.Open('/home/spolitan/alice/hcv2-prl-figures/ste/HEPData-ins1672822-v1-Table_81_kzero.root').Get(
                                                 'Table 23/Graph1D_y1')
    glambda = ROOT.TFile.Open('/home/spolitan/alice/hcv2-prl-figures/ste/HEPData-ins1672822-v1-Table_87_lambda.root').Get('Table 24/Graph1D_y1')
    SetObjectStyle(gkzero, linecolor=ROOT.kOrange+1, markercolor=ROOT.kOrange+1, markersize=1.2, markestyle=ROOT.kFullSquare)
    SetObjectStyle(glambda, linecolor=ROOT.kAzure+4, markercolor=ROOT.kAzure+4, markersize=1.2, markestyle=ROOT.kFullDiamond)
    gkzero.SetMarkerStyle(ROOT.kFullSquare)
    gkzero.SetMarkerSize(1.6)
    glambda.SetMarkerStyle(ROOT.kFullDiamond)
    glambda.SetMarkerSize(1.6)
    gkzeroempty = gkzero.Clone('gkzeroempty')
    SetObjectStyle(gkzeroempty, markercolor=ROOT.kBlack, linecolor=ROOT.kOrange+2, markersize=1.2, markestyle=ROOT.kOpenSquare)
    glambdaempty = glambda.Clone('glambdaempty')
    SetObjectStyle(glambdaempty, markercolor=ROOT.kBlack, linecolor=ROOT.kAzure+4, markersize=1.2, markestyle=ROOT.kOpenDiamond, linewdith=1)
    gkzeroempty.SetMarkerStyle(ROOT.kOpenSquare)
    gkzeroempty.SetMarkerSize(1.6)
    glambdaempty.SetMarkerStyle(ROOT.kOpenDiamond)
    glambdaempty.SetMarkerSize(1.6)

    #_____________________________________________________________________________________________________________________________
    # PASS 4
    # D+
    gist_dp_pass4_3050, gist_dp_pass4_syst_3050, gist_dp_pass4_fd_3050 = LoadGraphAndSyst('/home/spolitan/alice/hcv2-prl-figures/ste/finalresults/DplusLatest/v2_wsyst/v2_prompt_wsyst_Dplus_3050.root',
                                                       'gvn_prompt_stat',
                                                       'tot_syst',
                                                       'fd_syst',
                                                       ROOT.kGreen+2,
                                                       ROOT.kFullDiamond)
    gist_dp_pass4_3040, gist_dp_pass4_syst_3040, gist_dp_pass4_fd_3040 = LoadGraphAndSyst('/home/spolitan/alice/hcv2-prl-figures/ste/finalresults/DplusLatest/v2_wsyst/v2_prompt_wsyst_Dplus_3040.root',
                                                       'gvn_prompt_stat',
                                                       'tot_syst',
                                                       'fd_syst',
                                                       ROOT.kGreen+2,
                                                       ROOT.kFullDiamond)
    gist_dp_pass4_4050, gist_dp_pass4_syst_4050, gist_dp_pass4_fd_4050 = LoadGraphAndSyst('/home/spolitan/alice/hcv2-prl-figures/ste/finalresults/DplusLatest/v2_wsyst/v2_prompt_wsyst_Dplus_4050.root',
                                                       'gvn_prompt_stat',
                                                       'tot_syst',
                                                       'fd_syst',
                                                       ROOT.kGreen+2,
                                                       ROOT.kFullDiamond)
    gist_dp_pass4_6080, gist_dp_pass4_syst_6080, gist_dp_pass4_fd_6080 = LoadGraphAndSyst('/home/spolitan/alice/hcv2-prl-figures/ste/finalresults/DplusLatest/v2_wsyst/v2_prompt_wsyst_Dplus_6080.root',
                                                       'gvn_prompt_stat',
                                                       'tot_syst',
                                                       'fd_syst',
                                                       ROOT.kGreen+2,
                                                       ROOT.kFullDiamond)
    # D0
    gist_d0_pass4_3050, gist_d0_pass4_syst_3050, gist_d0_pass4_fd_3050 = LoadGraphAndSyst('/home/spolitan/alice/hcv2-prl-figures/ste/finalresults/Dzero/v2_wsyst/v2_prompt_wsyst_D0_3050_finer.root',
                                                       'gvn_prompt_stat',
                                                       'tot_syst',
                                                       'fd_syst',
                                                        ROOT.kRed-4,
                                                        ROOT.kFullCircle,
                                                        markersize=1.2)
    gist_npd0_pass4_3050, gist_npd0_pass4_syst_3050, gist_npd0_pass4_fd_3050 = LoadGraphAndSyst('/home/spolitan/alice/hcv2-prl-figures/ste/finalresults/Dzero/v2_wsyst/v2_np_wsyst_D0_3050_wider.root',
                                                       'gvn_np_stat',
                                                       'tot_syst',
                                                       'fd_syst',
                                                        ROOT.kAzure,
                                                        ROOT.kFullCircle,
                                                        markersize=1.2)
    gist_d0_pass4_3050_wider, gist_d0_pass4_syst_3050_wider, gist_d0_pass4_fd_3050_wider = LoadGraphAndSyst('/home/spolitan/alice/hcv2-prl-figures/ste/finalresults/Dzero/v2_wsyst/v2_prompt_wsyst_D0_3050_wider.root',
                                                       'gvn_prompt_stat',
                                                       'tot_syst',
                                                       'fd_syst',
                                                        ROOT.kRed-4,
                                                        ROOT.kFullCircle,
                                                        markersize=1.2)
    gist_d0_pass4_3040, gist_d0_pass4_syst_3040, gist_d0_pass4_fd_3040 = LoadGraphAndSyst('/home/spolitan/alice/hcv2-prl-figures/ste/finalresults/Dzero/v2_wsyst/v2_prompt_wsyst_D0_3040.root',
                                                       'gvn_prompt_stat',
                                                       'tot_syst',
                                                       'fd_syst',
                                                        ROOT.kRed-4,
                                                        ROOT.kFullCircle,
                                                        markersize=1.2)
    gist_d0_pass4_4050, gist_d0_pass4_syst_4050, gist_d0_pass4_fd_4050 = LoadGraphAndSyst('/home/spolitan/alice/hcv2-prl-figures/ste/finalresults/Dzero/v2_wsyst/v2_prompt_wsyst_D0_4050.root',
                                                       'gvn_prompt_stat',
                                                       'tot_syst',
                                                       'fd_syst',
                                                        ROOT.kRed-4,
                                                        ROOT.kFullCircle,
                                                        markersize=1.2)
    gist_d0_pass4_6080, gist_d0_pass4_syst_6080, gist_d0_pass4_fd_6080 = LoadGraphAndSyst('/home/spolitan/alice/hcv2-prl-figures/ste/finalresults/Dzero/v2_wsyst/v2_prompt_wsyst_D0_6080.root',
                                                       'gvn_prompt_stat',
                                                       'tot_syst',
                                                       'fd_syst',
                                                        ROOT.kRed-4,
                                                        ROOT.kFullCircle,
                                                        markersize=1.2)
    # Ds
    gist_ds_pass4_3040, gist_ds_pass4_syst_3040, gist_ds_pass4_fd_3040 = LoadGraphAndSyst('/home/spolitan/alice/hcv2-prl-figures/ste/finalresults/Ds/v2_wsyst/v2_prompt_wsyst_Ds_3040.root',
                                                       'gvn_prompt_stat',
                                                       'tot_syst',
                                                       'fd_syst',
                                                       ROOT.kAzure+2,
                                                       ROOT.kFullCross)
    gist_ds_pass4_4050, gist_ds_pass4_syst_4050, gist_ds_pass4_fd_4050 = LoadGraphAndSyst('/home/spolitan/alice/hcv2-prl-figures/ste/finalresults/Ds/v2_wsyst/v2_prompt_wsyst_Ds_4050.root',
                                                       'gvn_prompt_stat',
                                                       'tot_syst',
                                                       'fd_syst',
                                                       ROOT.kAzure+2,
                                                       ROOT.kFullCross)
    gist_ds_pass4_6080, gist_ds_pass4_syst_6080, gist_ds_pass4_fd_6080 = LoadGraphAndSyst('/home/spolitan/alice/hcv2-prl-figures/ste/finalresults/Ds/v2_wsyst/v2_prompt_wsyst_Ds_6080.root',
                                                       'gvn_prompt_stat',
                                                       'tot_syst',
                                                       'fd_syst',
                                                       ROOT.kAzure+2,
                                                       ROOT.kFullCross)
    gist_ds_pass4_3050, gist_ds_pass4_syst_3050, gist_ds_pass4_fd_3050 = LoadGraphAndSyst('/home/spolitan/alice/hcv2-prl-figures/ste/finalresults/Ds/v2_wsyst/v2_prompt_wsyst_Ds_3050.root',
                                                       'gvn_prompt_stat',
                                                       'tot_syst',
                                                       'fd_syst',
                                                       ROOT.kAzure+2,
                                                       ROOT.kFullCross)
    # JPsi
    # gist_jpsi_pass4_3050, gist_jpsi_pass4_syst_3050 = GetJPsiGraph('/Users/spolitan/cernbox/flowRun3/input/jpsi_run3/v2_vs_pt_30_50.txt', ROOT.kOrange+6, ROOT.kFullSquare)
    # Lc
    gist_lc_pass4_3050, gist_lc_pass4_syst_3050, gist_lc_pass4_fd_3050 = LoadGraphAndSyst('/home/spolitan/alice/hcv2-prl-figures/input/lc-prompt-allpt-wTotsyst.root',
                                                       'gV2PromptStat',
                                                       'gSystTotPrompt',
                                                       'gSystTotPrompt',
                                                       ROOT.kGreen+2,
                                                       ROOT.kFullCrossX)

#_____________________________________________________________________________________________________________________________
# Plot detail
axisname = ';#it{p}_{T} (GeV/#it{c}); #it{v}_{2}{SP, |#Delta#it{#eta}| > 1.3}'
axisnamevareta = ';#it{p}_{T} (GeV/#it{c}); #it{v}_{2}'
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


#_____________________________________________________________________________________________________________________________
# Plotting
#_____________________________________________________________________________________
# D meson mass fit in Run 3 pass4 3040
if plot_invmass_fit:
    hInvMassD0, fMassTotD0, fMassBkgD0, hV2VsMassD0, fV2TotD0, fV2BkgD0, fMassReflD0 = GetInvMassHistAndFit('/Users/spolitan/cernbox/flowRun3/output/pass4/finalresults/D0_raw_yields_uncorr_00.root', 4, 5, 7, hasreflections=True)
    hInvMassDp, fMassTotDp, fMassBkgDp, hV2VsMassDp, fV2TotDp, fV2BkgDp = GetInvMassHistAndFit('/Users/spolitan/cernbox/flowRun3/output/pass4/finalresults/Dp_raw_yields_combined_00_for_figure_final.root', 1.5, 2.0, 1)
    hInvMassDs, fMassTotDs, fMassBkgDs, hV2VsMassDs, fV2TotDs, fV2BkgDs = GetInvMassHistAndFit('/Users/spolitan/cernbox/flowRun3/output/pass4/finalresults/Ds_raw_yields_pass4_ds_full_new_3040_01.root', 4, 5, 3)
    xmins_mass = [1.72, 1.76, 1.76]
    xmaxs_mass = [2.0, 1.95, 2.08]
    ymins_mass = [1, 9000, -490]
    ymaxs_mass = [37252, 45252, 2650]
    ymins_v2 = [0.16, 0.095, 0.08]
    ymaxs_v2 = [0.285, 0.165, 0.52]
    axisnamebottoms = [';#it{M}(K#pi) (GeV/#it{c}^{2});#it{v}_{2}^{obs.}{SP, |#Delta#it{#eta}|>1.3}',
                       ';#it{M}(K#pi#pi) (GeV/#it{c}^{2});#it{v}_{2}^{obs.}{SP, |#Delta#it{#eta}|>1.3}',
                       ';#it{M}(KK#pi) (GeV/#it{c}^{2});#it{v}_{2}^{obs.}{SP, |#Delta#it{#eta}|>1.3}',
                       ]
    axisnametops = [f';;Counts per {hInvMassD0.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}',
                    f';;Counts per {hInvMassD0.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}',
                    f';;Counts per {hInvMassDs.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}']

    legD = GetLegend(xmax=0.5, ncolumns=1, ymin=0.05, ymax=0.3, textsize=0.055, xmin=0.22)
    legD.AddEntry(hInvMassDs, 'Data', 'p')
    legD.AddEntry(fV2TotDs, 'Total fit function', 'l')
    legD.AddEntry(fV2BkgDs, 'Combinatorial background', 'l')
    legD.AddEntry(fMassReflD0, 'K#minus#pi reflected', 'f')

    cD0v2run3, _ = GetCanvas2sub('cDv2run3', 
                                     xmins_mass[0],
                                     xmaxs_mass[0],
                                     ymins_mass[0],
                                     ymaxs_mass[0],
                                     ymins_v2[0],
                                     ymaxs_v2[0],
                                     axisnametops[0],
                                     axisnamebottoms[0]
                                     )
    cD0v2run3.cd(1)
    hInvMassD0.Draw('esame')
    fMassReflD0.Draw('same')
    fMassBkgD0.Draw('same')
    fMassTotD0.Draw('same')
    legD.Draw()
    latexlarge.DrawLatexNDC(0.22, 0.86, 'ALICE Preliminary')
    latexmedium.DrawLatexNDC(0.22, 0.78, 'Pb#minusPb, 30#minus40%, #sqrt{#it{s}_{NN}} = 5.36 TeV')
    latexmedium.DrawLatexNDC(0.22, 0.7, 'D^{0} #rightarrow K^{#font[122]{-}}#pi^{+} and charge conj.')
    latexmedium.DrawLatexNDC(0.22, 0.62, '4 < #it{p}_{T} < 5 GeV/#it{c}')
    cD0v2run3.cd(2)
    latexmedium.DrawLatexNDC(0.22, 0.92, '#it{v}_{2}^{obs.} = 0.165 #pm 0.003')
    hV2VsMassD0.Draw('esame')
    fV2BkgD0.Draw('same')
    fV2TotD0.Draw('same')
    SaveCanvas(cD0v2run3, f'{outdir}D0_invmassfit', suffix)
    input()

    cDpv2run3, _  = GetCanvas2sub('cDv2run3', 
                                     xmins_mass[1],
                                     xmaxs_mass[1],
                                     ymins_mass[1],
                                     ymaxs_mass[1],
                                     ymins_v2[1],
                                     ymaxs_v2[1],
                                     axisnametops[1],
                                     axisnamebottoms[1]
                                     )
    legD = GetLegend(xmax=0.5, ncolumns=1, ymin=0.05, ymax=0.24, textsize=0.06, xmin=0.22)
    legD.AddEntry(hInvMassDs, 'Data', 'p')
    legD.AddEntry(fV2TotDs, 'Total fit function', 'l')
    legD.AddEntry(fV2BkgDs, 'Combinatorial background', 'l')
    cDpv2run3.cd(1)
    hInvMassDp.Draw('esame')
    fMassBkgDp.Draw('same')
    fMassTotDp.Draw('same')
    legD.Draw()
    latexlarge.DrawLatexNDC(0.22, 0.86, 'ALICE Preliminary')
    latexmedium.DrawLatexNDC(0.22, 0.78, 'Pb#minusPb, 30#minus40%, #sqrt{#it{s}_{NN}} = 5.36 TeV')
    latexmedium.DrawLatexNDC(0.22, 0.7, 'D^{+} #rightarrow K^{#font[122]{-}}#pi^{+}#pi^{+} and charge conj.')  
    latexmedium.DrawLatexNDC(0.22, 0.62, '1.5 < #it{p}_{T} < 2.0 GeV/#it{c}')
    cDpv2run3.cd(2)
    latexmedium.DrawLatexNDC(0.22, 0.92, '#it{v}_{2}^{obs.} = 0.085 #pm 0.004')
    hV2VsMassDp.Draw('esame')
    fV2BkgDp.Draw('same')
    fV2TotDp.Draw('same')
    SaveCanvas(cDpv2run3, f'{outdir}Dp_invmassfit', suffix)
    input()

    cDsv2run3, _  = GetCanvas2sub('cDv2run3', 
                                     xmins_mass[2],
                                     xmaxs_mass[2],
                                     ymins_mass[2],
                                     ymaxs_mass[2],
                                     ymins_v2[2],
                                     ymaxs_v2[2],
                                     axisnametops[2],
                                     axisnamebottoms[2]
                                     )
    cDsv2run3.cd(1)
    hInvMassDs.Draw('esame')
    fMassBkgDs.Draw('same')
    fMassTotDs.Draw('same')
    legD.Draw()
    latexlarge.DrawLatexNDC(0.22, 0.86, 'ALICE Preliminary')
    latexmedium.DrawLatexNDC(0.22, 0.78, 'Pb#minusPb, 30#minus40%, #sqrt{#it{s}_{NN}} = 5.36 TeV')
    latexmedium.DrawLatexNDC(0.22, 0.7, 'D_{s}^{+} #rightarrow #phi#pi^{+} #rightarrow K^{+}K^{#font[122]{-}}#pi^{+} and charge conj.')
    latexmedium.DrawLatexNDC(0.22, 0.62, '4 < #it{p}_{T} < 5 GeV/#it{c}')
    cDsv2run3.cd(2)
    latexmedium.DrawLatexNDC(0.22, 0.92, '#it{v}_{2}^{obs.} = 0.143 #pm 0.045')
    hV2VsMassDs.Draw('esame')
    fV2BkgDs.Draw('same')
    fV2TotDs.Draw('same')
    SaveCanvas(cDsv2run3, f'{outdir}Ds_invmassfit', suffix)
    input()

#_____________________________________________________________________________________
# D meson v2 vs ffd fit in Run 3 pass4 3040
if plot_v2ffd_fit:
    gv2D0, hv2D0, tf1D0 = GetV2HistAndFit('/Users/spolitan/cernbox/flowRun3/output/pass4/finalresults/D0_V2VsFrac_combined_3040.root',
                                   'pt_35_40', 4, 5, 7)
    gv2Dp, hv2Dp, tf1Dp = GetV2HistAndFit('/Users/spolitan/cernbox/flowRun3/output/pass4/finalresults/Dp_V2VsFrac_combined_for_figure.root', 'pt_25_30', 3.0, 3.5, 4)
    inFile = ROOT.TFile.Open('/Users/spolitan/cernbox/flowRun3/output/pass4/finalresults/Ds_V2VsFrac_pass4_ds_full_new_3040.root')
    gv2Ds = inFile.Get('pt_40_50/gV2VsFrac')
    hv2Ds = inFile.Get('pt_40_50/hV2VsFrac')
    tf1Ds = gv2Ds.GetFunction("linear")
    SetObjectStyle(gv2Ds, color=ROOT.kBlack, markersize=1.2, linewidth=1)
    SetObjectStyle(hv2Ds, linewidth=1, linecolor=ROOT.kAzure+4, markersize=0, alpha=0.1)
    SetObjectStyle(hv2Dp, linewidth=1, linecolor=ROOT.kAzure+4, markersize=0, alpha=0.1)
    SetObjectStyle(hv2D0, linewidth=1, linecolor=ROOT.kAzure+4, markersize=0, alpha=0.1)
    SetObjectStyle(tf1Ds, linewidth=1, linecolor=ROOT.kRed-4, markersize=0, linestyle=9)
    cDv2run3, hframes = GetCanvas3sub('cDv2run3', ';Non-prompt fraction; #it{v}_{2}^{obs.}{SP, |#Delta#it{#eta}| > 1.3}')

    legD = GetLegend(xmax=0.8, ncolumns=1, ymin=0.202, ymax=0.402, textsize=0.045, xmin=0.5)
    legD.AddEntry(tf1Ds, 'Linear fit', 'l')
    hv2Dleg = hv2D0.Clone('hv2Dleg')
    hv2Dleg.SetFillColorAlpha(ROOT.kAzure+4, 0.4)
    legD.AddEntry(hv2Dleg, '68% confidence level', 'f')

    cDv2run3.cd(1)
    hframes[0].GetYaxis().SetRangeUser(0.0, 0.28)
    hframes[0].GetXaxis().SetRangeUser(-0.02, 1.02)
    hframes[0].GetYaxis().SetLabelSize(0.06)
    hframes[0].GetXaxis().SetLabelSize(0.06)
    hv2D0.Draw('same')
    gv2D0.Draw('pez same')
    tf1D0.Draw('same')
    latexlarge.DrawLatexNDC(0.22, 0.92, 'ALICE Preliminary')
    latexdetail2.DrawLatexNDC(0.22, 0.84, 'Pb#minusPb, 30#minus40%, #sqrt{#it{s}_{NN}} = 5.36 TeV')
    latexdetail2.DrawLatexNDC(0.22, 0.76, 'D^{0} #rightarrow K^{#font[122]{-}}#pi^{+} and charge conj.')
    latexdetail2.DrawLatexNDC(0.22, 0.68, '4 < #it{p}_{T} < 5 GeV/#it{c}')
    latexdetail2.DrawLatexNDC(0.22, 0.22, f'#chi^{{2}}/ndf = {tf1D0.GetChisquare()/tf1D0.GetNDF():.2f}')
    legD.Draw()

    cDv2run3.cd(2)
    hframes[1].GetYaxis().SetRangeUser(1.e-10, 0.28)
    hframes[1].GetXaxis().SetRangeUser(-0.02, 1.02)
    hframes[1].GetYaxis().SetLabelSize(0.06)
    hframes[1].GetXaxis().SetLabelSize(0.06)
    gv2Dp.Draw('pez same')
    hv2Dp.Draw('same')
    tf1Dp.Draw('same')
    latexdetail2.DrawLatexNDC(0.08, 0.92, 'D^{+} #rightarrow K^{#font[122]{-}}#pi^{+}#pi^{+} and charge conj.')
    latexdetail2.DrawLatexNDC(0.08, 0.86, '3.0 < #it{p}_{T} < 3.5 GeV/#it{c}')
    latexdetail2.DrawLatexNDC(0.08, 0.22, f'#chi^{{2}}/ndf = {tf1Dp.GetChisquare()/tf1Dp.GetNDF():.2f}')

    cDv2run3.cd(3)
    hframes[2].GetYaxis().SetRangeUser(0.0, 0.28)
    hframes[2].GetXaxis().SetRangeUser(-0.02, 1.02)
    hframes[2].GetYaxis().SetLabelSize(0.06)
    hframes[2].GetXaxis().SetLabelSize(0.06)
    hv2Ds.Draw('same')
    gv2Ds.Draw('pez same')
    tf1Ds.Draw('same')
    latexdetail2.DrawLatexNDC(0.08, 0.92, 'D_{s}^{+} #rightarrow #phi#pi^{+} #rightarrow K^{+}K^{#font[122]{-}}#pi^{+} and charge conj.')
    latexdetail2.DrawLatexNDC(0.08, 0.86, '4 < #it{p}_{T} < 5 GeV/#it{c}')
    latexdetail2.DrawLatexNDC(0.08, 0.22, f'#chi^{{2}}/ndf = {tf1Ds.GetChisquare()/tf1Ds.GetNDF():.2f}')

    cDv2run3.Update()
    SaveCanvas(cDv2run3, f'{outdir}D_v2ffd_fit', suffix)
    input()

#_____________________________________________________________________________________
# D0 cutvariation
if plot_cutvar:
    infile = ROOT.TFile.Open('/Users/spolitan/cernbox/flowRun3/output/pass4/finalresults/CutVarFrac_corr.root')
    hry = infile.Get('pt3.5_4.0/hRawYieldsVsCutPt_pT3.5_4.0')
    hry_prompt = infile.Get('pt3.5_4.0/hRawYieldPromptVsCut_pT3.5_4.0')
    hry_fd = infile.Get('pt3.5_4.0/hRawYieldFDVsCut_pT3.5_4.0')
    hry_sum = infile.Get('pt3.5_4.0/hRawYieldsVsCutReSum_pT3.5_4.0')
    heff_prompt = infile.Get('pt3.5_4.0/hEffPromptVsCut_pT3.5_4.0')
    heff_fd = infile.Get('pt3.5_4.0/hEffFDVsCut_pT3.5_4.0')
    hfrac_prompt = infile.Get('pt3.5_4.0/hPromptFracVsCut_pT3.5_4.0')
    hfrac_fd = infile.Get('pt3.5_4.0/hFDFracVsCut_pT3.5_4.0')
    infile2 = ROOT.TFile.Open('/Users/spolitan/cernbox/flowRun3/output/pass4/finalresults/cicio.root')
    hcorr = infile2.Get('hCorrMatrixCutSets_pT3.5_4.0')

    SetObjectStyle(hry, fillalpha=0., linecolor=ROOT.kBlack, markercolor=ROOT.kBlack, markerstyle=ROOT.kFullCircle, fillstyle=1001, markersize=1.6, linewidth=4)
    SetObjectStyle(hry_prompt, fillalpha=0.3, linecolor=ROOT.kRed+1, markercolor=ROOT.kRed+1, markerstyle=ROOT.kFullCircle, fillstyle=3245, markersize=0, fillcolor=ROOT.kRed+1)
    SetObjectStyle(hry_fd, fillalpha=0.3, linecolor=ROOT.kAzure+4, markercolor=ROOT.kAzure+4, markerstyle=ROOT.kFullCircle, fillstyle=3254, markersize=0, fillcolor=ROOT.kAzure+4)
    SetObjectStyle(hry_sum, fillalpha=0., linecolor=ROOT.kGreen+2, markercolor=ROOT.kGreen+2, markerstyle=ROOT.kFullCircle, fillstyle=1001, markersize=0)
    SetObjectStyle(heff_prompt, fillalpha=0., linecolor=ROOT.kRed+1, markercolor=ROOT.kRed+1, markerstyle=ROOT.kFullCircle, fillstyle=1001, markersize=1.6)
    SetObjectStyle(heff_fd, fillalpha=0., linecolor=ROOT.kAzure+4, markercolor=ROOT.kAzure+4, markerstyle=ROOT.kFullCircle, fillstyle=1001, markersize=1.6)
    SetObjectStyle(hfrac_prompt, fillalpha=0., linecolor=ROOT.kRed+1, markercolor=ROOT.kRed+1, markerstyle=ROOT.kFullCircle, fillstyle=1001, markersize=1.6)
    SetObjectStyle(hfrac_fd, fillalpha=0., linecolor=ROOT.kAzure+4, markercolor=ROOT.kAzure+4, markerstyle=ROOT.kFullCircle, fillstyle=1001, markersize=1.6)

    leg = GetLegend(header=' ', xmax=0.85, ncolumns=1, ymin=0.8, xmin=0.65, ymax=0.9, textsize=0.045)
    leg.AddEntry(heff_prompt, 'Prompt', 'p')
    leg.AddEntry(heff_fd, 'Non-prompt', 'p')

    legrun2 = GetLegend(header='3.5 < #it{p}_{T} < 4 GeV/#it{c}', xmax=0.85, ncolumns=1, ymin=0.4, xmin=0.4, ymax=0.7, textsize=0.045)
    legrun2.AddEntry(hry, 'Data', 'p')
    legrun2.AddEntry(hry_prompt, 'Prompt D^{0}', 'f')
    legrun2.AddEntry(hry_fd, 'Non-prompt D^{0}', 'f')
    legrun2.AddEntry(hry_sum, 'Total', 'l')

    cDv2run3 = ROOT.TCanvas('', '', 800, 800)
    cDv2run3.cd(1)
    pad = ROOT.gPad
    frame = pad.DrawFrame(0.5, 0, 26.5, 4.e+5, ';BDT-based selection;Raw yields') 
    frame.SetTitle("")
    frame.GetYaxis().SetTitleSize(0.06)
    frame.GetXaxis().SetTitleSize(0.06)
    frame.GetZaxis().SetTitleSize(0.06)
    frame.GetYaxis().SetTitleOffset(1.3)
    frame.GetXaxis().SetTitleOffset(1.)
    frame.GetYaxis().SetDecimals()
    frame.GetYaxis().SetMaxDigits(2)
    #cDv2run3.Divide(2, 2)  # Remove spacing between pads

    frames = []
    for i in range(4):
        cDv2run3.cd(i + 1)
        pad = ROOT.gPad
        #pad.SetLeftMargin(0.18 if i % 2 == 0 else 0.0)  # First pad in each row has a y-axis margin
        #pad.SetRightMargin(0.05 if i % 2 == 1 else 0.0)  # Last pad in each row has right margin
        
        if i == 0: # efficiency
           frame = pad.DrawFrame(0.5, 5.e-4, 26.5, 1.0, ';BDT cut;#varepsilon#timesAcc.') 
        elif i == 1: #ry
            frame = pad.DrawFrame(0.5, 0, 26.5, 3.e+5, ';BDT cut;Raw yields') 
        elif i == 2: #corr
            pad.SetRightMargin(0.18)
            frame = pad.DrawFrame(0.5, 0.5, 26.5, 26.5, ';BDT cut;BDT cut') 
        elif i == 3: #frac
            frame = pad.DrawFrame(0.5, 0, 26.5, 1.1, ';BDT cut;Non-prompt fraction') 
            
        #frame = pad.DrawFrame(-0.5, -0.20, 25, 0.62, axisname)
        frame.SetTitle("")
        
        frame.GetYaxis().SetTitleSize(0.06)
        frame.GetXaxis().SetTitleSize(0.06)
        frame.GetZaxis().SetTitleSize(0.06)
        frame.GetYaxis().SetTitleOffset(1.3) if i !=0 else frame.GetYaxis().SetTitleOffset(1.)
        frame.GetYaxis().SetDecimals()
        frame.GetYaxis().SetMaxDigits(2)
        if i == 2:
            frame.GetZaxis().SetDecimals()
            frame.GetZaxis().SetMaxDigits(1)
            frame.GetZaxis().SetMoreLogLabels()

        frames.append(frame)
    cDv2run3.cd(3)#.SetLogz()
    hcorr.GetZaxis().SetRangeUser(5.e-1, 1.)
    hcorr.GetZaxis().SetTitle('#rho')
    hcorr.GetZaxis().SetTitleOffset(1.4)
    hcorr.GetZaxis().SetDecimals()
    hcorr.Draw('colz same')


    cDv2run3.cd(2)
    hry_prompt.Draw('hist same')
    hry_fd.Draw('hist same')
    hry_sum.Draw('hist same')
    hry.Draw('same')
    legrun2.Draw('same')
    latex.DrawLatexNDC(0.20, 0.86, 'ALICE Preliminary')
    latexdetail.DrawLatexNDC(0.20, 0.8, 'Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV, 30#minus40%')
    latexdetail.DrawLatexNDC(0.2, 0.74, 'D^{0} #rightarrow K^{#font[122]{-}}#pi^{+} and charge conj.')
    #latexdetail.DrawLatexNDC(0.2, 0.68, '3.5 < #it{p}_{T} < 4 GeV/#it{c}')

    cDv2run3.cd(1).SetLogy()
    heff_prompt.Draw('same')
    heff_fd.Draw('same')

    cDv2run3.cd(4)
    hfrac_prompt.Draw('same')
    hfrac_fd.Draw('same')
    leg.Draw('same')


    cDv2run3.Update()
    SaveCanvas(cDv2run3, f'{outdir}D0_cutvar', suffix)
    input()

#_____________________________________________________________________________________
# prompt vs np D0 meson v2
if plot_d0_prompt_np:
    leg = GetLegend(header='D^{0}', xmax=0.5, ncolumns=1, ymin=0.58, xmin=0.196)
    legempty = leg.Clone('legempty')
    legempty.SetHeader(' ')
    leg.AddEntry(gist_d0_pass4_3050, 'Prompt', 'p')
    leg.AddEntry(gist_npd0_pass4_3050, 'Non-prompt', 'p')

    gtamuv2D0 = GetPrediction('/Users/spolitan/alice/DmesonAnalysis/models/tamu/PromptD_TAMU_v2_5TeV_3050.txt',
                            1, 10, 1, True)
    gtamuv2npD0 = GetPrediction('/Users/spolitan/alice/DmesonAnalysis/models/tamu/NonPromptD0_TAMU_v2_5TeV_3050.txt',
                            1, 148, 0.1, False)
    gtamuv2D0.SetLineColor(ROOT.kRed)
    gtamuv2D0.SetLineStyle(1)
    gtamuv2D0.SetLineWidth(0)
    gtamuv2D0.SetFillColorAlpha(ROOT.kRed, 0.2)
    gtamuv2npD0.SetLineColor(ROOT.kAzure-9)
    gtamuv2npD0.SetLineStyle(1)
    gtamuv2npD0.SetLineWidth(3)
    gtamuv2npD0.SetFillColorAlpha(ROOT.kAzure-4, 0.2)
    legmodels = GetLegend(header='TAMU  #sqrt{#it{s}_{NN}} = 5.02 TeV', xmax=0.9, ncolumns=1, ymin=0.58, xmin=0.496)
    legmodels.AddEntry(gtamuv2D0, 'Prompt D^{0}', 'f')
    legmodels.AddEntry(gtamuv2npD0, 'Non-prompt D^{0}', 'l')

    cDv2run3, hframe = GetCanvas('cDv2run3', axisname, ymin=-0.12, ymax=0.45, xmax=25)
    line = DrawLineAt0(0., 25)
    gtamuv2D0.Draw('3c same')
    gtamuv2npD0.Draw('3c same')
    DrawStatSystEmpty(gist_d0_pass4_3050, gist_d0_pass4_syst_3050, False, legempty, markersize=1.2)
    DrawStatSystEmpty(gist_npd0_pass4_3050, gist_npd0_pass4_syst_3050, False, legempty, markersize=1.2)

    leg.Draw('same')
    legempty.Draw('same')
    legmodels.Draw('same')
    latex.DrawLatexNDC(0.20, 0.86, 'ALICE Preliminary')
    latexdetail.DrawLatexNDC(0.20, 0.8, 'Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV, 30#minus50%')
    latexdetail.DrawLatexNDC(0.8, 0.86, '|#it{y}| < 0.8')
    cDv2run3.Update()
    SaveCanvas(cDv2run3, f'{outdir}D0_v2_prompt_np_run3', suffix)
    input()    

#_____________________________________________________________________________________
# np D0 run2 vs run3
if plot_npd0_run2_run3:
    cRun3D0vsModels, hframe = GetCanvas('cRun3D0vsModels', ';#it{p}_{T} (GeV/#it{c}); #it{v}_{2}{SP}', ymin=-0.1, ymax=0.35)

    _ = DrawLineAt0(0., 40)
    gnpd0_run2_alice, gsystnpd0_run2_alice = HepDataHandeler('/Users/spolitan/cernbox/flowRun3/input/npd0_run2/HEPData-ins2681666-v1-root.root',
                                                 'Table 1',
                                                 '1',
                                                  stat='e1',
                                                  syst='e2',
                                                 )
    gnpd0_run2_cms, gsystnpd0_run2_cms = HepDataHandeler('/Users/spolitan/cernbox/flowRun3/input/npd0_run2/HEPData-ins2610495-v1-root.root',
                                                 'Figure 3a',
                                                 '3',
                                                  syst='e2',
                                                 )    
    SetObjectStyle(gist_npd0_pass4_3050, linecolor=ROOT.kAzure+2, markercolor=ROOT.kAzure+2, markersize=1.6)
    SetObjectStyle(gist_npd0_pass4_syst_3050, linecolor=ROOT.kAzure+2, markercolor=ROOT.kAzure+2, markersize=1.6)
    SetObjectStyle(gnpd0_run2_alice, fillalpha=0., linecolor=ROOT.kGray+1, markercolor=ROOT.kGray+1, markerstyle=ROOT.kFullSquare, fillstyle=1001, markersize=1.6)
    SetObjectStyle(gsystnpd0_run2_alice, fillalpha=0., linecolor=ROOT.kGray+1, markercolor=ROOT.kGray+1, markerstyle=ROOT.kFullSquare, fillstyle=1001, markersize=1.6)
    SetObjectStyle(gnpd0_run2_cms, fillalpha=0., linecolor=ROOT.kAzure-9, markercolor=ROOT.kAzure-9, markerstyle=ROOT.kFullCrossX, fillstyle=1001, markersize=2)
    SetObjectStyle(gsystnpd0_run2_cms, fillalpha=0., linecolor=ROOT.kAzure-9, markercolor=ROOT.kAzure-9, markerstyle=ROOT.kFullCrossX, fillstyle=1001, markersize=2)

    legrun2 = GetLegend(header='Non-prompt D^{0}, #sqrt{#it{s}_{NN}} = 5.02 TeV',
                         xmax=0.5, ncolumns=1, ymin=0.54, ymax=0.68, textsize=0.03)
    legempty_run2 = legrun2.Clone('')
    legempty_run2.SetHeader(' ')
    #legrun2.AddEntry(gnpd0_run2_cms, 'CMS PLB 850 (2024) 138389', 'p')
    leg = GetLegend(header='Non-prompt D^{0}', xmax=0.5, ncolumns=1, ymin=0.58, textsize=0.04, xmin=0.2)
    legempty = leg.Clone('')
    legempty.SetHeader(' ')
    leg.AddEntry(gnpd0_run2_alice, 'ALICE #sqrt{#it{s}_{NN}} = 5.02 TeV |#Delta#it{#eta}| > 0.9', 'p') # #scale[0.8]{(EPJC 83 (2023) 12 1123)}
    leg.AddEntry(gist_npd0_pass4_3050, 'ALICE #sqrt{#it{s}_{NN}} = 5.36 TeV |#Delta#it{#eta}| > 1.3', 'p')
    
    DrawStatSystEmpty(gnpd0_run2_alice, gsystnpd0_run2_alice, False, legempty, markersize=1.6)
    #DrawStatSystEmpty(gnpd0_run2_cms, gsystnpd0_run2_cms, False, legempty_run2, markersize=2)
    DrawStatSystEmpty(gist_npd0_pass4_3050, gist_npd0_pass4_syst_3050, False, legempty, markersize=1.6)
    leg.Draw('same')
    legempty.Draw('same')
    
    latexlarge.DrawLatexNDC(0.20, 0.84, 'ALICE Preliminary')
    latexmedium.DrawLatexNDC(0.20, 0.78, 'Pb#minusPb, 30#minus50%')
    latexdetail.DrawLatexNDC(0.8, 0.86, '|#it{y}| < 0.8')
    cRun3D0vsModels.Update()
    SaveCanvas(cRun3D0vsModels, f'{outdir}npd0_run3_vs_run2', suffix)
    input()  

#_____________________________________________________________________________________
# np D0 meson v2 vsmodels
if plot_npd0_vsmodel:
    cRun3D0vsModels, hframe = GetCanvas('cRun3D0vsModels', axisname, ymin=-0.1, ymax=0.3)

    gtamuv2npD0 = GetPrediction('/Users/spolitan/alice/DmesonAnalysis/models/tamu/NonPromptD0_TAMU_v2_5TeV_3050.txt',
                            1, 148, 0.1, False)
    glidov2npD0 = GetPrediction('/Users/spolitan/cernbox/flowRun3/input/npD0_models/theory-driven/Lido-B2D/vn/cen-30-50/b2D-vn.dat',
                            2, 26, 1, True, model='lido')
    glidov2npD0.RemovePoint(0)
    glidov2npD0.RemovePoint(0)
    glgrv2npD0 = GetPrediction('/Users/spolitan/cernbox/flowRun3/input/npD0_models/theory-driven/LGR/V2_NonPrompt_D0_C30_50_PbPb5020GeV.dat',
                            225, 2025, 0.01, True, model='lgr')
    glgrv2npD0.RemovePoint(0)
    glgrv2npD0.RemovePoint(0)

    gtamuv2npD0.SetLineColor(ROOT.kAzure-9)
    gtamuv2npD0.SetLineStyle(1)
    gtamuv2npD0.SetLineWidth(3)
    gtamuv2npD0.SetFillColorAlpha(ROOT.kAzure-4, 0.2)
    glidov2npD0.SetLineColor(ROOT.kRed-4)
    glidov2npD0.SetLineStyle(1)
    glidov2npD0.SetLineWidth(0)
    glidov2npD0.SetFillColorAlpha(ROOT.kRed-4, 0.2)
    glgrv2npD0.SetLineColor(ROOT.kGreen+2)
    glgrv2npD0.SetLineStyle(1)
    glgrv2npD0.SetLineWidth(0)
    glgrv2npD0.SetFillColorAlpha(ROOT.kGreen+2, 0.2)

    legmodels = GetLegend(ncolumns=1, header='Transport models, #sqrt{#it{s}_{NN}} = 5.02 TeV', xmax=0.5, ymin=0.52, ymax=0.7)
    legmodels.AddEntry(gtamuv2npD0, 'TAMU', 'l')
    legmodels.AddEntry(glidov2npD0, 'LIDO', 'f')
    legmodels.AddEntry(glgrv2npD0, 'LGR', 'f')

    _ = DrawLineAt0(0., 40)
    gtamuv2npD0.Draw('3c same')
    glidov2npD0.Draw('3c same')
    glgrv2npD0.Draw('3c same')
    SetObjectStyle(gist_npd0_pass4_3050, linecolor=ROOT.kBlack, markercolor=ROOT.kBlack)
    SetObjectStyle(gist_npd0_pass4_syst_3050, linecolor=ROOT.kBlack, markercolor=ROOT.kBlack)
    leg = GetLegend(header='', xmax=0.5, ncolumns=1, ymin=0.73)
    legempty = leg.Clone('')
    legempty.SetHeader(' ')
    leg.AddEntry(gist_npd0_pass4_3050, 'Non-prompt D^{0}', 'p')
    DrawStatSystEmpty(gist_npd0_pass4_3050, gist_npd0_pass4_syst_3050, False, legempty, markersize=1)
    leg.Draw('same')
    legempty.Draw('same')
    legmodels.Draw('same')
    
    latex.DrawLatexNDC(0.20, 0.86, 'ALICE Preliminary')
    latexdetail.DrawLatexNDC(0.20, 0.8, 'Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV, 30#minus50%')
    latexdetail.DrawLatexNDC(0.8, 0.86, '|#it{y}| < 0.8')
    cRun3D0vsModels.Update()
    SaveCanvas(cRun3D0vsModels, f'{outdir}d0_run3_vs_model', suffix)
    input()

#_____________________________________________________________________________________
# D meson v2 in Run 3 pass4 3050
if plot_d_v2_run3_3050:
    leg = GetLegend(header='30#minus50%, |#it{y}| < 0.8, |#Delta#it{#eta}| > 1.3', xmax=0.5, ncolumns=1, ymin=0.56, textsize=0.033, xmin=0.186, ymax=0.8)
    legempty = leg.Clone('legempty')
    legempty.SetHeader(' ')
    leg.AddEntry(gist_d0_pass4_3050, 'Prompt D^{0}', 'p')
    leg.AddEntry(gist_dp_pass4_3040, 'Prompt D^{+}', 'p')
    leg.AddEntry(gist_ds_pass4_3040, 'Prompt D_{s}^{+}', 'p')
    leg.AddEntry(gist_jpsi_pass4_3050, 'Inclusive J/#psi 2.5 < #it{y} < 4, |#Delta#it{#eta}| > 1.5', 'p')
    legrun2 = GetLegend(header='#sqrt{#it{s}_{NN}} = 5.02 TeV, 30#minus40%', xmax=0.85, ncolumns=1, ymin=0.64, xmin=0.58, textsize=0.033, ymax=0.8)
    legrun2empty = legrun2.Clone('legempty')
    legrun2empty.SetHeader(' ')
    legrun2.AddEntry(gh_run2_alice_3040, '#pi^{+} |#it{y}| < 0.5,', 'p')
    legrun2.AddEntry(gh_run2_alice_3040, '|#Delta#it{#eta}| > 2.0', '')
    legrun2.AddEntry(gh_run2_alice_3040, 'JHEP 09 (2018) 006', '')

    cDv2run3, hframe = GetCanvas('cDv2run3', axisnamevareta, ymax=0.5, ymin=-0, xmax=25) 
    cDv2run3.cd(1).SetTopMargin(0.02)
    DrawStatSystEmpty(gh_run2_alice_3040, ghsyst_run2_alice_3040, False, legrun2empty, markersize=1.2)
    DrawStatSystEmpty(gist_d0_pass4_3050, gist_d0_pass4_syst_3050, False, legempty, markersize=1.2)
    DrawStatSystEmpty(gist_dp_pass4_3050, gist_dp_pass4_syst_3050, False, legempty)
    DrawStatSystEmpty(gist_ds_pass4_3050, gist_ds_pass4_syst_3050, False, legempty)
    DrawStatSystEmpty(gist_jpsi_pass4_3050, gist_jpsi_pass4_syst_3050, False, legempty, markersize=1)

    leg.Draw('same')
    legempty.Draw('same')
    legrun2.Draw('same')
    legrun2empty.AddEntry(gh_run2_alice_3040, ' ', '')
    legrun2empty.AddEntry(gh_run2_alice_3040, ' ', '')
    legrun2empty.Draw('same')
    
    latexmedium.DrawLatexNDC(0.2, 0.89, 'ALICE Preliminary')
    latex.DrawLatexNDC(0.2, 0.83, 'Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV')
    cDv2run3.Update()
    SaveCanvas(cDv2run3, f'{outdir}Dmeson_v2_run3_3050', suffix)
    input()

#_____________________________________________________________________________________
# Charm v2 in Run 3 pass4 3050
if plot_charm_3050:
    leg = GetLegend(header='', xmax=0.6, ncolumns=2, ymin=0.60, textsize=0.033, xmin=0.196, ymax=0.78)
    legempty = GetLegend(header='', xmax=0.9, ncolumns=4, ymin=0.74, textsize=0.033, xmin=0.196, ymax=0.88)

    leg.AddEntry(gist_lc_pass4_3050, 'Prompt #Lambda_{c}^{+}', 'p')
    leg.AddEntry(gist_d0_pass4_3050, 'Prompt D^{0}', 'p')
    leg.AddEntry(gist_dp_pass4_3040, 'Prompt D^{+}', 'p')
    leg.AddEntry(gist_ds_pass4_3040, 'Prompt D_{s}^{+}', 'p')
    #leg.AddEntry(gist_jpsi_pass4_3050, 'Inclusive J/#psi', 'p')

    gist_lc_pass4_3050_empty = gist_lc_pass4_3050.Clone('empty')
    gist_d0_pass4_3050_empty = gist_d0_pass4_3050.Clone('empty')
    gist_dp_pass4_3050_empty = gist_dp_pass4_3050.Clone('empty')
    gist_ds_pass4_3050_empty = gist_ds_pass4_3050.Clone('empty')
    #gist_jpsi_pass4_3050_empty = gist_jpsi_pass4_3050.Clone('empty')

    SetObjectStyle(gist_lc_pass4_3050_empty, markercolor=ROOT.kBlack, linewidth=0, markerstyle=ROOT.kOpenCrossX, markersize=1.8)
    SetObjectStyle(gist_d0_pass4_3050_empty, markercolor=ROOT.kBlack, linewidth=0, markerstyle=ROOT.kOpenCircle)
    SetObjectStyle(gist_dp_pass4_3050_empty, markercolor=ROOT.kBlack, linewidth=0, markerstyle=ROOT.kOpenDiamond)
    SetObjectStyle(gist_ds_pass4_3050_empty, markercolor=ROOT.kBlack, linewidth=0, markerstyle=ROOT.kOpenCross, markersize=1.8)
    #ßSetObjectStyle(gist_jpsi_pass4_3050_empty, markercolor=ROOT.kBlack, linewidth=0, markerstyle=ROOT.kOpenSquare)


    gtamuv2D0 = GetPrediction('/home/spolitan/alice/hcv2-prl-figures/ste/models/tamu/PromptD_TAMU_v2_5TeV_3050.txt',
                            1, 10, 1, True)
    #gtamuv2Ds = GetPrediction('/Users/spolitan/alice/DmesonAnalysis/models/tamu/PromptDs_TAMU_v2_5TeV_3050.txt',
    #                        1, 10, 1, True)
    #gtamuv2Jpsi = GetPrediction('/Users/spolitan/cernbox/flowRun3/input/tamu_jpsi_3050.txt',
    #                        35, 110, 0.1, False)
    gtamuv2lcup = GetPrediction('/home/spolitan/alice/hcv2-prl-figures/ste/lc-up.dat',
                            11, 110, 0.1, False)
    gtamuv2lclow = GetPrediction('/home/spolitan/alice/hcv2-prl-figures/ste/lc-low.dat',
                            12, 110, 0.1, False)
    gtamuv2lc = compute_central_graph(gtamuv2lcup, gtamuv2lclow)

    #gtamuv2Ds.SetLineColor(ROOT.kAzure)
    #gtamuv2Ds.SetLineStyle(1)
    #gtamuv2Ds.SetLineWidth(0)
    #gtamuv2Ds.SetFillColorAlpha(ROOT.kAzure, 0.2)
    #gtamuv2D0.SetFillStyle(3145)
    gtamuv2D0.SetLineColor(ROOT.kRed)
    gtamuv2D0.SetLineStyle(1)
    gtamuv2D0.SetLineWidth(0)
    gtamuv2D0.SetFillColorAlpha(ROOT.kRed, 0.2)
    #gtamuv2D0.SetFillStyle(3145)
    #ßgtamuv2Jpsi.SetLineColor(ROOT.kOrange-3)
    #ßgtamuv2Jpsi.SetLineStyle(9)
    #ßgtamuv2Jpsi.SetLineWidth(3)
    #ßgtamuv2Jpsi.SetFillColorAlpha(ROOT.kOrange-3, 0.2)
    gtamuv2lcup.SetLineColor(ROOT.kGreen+2)
    gtamuv2lcup.SetLineStyle(9)
    gtamuv2lcup.SetLineWidth(3)
    gtamuv2lcup.SetFillColorAlpha(ROOT.kGreen+2, 0.2)
    gtamuv2lclow.SetLineColor(ROOT.kGreen+2)
    gtamuv2lclow.SetLineStyle(9)
    gtamuv2lclow.SetLineWidth(3)
    gtamuv2lc.SetFillColorAlpha(ROOT.kGreen+2, 0.2)
    gtamuv2lc.SetLineColor(ROOT.kGreen+2)
    gtamuv2lc.SetLineStyle(9)
    gtamuv2lc.SetLineWidth(0)
    gtamuv2lc.SetFillColorAlpha(ROOT.kGreen+2, 0.2)
    #gtamuv2Ds.SetFillStyle(3145)
    legmodels = GetLegend(header='TAMU', xmin=0.65, ncolumns=2, ymin=0.65, textsize=0.033, xmax=0.96, ymax=0.8)
    legmodels.AddEntry(gtamuv2lc, '#Lambda_{c}^{+}', 'f')
    legmodels.AddEntry(gtamuv2D0, 'D^{0}', 'f')
    #legmodels.AddEntry(gtamuv2Ds, 'D_{s}^{+}', 'f')
    #legmodels.AddEntry(gtamuv2Jpsi, 'J/#psi', 'l')

    
    cDv2run3, hframe = GetCanvas('cDv2run3', axisnamevareta, ymax=0.40, ymin=-0.04, xmax=26)
    cDv2run3.SetTopMargin(0.05)
    hframe.GetYaxis().SetTitleOffset(1.4)
    hframe.GetYaxis().SetLabelSize(0.04)
    hframe.GetXaxis().SetLabelSize(0.04)
    gtamuv2lc.Draw('3c same')
    gtamuv2D0.Draw('3c same')
    gtamuv2Ds.Draw('3c same')
    gtamuv2Jpsi.Draw('3c same')
    DrawStatSystEmpty(gist_lc_pass4_3050, gist_lc_pass4_syst_3050, False, False)
    DrawStatSystEmpty(gist_d0_pass4_3050, gist_d0_pass4_syst_3050, False, False)
    DrawStatSystEmpty(gist_dp_pass4_3050, gist_dp_pass4_syst_3050, False, False)
    DrawStatSystEmpty(gist_ds_pass4_3050, gist_ds_pass4_syst_3050, False, False)
    DrawStatSystEmpty(gist_jpsi_pass4_3050, gist_jpsi_pass4_syst_3050, False, False, markersize=1)

    leg.Draw('same')
    #legempty.Draw('same')
    legmodels.Draw('same')
    line = DrawLineAt0(0., 26)

    latex.DrawLatexNDC(0.20, 0.88, 'ALICE Preliminary')
    latexdetail2.DrawLatexNDC(0.20, 0.81, 'Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV, 30#minus50%')
    latexsmall.DrawLatexNDC(0.20, 0.18, 'D, #Lambda_{c}^{+}: |#it{y}| < 0.8 |#Delta#it{#eta}| > 1.3, J/#psi: 2.5 < #it{y} < 4 |#Delta#it{#eta}| > 1.5')
    #latexsmall.DrawLatexNDC(0.79, 0.765, '2.5 < #it{y} < 4')
    cDv2run3.Update()
    SaveCanvas(cDv2run3, f'{outdir}charm_v2_3050', suffix)
    input()

    # no log
    cDv2run3, hframe = GetCanvas('cDv2run3', axisnamevareta, ymax=0.40, ymin=-0.04, xmax=26, setlogx=False)
    cDv2run3.SetTopMargin(0.05)
    hframe.GetYaxis().SetTitleOffset(1.4)
    hframe.GetYaxis().SetLabelSize(0.04)
    hframe.GetXaxis().SetLabelSize(0.04)
    gtamuv2lc.Draw('3c same')
    gtamuv2D0.Draw('3c same')
    gtamuv2Ds.Draw('3c same')
    gtamuv2Jpsi.Draw('3c same')
    DrawStatSystEmpty(gist_lc_pass4_3050, gist_lc_pass4_syst_3050, False, False)
    DrawStatSystEmpty(gist_d0_pass4_3050, gist_d0_pass4_syst_3050, False, False)
    DrawStatSystEmpty(gist_dp_pass4_3050, gist_dp_pass4_syst_3050, False, False)
    DrawStatSystEmpty(gist_ds_pass4_3050, gist_ds_pass4_syst_3050, False, False)
    DrawStatSystEmpty(gist_jpsi_pass4_3050, gist_jpsi_pass4_syst_3050, False, False, markersize=1)

    leg.Draw('same')
    #legempty.Draw('same')
    legmodels.Draw('same')
    line = DrawLineAt0(0.4, 26)

    latex.DrawLatexNDC(0.20, 0.88, 'ALICE Preliminary')
    latexdetail2.DrawLatexNDC(0.20, 0.81, 'Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV, 30#minus50%')
    latexsmall.DrawLatexNDC(0.20, 0.18, 'D, #Lambda_{c}^{+}: |#it{y}| < 0.8 |#Delta#it{#eta}| > 1.3, J/#psi: 2.5 < #it{y} < 4 |#Delta#it{#eta}| > 1.5')
    #latexsmall.DrawLatexNDC(0.79, 0.765, '2.5 < #it{y} < 4')
    cDv2run3.Update()
    SaveCanvas(cDv2run3, f'{outdir}charm_v2_3050_nolog', suffix)
    input()
    

    # d0 lc
    leg = GetLegend(header='#it{v}_{2}{SP, |#Delta#it{#eta}| > 1.3}, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV, 30#minus50%', xmax=0.6, ncolumns=2, ymin=0.75, textsize=0.028, xmin=0.196, ymax=0.85)
    legempty = GetLegend(header='', xmax=0.9, ncolumns=4, ymin=0.74, textsize=0.035, xmin=0.196, ymax=0.8)

    leg.AddEntry(gist_d0_pass4_3050, 'Prompt D^{0}', 'p')
    leg.AddEntry(gist_lc_pass4_3050, 'Prompt #Lambda_{c}^{+}', 'p')

    legLF = GetLegend(header='#it{v}_{2}{2, |#Delta#it{#eta}| > 0.8}, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV, 30#minus40%',  xmax=0.55, ncolumns=2, ymin=0.63, textsize=0.028, xmin=0.196, ymax=0.73)
    legLF.AddEntry(gkzero, 'K_{S}^{0}', 'p')
    legLF.AddEntry(glambda, '#Lambda', 'p')
    #leg.AddEntry(gist_dp_pass4_3040, 'D^{+}', 'p')
    #ßleg.AddEntry(gist_ds_pass4_3040, 'Prompt D_{s}^{+}', 'p')
    #ßleg.AddEntry(gist_jpsi_pass4_3050, 'Inclusive J/#psi', 'p')

    gist_lc_pass4_3050_empty = gist_lc_pass4_3050.Clone('empty')
    gist_d0_pass4_3050_empty = gist_d0_pass4_3050.Clone('empty')
    gist_dp_pass4_3040_empty = gist_dp_pass4_3040.Clone('empty')
    gist_ds_pass4_3050_empty = gist_ds_pass4_3050.Clone('empty')
    #gist_jpsi_pass4_3050_empty = gist_jpsi_pass4_3050.Clone('empty')

    SetObjectStyle(gist_lc_pass4_3050_empty, markercolor=ROOT.kBlack, linewidth=0, markerstyle=ROOT.kOpenCrossX, markersize=1.8)
    SetObjectStyle(gist_d0_pass4_3050_empty, markercolor=ROOT.kBlack, linewidth=0, markerstyle=ROOT.kOpenCircle)
    SetObjectStyle(gist_dp_pass4_3040_empty, markercolor=ROOT.kBlack, linewidth=0, markerstyle=ROOT.kOpenDiamond)
    SetObjectStyle(gist_ds_pass4_3050_empty, markercolor=ROOT.kBlack, linewidth=0, markerstyle=ROOT.kOpenCross, markersize=1.8)
    #ßSetObjectStyle(gist_jpsi_pass4_3050_empty, markercolor=ROOT.kBlack, linewidth=0, markerstyle=ROOT.kOpenSquare)

    legempty.AddEntry(gist_lc_pass4_3050_empty, '#Lambda_{c}^{+}', 'p')
    legempty.AddEntry(gist_d0_pass4_3050_empty, 'D^{0}', 'p')
    #legempty.AddEntry(gist_dp_pass4_3040_empty, ' ', 'p')

    cDv2run3, hframe = GetCanvas('cDv2run3', axisnamevareta, ymax=0.50, ymin=-0.04, xmax=26)
    cDv2run3.SetTopMargin(0.05)
    hframe.GetYaxis().SetTitleOffset(1.4)
    hframe.GetYaxis().SetLabelSize(0.04)
    hframe.GetXaxis().SetLabelSize(0.04)
    #gtamuv2lc.Draw('3c same')
    #gtamuv2D0.Draw('3c same')
    #gtamuv2Ds.Draw('3c same')
    #gtamuv2Jpsi.Draw('3c same')
    glambda.Draw('EPZ same')
    gkzero.Draw('EPZ same')
    glambdaempty.Draw('PZ same')
    gkzeroempty.Draw('PZ same')
    DrawStatSystEmpty(gist_lc_pass4_3050, gist_lc_pass4_syst_3050, False, False)
    DrawStatSystEmpty(gist_d0_pass4_3050, gist_d0_pass4_syst_3050, False, False)
    #DrawStatSystEmpty(gist_dp_pass4_3040, gist_dp_pass4_syst_3050, False, False)
    #DrawStatSystEmpty(gist_ds_pass4_3050, gist_ds_pass4_syst_3050, False, False)
    #DrawStatSystEmpty(gist_jpsi_pass4_3050, gist_jpsi_pass4_syst_3050, False, False, markersize=1)

    leg.Draw('same')
    legLF.Draw('same')
    # Dirty trick to create a white box over the legend
    x1 = 0.5
    x2 = 4.78   
    y1 = 0.23
    y2 = 0.26
    box = ROOT.TBox(x1, y1, x2, y2)
    box.SetFillColor(ROOT.kWhite)  # White box
    box.SetLineColor(ROOT.kWhite)  # Hide border
    box.Draw()
    #legempty.Draw('same')
    #legmodels.Draw('same')
    line = DrawLineAt0(0.4, 26)

    latex.DrawLatexNDC(0.20, 0.88, 'ALICE')
    #latexsmall.DrawLatexNDC(0.79, 0.765, '2.5 < #it{y} < 4')
    cDv2run3.Update()
    SaveCanvas(cDv2run3, f'{outdir}charm_v2_3050_d0_lc', suffix)
    input()

    # no log d0 lc
    
    leg = GetLegend(header='', xmax=0.6, ncolumns=2, ymin=0.65, textsize=0.033, xmin=0.196, ymax=0.8)
    legempty = GetLegend(header='', xmax=0.9, ncolumns=4, ymin=0.74, textsize=0.035, xmin=0.196, ymax=0.8)

    leg.AddEntry(gist_lc_pass4_3050, 'Prompt #Lambda_{c}^{+}', 'p')
    leg.AddEntry(gist_d0_pass4_3050, 'Prompt D^{0}', 'p')
    #leg.AddEntry(gist_dp_pass4_3040, 'D^{+}', 'p')
    leg.AddEntry(gist_ds_pass4_3040, 'Prompt D_{s}^{+}', 'p')
    leg.AddEntry(gist_jpsi_pass4_3050, 'Inclusive J/#psi', 'p')
    

    gist_lc_pass4_3050_empty = gist_lc_pass4_3050.Clone('empty')
    gist_d0_pass4_3050_empty = gist_d0_pass4_3050.Clone('empty')
    gist_dp_pass4_3040_empty = gist_dp_pass4_3040.Clone('empty')
    gist_ds_pass4_3050_empty = gist_ds_pass4_3050.Clone('empty')
    gist_jpsi_pass4_3050_empty = gist_jpsi_pass4_3050.Clone('empty')

    SetObjectStyle(gist_lc_pass4_3050_empty, markercolor=ROOT.kBlack, linewidth=0, markerstyle=ROOT.kOpenCrossX, markersize=1.8)
    SetObjectStyle(gist_d0_pass4_3050_empty, markercolor=ROOT.kBlack, linewidth=0, markerstyle=ROOT.kOpenCircle)
    SetObjectStyle(gist_dp_pass4_3040_empty, markercolor=ROOT.kBlack, linewidth=0, markerstyle=ROOT.kOpenDiamond)
    SetObjectStyle(gist_ds_pass4_3050_empty, markercolor=ROOT.kBlack, linewidth=0, markerstyle=ROOT.kOpenCross, markersize=1.8)
    SetObjectStyle(gist_jpsi_pass4_3050_empty, markercolor=ROOT.kBlack, linewidth=0, markerstyle=ROOT.kOpenSquare)

    legempty.AddEntry(gist_lc_pass4_3050_empty, '#Lambda_{c}^{+}', 'p')
    legempty.AddEntry(gist_d0_pass4_3050_empty, 'D^{0}', 'p')
    #legempty.AddEntry(gist_dp_pass4_3040_empty, ' ', 'p')

    cDv2run3, hframe = GetCanvas('cDv2run3', axisnamevareta, ymax=0.40, ymin=-0.04, xmax=26, setlogx=False)
    cDv2run3.SetTopMargin(0.05)
    hframe.GetYaxis().SetTitleOffset(1.4)
    hframe.GetYaxis().SetLabelSize(0.04)
    hframe.GetXaxis().SetLabelSize(0.04)
    #gtamuv2lc.Draw('3c same')
    #gtamuv2D0.Draw('3c same')
    #gtamuv2Ds.Draw('3c same')
    #gtamuv2Jpsi.Draw('3c same')
    DrawStatSystEmpty(gist_lc_pass4_3050, gist_lc_pass4_syst_3050, False, False)
    DrawStatSystEmpty(gist_d0_pass4_3050, gist_d0_pass4_syst_3050, False, False)
    #DrawStatSystEmpty(gist_dp_pass4_3040, gist_dp_pass4_syst_3050, False, False)
    #DrawStatSystEmpty(gist_ds_pass4_3050, gist_ds_pass4_syst_3050, False, False)
    #DrawStatSystEmpty(gist_jpsi_pass4_3050, gist_jpsi_pass4_syst_3050, False, False, markersize=1)

    leg.Draw('same')
    # Dirty trick to create a white box over the legend
    x1 = 0.5
    x2 = 15 
    y1 = 0.237
    y2 = 0.259
    box = ROOT.TBox(x1, y1, x2, y2)
    box.SetFillColor(ROOT.kWhite)  # White box
    box.SetLineColor(ROOT.kWhite)  # Hide border
    box.Draw()
    #legempty.Draw('same')
    #legmodels.Draw('same')
    line = DrawLineAt0(0.4, 26)

    latex.DrawLatexNDC(0.20, 0.88, 'ALICE Preliminary')
    latexdetail2.DrawLatexNDC(0.20, 0.81, 'Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV, 30#minus50%')
    latexsmall.DrawLatexNDC(0.20, 0.18, 'D^{0}, #Lambda_{c}^{+}: |#it{y}| < 0.8 |#Delta#it{#eta}| > 1.3')
    #latexsmall.DrawLatexNDC(0.79, 0.765, '2.5 < #it{y} < 4')
    cDv2run3.Update()
    SaveCanvas(cDv2run3, f'{outdir}charm_v2_3050_nolog_d0_lc', suffix)

    input()

#_____________________________________________________________________________________
# Ds vs D0 meson v2 in Run 3 pass4 3050
if plot_d0_ds_v2_run3_3050:
    leg = GetLegend(header='Prompt D', xmax=0.5, ncolumns=1, ymin=0.6, xmin=0.195, ymax=0.76)
    legempty = leg.Clone('legempty')
    legempty.SetHeader(' ')
    leg.AddEntry(gist_d0_pass4_3050, 'D^{0}', 'p')
    leg.AddEntry(gist_ds_pass4_3040, 'D_{s}^{+}', 'p')

    gtamuv2D0 = GetPrediction('/Users/spolitan/alice/DmesonAnalysis/models/tamu/PromptD_TAMU_v2_5TeV_3050.txt',
                            1, 10, 1, True)
    gtamuv2Ds = GetPrediction('/Users/spolitan/alice/DmesonAnalysis/models/tamu/PromptDs_TAMU_v2_5TeV_3050.txt',
                            1, 10, 1, True)
    gdiff, gdiffsyst = compute_ds_d0_v2_diff(gist_ds_pass4_3050, gist_d0_pass4_3050_wider, gist_ds_pass4_syst_3050, gist_d0_pass4_syst_3050_wider)
    gtamuv2Ds.SetLineColor(ROOT.kAzure)
    gtamuv2Ds.SetLineStyle(1)
    gtamuv2Ds.SetLineWidth(0)
    gtamuv2Ds.SetFillColorAlpha(ROOT.kAzure, 0.2)
    #gtamuv2D0.SetFillStyle(3145)
    gtamuv2D0.SetLineColor(ROOT.kRed)
    gtamuv2D0.SetLineStyle(1)
    gtamuv2D0.SetLineWidth(0)
    gtamuv2D0.SetFillColorAlpha(ROOT.kRed, 0.2)
    #gtamuv2Ds.SetFillStyle(3145)
    legmodels = GetLegend(header='TAMU #sqrt{#it{s}_{NN}} = 5.02 TeV', xmax=0.58, ncolumns=1, ymin=0.6, xmin=0.38, ymax=0.76)
    legmodels.AddEntry(gtamuv2D0, 'D^{0}', 'f')
    legmodels.AddEntry(gtamuv2Ds, 'D_{s}^{+}', 'f')

    cDv2run3, hframe = GetCanvas('cDv2run3', axisname, xmax=25, ymin=0, ymax=0.3)
    gtamuv2D0.Draw('3c same')
    gtamuv2Ds.Draw('3c same')
    DrawStatSystEmpty(gist_d0_pass4_3050_wider, gist_d0_pass4_syst_3050_wider, False, legempty, markersize=1.2)
    DrawStatSystEmpty(gist_ds_pass4_3050, gist_ds_pass4_syst_3050, False, legempty)
    latex.DrawLatexNDC(0.20, 0.86, 'ALICE Preliminary')
    latexdetail.DrawLatexNDC(0.20, 0.8, 'Pb#minusPb, 30#minus50%, #sqrt{#it{s}_{NN}} = 5.36 TeV')
    latexdetail.DrawLatexNDC(0.8, 0.86, '|#it{y}| < 0.8')
    leg.Draw('same')
    legempty.Draw('same')
    legmodels.Draw('same')

    cDv2run3.Update()
    SaveCanvas(cDv2run3, f'{outdir}Ds_vs_D0_v2_run3', suffix)
    input()

#_____________________________________________________________________________________
# Ds vs D0 meson v2 in Run 3 pass4 3050
if plot_d0_ds_v2_run3_3050_wdiff:
    leg = GetLegend(header='Prompt D', xmax=0.5, ncolumns=1, ymin=0.54, xmin=0.22, ymax=0.7)
    legempty = leg.Clone('legempty')
    legempty.SetHeader(' ')
    leg.AddEntry(gist_d0_pass4_3050, 'D^{0}', 'p')
    leg.AddEntry(gist_ds_pass4_3040, 'D_{s}^{+}', 'p')

    legdiff = GetLegend(header=' ', xmax=0.5, ncolumns=1, ymin=0.8, xmin=0.22, ymax=0.9)

    gtamuv2D0 = GetPrediction('/Users/spolitan/alice/DmesonAnalysis/models/tamu/PromptD_TAMU_v2_5TeV_3050.txt',
                            1, 10, 1, True)
    gtamuv2Ds = GetPrediction('/Users/spolitan/alice/DmesonAnalysis/models/tamu/PromptDs_TAMU_v2_5TeV_3050.txt',
                            1, 10, 1, True)
    gdiff, gdiffsyst = compute_ds_d0_v2_diff(gist_ds_pass4_3050, gist_d0_pass4_3050_wider, gist_ds_pass4_syst_3050, gist_d0_pass4_syst_3050_wider)
    gtamuv2Ds.SetLineColor(ROOT.kAzure)
    gtamuv2Ds.SetLineStyle(1)
    gtamuv2Ds.SetLineWidth(0)
    gtamuv2Ds.SetFillColorAlpha(ROOT.kAzure, 0.2)
    #gtamuv2D0.SetFillStyle(3145)
    gtamuv2D0.SetLineColor(ROOT.kRed)
    gtamuv2D0.SetLineStyle(1)
    gtamuv2D0.SetLineWidth(0)
    gtamuv2D0.SetFillColorAlpha(ROOT.kRed, 0.2)
    #gtamuv2Ds.SetFillStyle(3145)
    legmodels = GetLegend(header='TAMU #sqrt{#it{s}_{NN}} = 5.02 TeV', xmax=0.58, ncolumns=1, ymin=0.54, xmin=0.38, ymax=0.7)
    legmodels.AddEntry(gtamuv2D0, 'D^{0}', 'f')
    legmodels.AddEntry(gtamuv2Ds, 'D_{s}^{+}', 'f')

    cD0v2run3, _ = GetCanvas2sub('cDv2run3', 
                                     1,
                                     10,
                                     0.001,
                                     0.4,
                                     -0.082,
                                     0.064,
                                     axisname,
                                     ';#it{p}_{T} (GeV/#it{c}); #it{v}_{2}(D_{s}^{+}) - #it{v}_{2}(D^{0})'
                                     )
    cD0v2run3.cd(1)
    gtamuv2D0.Draw('3c same')
    gtamuv2Ds.Draw('3c same')
    DrawStatSystEmpty(gist_d0_pass4_3050_wider, gist_d0_pass4_syst_3050_wider, False, legempty, markersize=1.2)
    DrawStatSystEmpty(gist_ds_pass4_3050, gist_ds_pass4_syst_3050, False, legempty)
    latexlarge.DrawLatexNDC(0.22, 0.82, 'ALICE Preliminary')
    latexmedium.DrawLatexNDC(0.22, 0.76, 'Pb#minusPb, 30#minus50%, #sqrt{#it{s}_{NN}} = 5.36 TeV')
    latexmedium.DrawLatexNDC(0.8, 0.82, '|#it{y}| < 0.8')
    leg.Draw('same')
    legempty.Draw('same')
    legmodels.Draw('same')
    cD0v2run3.cd(2)
    legdiff.AddEntry(gdiff, '#sqrt{stat.^{2} + syst.^{2}}', 'l')
    SetObjectStyle(gdiff, markercolor=ROOT.kBlack, linecolor=ROOT.kBlack, markerstyle=ROOT.kFullCircle,linewidth=3, linestyle=9)
    _ = DrawLineAt0(1, 10)
    legdiff.Draw()
    gdiff.Draw('same pe')
    gdiff.Draw('same ||')

    cD0v2run3.Update()
    SaveCanvas(cD0v2run3, f'{outdir}Ds_vs_D0_v2_run3_diff', suffix)
    input()

#_____________________________________________________________________________________
# D meson v2 in Run 3 pass4 3040/40500/6080
if plot_d_v2_run3_vscent:
    leg = GetLegend(header='#sqrt{#it{s}_{NN}} = 5.36 TeV, |#it{y}| < 0.8, |#Delta#it{#eta}| > 1.3', xmax=0.5, xmin=0.215, ncolumns=1, ymin=0.64, ymax=0.84, textsize=0.042)
    legempty = leg.Clone('legempty')
    legempty.SetHeader(' ')
    leg.AddEntry(gist_d0_pass4_3040, 'Prompt D^{0}', 'p')
    leg.AddEntry(gist_dp_pass4_3040, 'Prompt D^{+}', 'p')
    leg.AddEntry(gist_ds_pass4_3040, 'Prompt D_{s}^{+}', 'p')
    legrun2 = GetLegend(header='#sqrt{#it{s}_{NN}} = 5.02 TeV, |#it{y}| < 0.5, |#Delta#it{#eta}| > 2.0', xmax=0.3, xmin=0.065, ncolumns=1, ymin=0.82, ymax=0.96, textsize=0.042)
    legrun2empty = legrun2.Clone('legempty')
    legrun2empty.SetHeader(' ')
    legrun2.AddEntry(gh_run2_alice_3040, '#pi^{+} (JHEP 09 (2018) 006)', 'p')

    yrange = [-0.05, 0.41]
    cDv2run3, hframes = GetCanvas3sub('cDv2run3', axisnamevareta)
    cDv2run3.cd(1).SetLogx()
    hframes[0].GetYaxis().SetRangeUser(yrange[0], yrange[1])
    hframes[0].GetXaxis().SetRangeUser(0.48, 26)
    hframes[0].GetXaxis().SetMoreLogLabels()
    line1 = DrawLineAt0(-0., 25)
    line1.Draw()
    leg.Draw('same')
    legempty.Draw('same')
    latexmedium.DrawLatexNDC(0.22, 0.92, 'ALICE Preliminary')
    latexmedium.DrawLatexNDC(0.22, 0.86, 'Pb#minusPb')
    latexdetail2.DrawLatexNDC(0.84, 0.92, '30#minus40%')
    DrawStatSystEmpty(gh_run2_alice_3040, ghsyst_run2_alice_3040, False, legrun2empty, markersize=1.2)
    DrawStatSystEmpty(gist_d0_pass4_3040, gist_d0_pass4_syst_3040, False, legempty, markersize=1.2)
    DrawStatSystEmpty(gist_dp_pass4_3040, gist_dp_pass4_syst_3040, False, legempty)
    DrawStatSystEmpty(gist_ds_pass4_3040, gist_ds_pass4_syst_3040, False, legempty)
    cDv2run3.Update()
    
    cDv2run3.cd(2).SetLogx()
    hframes[1].GetYaxis().SetRangeUser(yrange[0], yrange[1])
    hframes[1].GetXaxis().SetRangeUser(0.51, 26)
    hframes[1].GetXaxis().SetMoreLogLabels()
    legrun2.Draw('same')
    legrun2empty.Draw('same')
    latexdetail2.DrawLatexNDC(0.82, 0.92, '40#minus50%')
    line2 = DrawLineAt0(0.1, 25)
    line2.Draw()

    DrawStatSystEmpty(gh_run2_alice_4050, ghsyst_run2_alice_4050, False, False, markersize=1.2)
    DrawStatSystEmpty(gist_d0_pass4_4050, gist_d0_pass4_syst_4050, False, False, markersize=1.2)
    DrawStatSystEmpty(gist_dp_pass4_4050, gist_dp_pass4_syst_4050, False, False)
    DrawStatSystEmpty(gist_ds_pass4_4050, gist_ds_pass4_syst_4050, False, False)
    
    # Dirty trick to create a white box over the label
    x1 = 0.478
    x2 = 0.665
    y1 = cDv2run3.cd(2).GetUymin()  # Get bottom of the pad
    y2 = cDv2run3.cd(2).GetUymin() + (cDv2run3.cd(2).GetUymax() - cDv2run3.cd(2).GetUymin()) * 0.172  # Small height to cover only label
    box = ROOT.TBox(x1, y1, x2, y2)
    box.SetFillColor(ROOT.kWhite)  # White box
    box.SetLineColor(ROOT.kWhite)  # Hide border
    box.Draw()
    
    cDv2run3.cd(3).SetLogx()
    hframes[2].GetYaxis().SetRangeUser(yrange[0], yrange[1])
    hframes[2].GetXaxis().SetRangeUser(0.22, 12)
    hframes[2].GetXaxis().SetMoreLogLabels()
    latexdetail2.DrawLatexNDC(0.68, 0.92, '#pi^{+}: 60#minus70%')
    latexdetail2.DrawLatexNDC(0.68, 0.86, 'D: 60#minus80%')
    line3 = DrawLineAt0(-0.5, 12)
    line3.Draw()
    DrawStatSystEmpty(gh_run2_alice_6070, ghsyst_run2_alice_6070, False, False, markersize=1.2)
    DrawStatSystEmpty(gist_d0_pass4_6080, gist_d0_pass4_syst_6080, False, False, markersize=1.2)
    DrawStatSystEmpty(gist_dp_pass4_6080, gist_dp_pass4_syst_6080, False, False)
    DrawStatSystEmpty(gist_ds_pass4_6080, gist_ds_pass4_syst_6080, False, False)
    # Dirty trick to create a white box over the label
    x1 = 0.18
    x2 = 0.30
    y1 = -0.085  # Get bottom of the pad
    y2 = -0.055  # Small height to cover only label
    box2 = ROOT.TBox(x1, y1, x2, y2)
    cDv2run3.cd(3)
    box2.SetFillColor(ROOT.kWhite)  # White box
    box2.SetLineColor(ROOT.kWhite)  # Hide border
    box2.Draw()
    
    cDv2run3.Update()
    SaveCanvas(cDv2run3, f'{outdir}Dmeson_v2_run3_vscent', suffix)
    input()

#_____________________________________________________________________________________
# D meson v2 in Run 3 pass4 3040/40500/6080 vs phsd
if plot_d_v2_run3_vscent_vmodel:
    leg = GetLegend(header='Prompt D', xmax=0.4, xmin=0.065, ncolumns=1, ymin=0.74, ymax=0.96, textsize=0.05)
    legempty = leg.Clone('legempty')
    legempty.SetHeader(' ')
    leg.AddEntry(gist_d0_pass4_3040, 'D^{0}', 'p')
    leg.AddEntry(gist_dp_pass4_3040, 'D^{+}', 'p')
    leg.AddEntry(gist_ds_pass4_3040, 'D_{s}^{+}', 'p')
    legmodel = GetLegend(header='PHSD',xmax=0.65, xmin=0.275, ncolumns=1, ymin=0.74, ymax=0.96, textsize=0.05)

    gphsdv2_dp_3040, gphsdv2_dp_4050, gphsdv2_dp_6080 = GetPrediction('/Users/spolitan/cernbox/flowRun3/input/phsd/Dplus.txt',15, 100, 0.1, True, model='phsd')
    gphsdv2_ds_3040, gphsdv2_ds_4050, gphsdv2_ds_6080 = GetPrediction('/Users/spolitan/cernbox/flowRun3/input/phsd/Ds.txt',15, 100, 0.1, True, model='phsd')
    gphsdv2_d0_3040, gphsdv2_d0_4050, gphsdv2_d0_6080 = GetPrediction('/Users/spolitan/cernbox/flowRun3/input/phsd/D0.txt',15, 100, 0.1, True, model='phsd')
    gphsdv2_dp_3040.SetLineColor(ROOT.kGreen+2)
    gphsdv2_dp_3040.SetLineStyle(1)
    gphsdv2_dp_3040.SetLineWidth(0)
    gphsdv2_dp_3040.SetFillColorAlpha(ROOT.kGreen+2, 0.2)
    gphsdv2_dp_4050.SetLineColor(ROOT.kGreen+2)
    gphsdv2_dp_4050.SetLineStyle(1)
    gphsdv2_dp_4050.SetLineWidth(0)
    gphsdv2_dp_4050.SetFillColorAlpha(ROOT.kGreen+2, 0.2)
    gphsdv2_dp_6080.SetLineColor(ROOT.kGreen+2)
    gphsdv2_dp_6080.SetLineStyle(1)
    gphsdv2_dp_6080.SetLineWidth(0)
    gphsdv2_dp_6080.SetFillColorAlpha(ROOT.kGreen+2, 0.2)
    gphsdv2_ds_3040.SetLineColor(ROOT.kAzure+4)
    gphsdv2_ds_3040.SetLineStyle(1)
    gphsdv2_ds_3040.SetLineWidth(0)
    gphsdv2_ds_3040.SetFillColorAlpha(ROOT.kAzure+4, 0.4)
    gphsdv2_ds_4050.SetLineColor(ROOT.kAzure+4)
    gphsdv2_ds_4050.SetLineStyle(1)
    gphsdv2_ds_4050.SetLineWidth(0)
    gphsdv2_ds_4050.SetFillColorAlpha(ROOT.kAzure+4, 0.4)
    gphsdv2_ds_6080.SetLineColor(ROOT.kAzure+4)
    gphsdv2_ds_6080.SetLineStyle(1)
    gphsdv2_ds_6080.SetLineWidth(0)
    gphsdv2_ds_6080.SetFillColorAlpha(ROOT.kAzure+4, 0.4)
    gphsdv2_d0_3040.SetLineColor(ROOT.kRed+1)
    gphsdv2_d0_3040.SetLineStyle(1)
    gphsdv2_d0_3040.SetLineWidth(0)
    gphsdv2_d0_3040.SetFillColorAlpha(ROOT.kRed+1, 0.2)
    gphsdv2_d0_4050.SetLineColor(ROOT.kRed+1)
    gphsdv2_d0_4050.SetLineStyle(1)
    gphsdv2_d0_4050.SetLineWidth(0)
    gphsdv2_d0_4050.SetFillColorAlpha(ROOT.kRed+1, 0.2)
    gphsdv2_d0_6080.SetLineColor(ROOT.kRed+1)
    gphsdv2_d0_6080.SetLineStyle(1)
    gphsdv2_d0_6080.SetLineWidth(0)
    gphsdv2_d0_6080.SetFillColorAlpha(ROOT.kRed+1, 0.2)

    legmodel.AddEntry(gphsdv2_d0_3040, 'D^{0}', 'f')    
    legmodel.AddEntry(gphsdv2_dp_3040, 'D^{+}', 'f')
    legmodel.AddEntry(gphsdv2_ds_3040, 'D_{s}^{+}', 'f')

    cDv2run3, hframes = GetCanvas3sub('cDv2run3', axisname)
    cDv2run3.cd(1).SetLogx()
    yranges = [-0.05, 0.26]
    hframes[0].GetYaxis().SetRangeUser(yranges[0], yranges[1])
    hframes[0].GetXaxis().SetRangeUser(0.48, 26)
    hframes[0].GetXaxis().SetMoreLogLabels()
    line1 = DrawLineAt0(-0., 25)
    line1.Draw()
    latexlarge.DrawLatexNDC(0.22, 0.92, 'ALICE Preliminary')
    latexmedium.DrawLatexNDC(0.22, 0.86, 'Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV, |#it{y}| < 0.8')
    latexdetail2.DrawLatexNDC(0.82, 0.92, '30#minus40%')
    gphsdv2_dp_3040.Draw('3c same')
    gphsdv2_ds_3040.Draw('3c same')
    gphsdv2_d0_3040.Draw('3c same')
    DrawStatSystEmpty(gist_d0_pass4_3040, gist_d0_pass4_syst_3040, False, legempty, markersize=1.2)
    DrawStatSystEmpty(gist_dp_pass4_3040, gist_dp_pass4_syst_3040, False, legempty)
    DrawStatSystEmpty(gist_ds_pass4_3040, gist_ds_pass4_syst_3040, False, legempty)
    cDv2run3.Update()
    
    cDv2run3.cd(2).SetLogx()
    hframes[1].GetXaxis().SetLabelSize(0.055)
    hframes[1].GetXaxis().SetTitleSize(0.055)
    hframes[1].GetXaxis().SetLabelOffset(0.0002)
    hframes[1].GetXaxis().SetTitleOffset(1.08)
    hframes[1].GetYaxis().SetRangeUser(yranges[0], yranges[1])
    hframes[1].GetXaxis().SetRangeUser(0.51, 26)
    hframes[1].GetXaxis().SetMoreLogLabels()
    hframes[1].GetYaxis().SetRangeUser(yranges[0], yranges[1])
    hframes[1].GetXaxis().SetRangeUser(0.51, 26)
    hframes[1].GetXaxis().SetMoreLogLabels()
    latexdetail2.DrawLatexNDC(0.82, 0.92, '40#minus50%')
    line2 = DrawLineAt0(0.1, 25)
    line2.Draw()
    leg.Draw('same')
    legempty.Draw('same')
    legmodel.Draw('same')

    gphsdv2_dp_4050.Draw('3c same')
    gphsdv2_ds_4050.Draw('3c same')
    gphsdv2_d0_4050.Draw('3c same')
    DrawStatSystEmpty(gist_d0_pass4_4050, gist_d0_pass4_syst_4050, False, False, markersize=1.2)
    DrawStatSystEmpty(gist_dp_pass4_4050, gist_dp_pass4_syst_4050, False, False)
    DrawStatSystEmpty(gist_ds_pass4_4050, gist_ds_pass4_syst_4050, False, False)
    
    # Dirty trick to create a white box over the label
    x1 = 0.478
    x2 = 0.665
    y1 = cDv2run3.cd(2).GetUymin()  # Get bottom of the pad
    y2 = cDv2run3.cd(2).GetUymin() + (cDv2run3.cd(2).GetUymax() - cDv2run3.cd(2).GetUymin()) * 0.178  # Small height to cover only label
    box = ROOT.TBox(x1, y1, x2, y2)
    box.SetFillColor(ROOT.kWhite)  # White box
    box.SetLineColor(ROOT.kWhite)  # Hide border
    box.Draw()
    
    cDv2run3.cd(3).SetLogx()
    hframes[2].GetXaxis().SetLabelSize(0.055)
    hframes[2].GetXaxis().SetTitleSize(0.055)
    hframes[2].GetXaxis().SetLabelOffset(0.0002)
    hframes[2].GetXaxis().SetTitleOffset(1.08)
    hframes[2].GetYaxis().SetRangeUser(yranges[0], yranges[1])
    hframes[2].GetXaxis().SetRangeUser(0.51, 26)
    hframes[2].GetXaxis().SetMoreLogLabels()
    hframes[2].GetYaxis().SetRangeUser(yranges[0], yranges[1])
    hframes[2].GetXaxis().SetRangeUser(0.51, 26)
    hframes[2].GetYaxis().SetRangeUser(yranges[0], yranges[1])
    hframes[2].GetXaxis().SetRangeUser(0.22, 12)
    hframes[2].GetXaxis().SetMoreLogLabels()
    latexdetail2.DrawLatexNDC(0.05, 0.92, '60#minus80%')
    line3 = DrawLineAt0(-0.5, 12)
    line3.Draw()
    gphsdv2_dp_6080.Draw('3c same')
    gphsdv2_ds_6080.Draw('3c same')
    gphsdv2_d0_6080.Draw('3c same')
    DrawStatSystEmpty(gist_d0_pass4_6080, gist_d0_pass4_syst_6080, False, False, markersize=1.2)
    DrawStatSystEmpty(gist_dp_pass4_6080, gist_dp_pass4_syst_6080, False, False)
    DrawStatSystEmpty(gist_ds_pass4_6080, gist_ds_pass4_syst_6080, False, False)
    # Dirty trick to create a white box over the label
    x1 = 0.18
    x2 = 0.30
    y1 = -0.075  # Get bottom of the pad
    y2 = -0.055  # Small height to cover only label
    box2 = ROOT.TBox(x1, y1, x2, y2)
    cDv2run3.cd(3)
    box2.SetFillColor(ROOT.kWhite)  # White box
    box2.SetLineColor(ROOT.kWhite)  # Hide border
    box2.Draw()
    
    cDv2run3.Update()
    SaveCanvas(cDv2run3, f'{outdir}Dmeson_v2_run3_vscent_vsmodel', suffix)
    input()

#_____________________________________________________________________________________
# D0 meson v2 in Run 3 pass4 6080 vs models
if plot_d0_v2_run3_6080_vsmodel:
    cRun3D0vsModels, hframe = GetCanvas('cRun3D0vsModels', axisname, ymin=-0.05, ymax=0.36, xmax=12, xmin=0.2)

    gccnuD0 = GetPrediction('/Users/spolitan/cernbox/flowRun3/input/ccnu/Langevin-results-to pxy25.2.13/D-v2-pbpb5.02-60-80-data/D0v2fnaver60-80.dat', 1, 27, 1, False)
    gccnuD02 = GetPrediction('/Users/spolitan/cernbox/flowRun3/input/ccnu/Langevin-results-to pxy25.2.13/D-v2-pbpb5.02-60-80-data/D0v2fnwsnloaver60-80.dat',1, 27, 1, False)
    gccnuD0.SetLineColor(ROOT.kAzure-9)
    gccnuD0.SetLineStyle(1)
    gccnuD0.SetLineWidth(3)
    gccnuD02.SetLineColor(ROOT.kOrange+6)
    gccnuD02.SetLineStyle(1)
    gccnuD02.SetLineWidth(3)

    legmodels = GetLegend(ncolumns=1, header='Langevin #sqrt{#it{s}_{NN}} = 5.02 TeV', xmax=0.5, ymin=0.48, ymax=0.66, textsize=0.045, xmin=0.193)
    legmodels.AddEntry(gccnuD0, 'w/o shadowing', 'l')
    legmodels.AddEntry(gccnuD02, 'w NLO shadowing', 'l')
    gccnuD0.Draw('3c same')
    gccnuD02.Draw('3c same')

    _ = DrawLineAt0(0., 12)
    leg = GetLegend(header='', xmax=0.5, ncolumns=1, ymin=0.66, ymax=0.82, xmin=0.193, textsize=0.045)
    legempty = leg.Clone('')
    legempty.SetHeader(' ')
    leg.AddEntry(gist_d0_pass4_6080, 'Prompt D^{0}', 'p')
    leg.AddEntry(gist_dp_pass4_6080, 'Prompt D^{+}', 'p')
    DrawStatSystEmpty(gist_d0_pass4_6080, gist_d0_pass4_syst_6080, False, legempty, markersize=1)
    DrawStatSystEmpty(gist_dp_pass4_6080, gist_dp_pass4_syst_6080, False, legempty)
    leg.Draw('same')
    legempty.Draw('same')
    legmodels.Draw('same')
    
    latexlarge.DrawLatexNDC(0.20, 0.84, 'ALICE Preliminary')
    latexdetail2.DrawLatexNDC(0.20, 0.78, 'Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV, 60#minus80%')
    latexdetail.DrawLatexNDC(0.8, 0.84, '|#it{y}| < 0.8')
    cRun3D0vsModels.Update()
    SaveCanvas(cRun3D0vsModels, f'{outdir}d0_v2_run3_6080_vsmodel', suffix)
    input()

#_____________________________________________________________________________________
# nonstrange D Run 2 vs Run 3
if plot_d_v2_run2_run3:
    cRun2VsRun3, hframe = GetCanvas('cRun2VsRun3', axisname, ymin=-0.05, ymax=0.46, xmax=42)
    hframe.GetYaxis().SetTitle('#it{v}_{2}{SP}')
    leg = GetLegend(header='#sqrt{#it{s}_{NN}} = 5.36 TeV, |#it{y}| < 0.8, |#Delta#it{#eta}| > 1.3', ymin=0.44,
                    ymax=0.58, ncolumns=1, textsize=0.03, xmax=0.5)
    legempty = leg.Clone('legempty')
    legempty.SetHeader(' ')
    leg.AddEntry(gist_dp_pass4_3040, 'Prompt D^{+}', 'p')
    leg.AddEntry(gist_d0_pass4_3040, 'Prompt D^{0}', 'p')

    legrun2 = GetLegend(header='#sqrt{#it{s}_{NN}} = 5.02 TeV',
                         ncolumns=1, textsize=0.03, ymin=0.58, xmax=0.5, ymax=0.76)
    legempty_run2 = legrun2.Clone('')
    legempty_run2.SetHeader(' ')
    legrun2.AddEntry(gist_run2_av, 'ALICE Prompt D^{0}, D^{+}, D*^{+} average, |#it{y}| < 0.8, |#Delta#it{#eta}| > 0.9', 'p')
    legrun2.AddEntry(gist_run2_av, 'JHEP 01 (2022) 174', '')
    legrun2.AddEntry(gist_d0_cms_run2, 'CMS Prompt D^{0}, |#it{y}| < 1, |#Delta#it{#eta}| > 2.0', 'p')
    legrun2.AddEntry(gist_d0_cms_run2, 'PLB 816 (2021) 136253', '')

    line = DrawLineAt0(0., 41)
    DrawStatSystEmpty(gist_run2_av, gist_run2_av_syst, False, legempty_run2)
    legempty_run2.AddEntry(gist_d0_cms_run2, ' ', '')
    DrawStatSystEmpty(gist_d0_cms_run2, gist_d0_cms_run2_syst, False, legempty_run2)
    legempty_run2.AddEntry(gist_d0_cms_run2, ' ', '')
    DrawStatSystEmpty(gist_dp_pass4_3040, gist_dp_pass4_syst_3040, False, legempty)
    DrawStatSystEmpty(gist_d0_pass4_3050, gist_d0_pass4_syst_3050, False, legempty, markersize=1.2)

    leg.Draw('same')
    legempty.Draw('same')
    legrun2.Draw('same')
    legempty_run2.Draw('same')
    latexlarge.DrawLatexNDC(0.20, 0.84, 'ALICE Preliminary')
    latexdetail2.DrawLatexNDC(0.20, 0.78, 'Pb#minusPb, 30#minus50%')
    cRun2VsRun3.Update()
    SaveCanvas(cRun2VsRun3, f'{outdir}run2_vs_run3_d0_dp', suffix)
    input()

#_____________________________________________________________________________________
# Ds Run 2 vs Run 3 
if plot_ds_v2_run2_run3:
    cRun2VsRun3Ds, hframe = GetCanvas('cRun2VsRun3Ds', axisnamevareta, xmin=0.91, ymin=-0.0, ymax=0.45, xmax=27)

    leg = GetLegend(header='Prompt D_{s}^{+}', ymin=0.53,
                    ymax=0.78, ncolumns=1, textsize=0.035, xmin=0.20, xmax=0.5)
    legempty = leg.Clone('')
    legempty.SetHeader(' ')
    leg.AddEntry(gist_run2_ds, '#sqrt{#it{s}_{NN}} = 5.02 TeV, |#Delta#it{#eta}| > 0.9', 'p')
    leg.AddEntry(gist_run2_ds, 'PLB 827 (2022) 136986', '')
    leg.AddEntry(gist_ds_pass4_3050, '#sqrt{#it{s}_{NN}} = 5.36 TeV, |#Delta#it{#eta}| > 1.3', 'p')

    DrawStatSystEmpty(gist_run2_ds, gist_run2_ds_syst, False, legempty)
    legempty.AddEntry(gist_run2_ds, ' ', '')
    DrawStatSystEmpty(gist_ds_pass4_3050, gist_ds_pass4_syst_3050, False, legempty)

    leg.Draw('same')
    legempty.Draw('same')
    latexlarge.DrawLatexNDC(0.20, 0.84, 'ALICE Preliminary') 
    latexmedium.DrawLatexNDC(0.20, 0.78, 'Pb#minusPb, 30#minus50%')
    latexdetail.DrawLatexNDC(0.8, 0.84, '|#it{y}| < 0.8')
    cRun2VsRun3Ds.Update()
    SaveCanvas(cRun2VsRun3Ds, f'{outdir}run2_vs_run3_ds', suffix)
    input()

#_____________________________________________________________________________________
# Ds Run 2 vs Run 3 vs models
if plot_ds_v2_vsmodel:
    cRun3DsvsModels, hframe = GetCanvas('cRun3DsvsModels', axisname, ymin=0, ymax=0.40, xmax=26, xmin=0.6)

    gtamuv2 = GetPrediction('/Users/spolitan/alice/DmesonAnalysis/models/tamu/PromptDs_TAMU_v2_5TeV_3050.txt',
                                       1, 10, 1, True)    
    gphsdv2 = GetPrediction('/Users/spolitan/flowD/input/models/phsd/PromptDs_v2_3050_frag_coal_smooth.txt',
                            5, 150, 0.1, False)
    gpowlangdv2 = GetPrediction('/Users/spolitan/cernbox/flowRun3/input/pwlang_Ds.txt',
                            125, 933, 0.01, False)
    gtamuv2.SetLineStyle(1)
    gtamuv2.SetLineWidth(0)
    gtamuv2.SetFillColorAlpha(ROOT.kAzure-4, 0.6)
    gphsdv2.SetLineColor(ROOT.kOrange+6)
    gphsdv2.SetLineWidth(4)
    gpowlangdv2.SetLineColor(ROOT.kAzure)
    gpowlangdv2.SetLineStyle(3)
    gpowlangdv2.SetLineWidth(8)

    legmodels = GetLegend(ncolumns=1, header='Transport models #sqrt{#it{s}_{NN}} = 5.02 TeV', xmax=0.5, ymin=0.46, ymax=0.68, textsize=0.04)
    legmodels.AddEntry(gtamuv2, 'TAMU', 'f')
    legmodels.AddEntry(gphsdv2, 'PHSD', 'l')
    legmodels.AddEntry(gpowlangdv2, 'POWLANG', 'l')

    leg = GetLegend(header='', xmax=0.5, ncolumns=1, ymin=0.70, textsize=0.045)
    legempty = leg.Clone('')
    legempty.SetHeader(' ')
    leg.AddEntry(gist_ds_pass4_3050, 'Prompt D_{s}^{+}', 'p')
    
    gtamuv2.Draw('3c same')
    gphsdv2.Draw('c same')
    gpowlangdv2.Draw('c same')
    SetObjectStyle(gist_ds_pass4_3050, color=ROOT.kBlack)
    SetObjectStyle(gist_ds_pass4_syst_3050, color=ROOT.kBlack, markersize=1)
    DrawStatSystEmpty(gist_ds_pass4_3050, gist_ds_pass4_syst_3050, False, legempty, markersize=1)

    leg.Draw('same')
    legempty.Draw('same')
    legmodels.Draw('same')
    
    latexlarge.DrawLatexNDC(0.20, 0.84, 'ALICE Preliminary')
    latexdetail2.DrawLatexNDC(0.20, 0.78, 'Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV, 30#minus50%')
    latexdetail.DrawLatexNDC(0.8, 0.84, '|#it{y}| < 0.8')
    cRun3DsvsModels.Update()
    SaveCanvas(cRun3DsvsModels, f'{outdir}ds_run2_vs_run3_vs_model', suffix)
    input()

#_____________________________________________________________________________________
# Run 2 vs Run 3 nonstrange vs models
if plot_d_v2_vsmodel:
    cRun2VsRun3, _ = GetCanvas('cRun2VsRun3', axisname, ymin=-0.0, ymax=0.42, xmax=26)
    leg = GetLegend(header='', xmax=0.5, ncolumns=1, ymin=0.71, textsize=0.045)
    legempty = leg.Clone('')
    legempty.SetHeader(' ')
    leg.AddEntry(gist_d0_pass4_3050, 'Prompt D^{0}', 'p')

    gtamuv2 = GetPrediction('/Users/spolitan/alice/DmesonAnalysis/models/tamu/PromptD_TAMU_v2_5TeV_3050.txt',
                            1, 10, 1, True)
    gcataniav2 = GetPrediction('/Users/spolitan/flowD/input/models/catania/v2_D0_502_3050_Catania_band.dat',
                            9, 120, 0.1, True)
    gphsdv2 = GetPrediction('/Users/spolitan/flowD/input/models/phsd/PromptD_v2_3050_frag_coal_smooth.txt',
                            5, 150, 0.1, False)
    gpowlangdv2 = GetPrediction('/Users/spolitan/cernbox/flowRun3/input/pwlang_D0.txt',
                            125, 940, 0.01, False)
    glangdv2 = GetPrediction('/Users/spolitan/cernbox/flowRun3/input/langevin-d4-results-to-pxy_25.3.21/Dv2fnwsnlo30-50.dat',
                            1, 35, 1, False)
    glbtdv2 = GetPrediction('/Users/spolitan/cernbox/flowRun3/input/v2_PbPb5020/v2_PbPb5360/v2_D_30-50.dat',
                            25, 260, 0.1, False)
    gepos4hqv2 = GetPrediction('//Users/spolitan/cernbox/flowRun3/input/epos4hq/v2pt_D0_PbPb5.02TeV_30-50.dat',
                            3, 110, 0.1, False)
    splinev2, _, _, _ = ReadLGR('/Users/spolitan/flowD/input/models/lgr/V2_Prompt_Dmeson_30-50__ColRad_FragCoal.dat')
    glgrv2 = ROOT.TGraphAsymmErrors(1)
    for ipt in range(8, 200):
        pt = ipt * 0.1
        glgrv2.AddPoint(pt, splinev2['yCent'](pt))
        glgrv2.SetPointError(ipt, 0., 0.,
                             splinev2['yCent'](pt) - splinev2['yMin'](pt),
                             splinev2['yMax'](pt) - splinev2['yCent'](pt))
    glgrv2.RemovePoint(0)
    for i in range(10):
        gcataniav2.RemovePoint(0)

    gtamuv2.SetLineColor(ROOT.kAzure-4)
    gtamuv2.SetLineStyle(1)
    gtamuv2.SetLineWidth(0)
    gtamuv2.SetFillColorAlpha(ROOT.kAzure-4, 0.6)
    gcataniav2.SetLineWidth(0)
    gcataniav2.SetFillColorAlpha(ROOT.kOrange+10, 0.6)
    gphsdv2.SetLineColor(ROOT.kGreen+3)
    gphsdv2.SetLineStyle(1)
    gphsdv2.SetLineWidth(3)
    glgrv2.SetLineColor(ROOT.kGreen+2)
    glgrv2.SetLineWidth(0)
    glgrv2.SetFillColorAlpha(ROOT.kGreen+2, 0.2)
    gpowlangdv2.SetLineColor(ROOT.kAzure)
    gpowlangdv2.SetLineStyle(3)
    gpowlangdv2.SetLineWidth(8)
    glangdv2.SetLineColor(ROOT.kOrange+1)
    glangdv2.SetLineStyle(7)
    glangdv2.SetLineWidth(3)
    glbtdv2.SetLineColor(ROOT.kCyan+1)
    glbtdv2.SetLineStyle(2)
    glbtdv2.SetLineWidth(3)
    gepos4hqv2.SetLineColor(ROOT.kViolet)
    gepos4hqv2.SetLineStyle(5)
    gepos4hqv2.SetLineWidth(3)

    legmodels = GetLegend(ncolumns=3, header='Transport models #sqrt{#it{s}_{NN}} = 5.02 TeV', xmax=0.98, ymin=0.49, ymax=0.7, textsize=0.038)
    legmodels.AddEntry(gtamuv2, 'TAMU', 'f')
    legmodels.AddEntry(glgrv2, 'LGR', 'f')
    legmodels.AddEntry(gphsdv2, 'PHSD', 'l')
    legmodels.AddEntry(gcataniav2, 'Catania', 'f')
    legmodels.AddEntry(gpowlangdv2, 'POWLANG', 'l')
    legmodels.AddEntry(glangdv2, 'Langevin', 'l')
    legmodels.AddEntry(glbtdv2, 'LBT+PNP', 'l')
    legmodels.AddEntry(gepos4hqv2, 'EPOS4HQ', 'l')

    glgrv2.Draw('3c same')
    gtamuv2.Draw('3c same')
    gphsdv2.Draw('c same')
    gcataniav2.Draw('3c same')
    glangdv2.Draw('c same')
    glbtdv2.Draw('same')
    gpowlangdv2.Draw('c same')
    gepos4hqv2.Draw('same')
    SetObjectStyle(gist_d0_pass4_3050, color=ROOT.kBlack)
    SetObjectStyle(gist_d0_pass4_syst_3050, color=ROOT.kBlack, markersize=1)
    DrawStatSystEmpty(gist_d0_pass4_3050, gist_d0_pass4_syst_3050, False, legempty, markersize=1.)

    leg.Draw('same')
    legempty.Draw('same')
    legmodels.Draw('same')
    latexlarge.DrawLatexNDC(0.20, 0.84, 'ALICE Preliminary')
    latexdetail2.DrawLatexNDC(0.20, 0.78, 'Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV, 30#minus50%')
    latexdetail.DrawLatexNDC(0.8, 0.84, '|#it{y}| < 0.8')
    cRun2VsRun3.Update()
    SaveCanvas(cRun2VsRun3, f'{outdir}run2_vs_run3_vs_model', suffix)
    input()

#_____________________________________________________________________________________
# c, b, h
if plot_v2_compilation:
    cDv2run3, hframe = GetCanvas('cDv2run3', axisnamevareta, ymax=0.5, ymin=-0.2)
    line = DrawLineAt0(0., 40)
    gist_ds_pass4_3050.SetLineColor(ROOT.kRed+2)
    gist_ds_pass4_3050.SetMarkerColor(ROOT.kRed+2)
    gist_ds_pass4_syst_3050.SetLineColor(ROOT.kRed+2)
    gist_ds_pass4_syst_3050.SetMarkerColor(ROOT.kRed+2)
    gist_npd0_pass4_3050.SetLineColor(ROOT.kAzure+4)
    gist_npd0_pass4_3050.SetMarkerColor(ROOT.kAzure+4)
    gist_npd0_pass4_syst_3050.SetLineColor(ROOT.kAzure+4)
    gist_npd0_pass4_syst_3050.SetMarkerColor(ROOT.kAzure+4)
    
    leg = GetLegend(header='Run 3 (30#minus50%)', xmax=0.5, ncolumns=1, ymin=0.54, textsize=0.035, xmin=0.196)
    legempty = leg.Clone('legempty')
    legempty.SetHeader(' ')
    leg.AddEntry(gist_d0_pass4_3050, 'D^{0}', 'p')
    leg.AddEntry(gist_ds_pass4_3050, 'D_{s}^{+}', 'p')
    leg.AddEntry(gist_jpsi_pass4_3050, 'J/#psi', 'p')
    leg.AddEntry(gist_npd0_pass4_3050, 'Non-prompt D^{0}', 'p')
    legrun2 = GetLegend(header='Run 2 (30#minus40%)', xmax=0.7, ncolumns=1, ymin=0.68, xmin=0.50, textsize=0.035)
    legrun2empty = legrun2.Clone('legempty')
    legrun2empty.SetHeader(' ')
    legrun2.AddEntry(gh_run2_alice_3040, '#pi^{+} (JHEP 09 (2018) 006)', 'p')

    DrawStatSystEmpty(gh_run2_alice_3040, ghsyst_run2_alice_3040, False, legrun2empty, markersize=1.15)
    DrawStatSystEmpty(gist_d0_pass4_3050, gist_d0_pass4_syst_3050, False, legempty, markersize=1.2)
    DrawStatSystEmpty(gist_ds_pass4_3050, gist_ds_pass4_syst_3050, False, legempty)
    DrawStatSystEmpty(gist_jpsi_pass4_3050, gist_jpsi_pass4_syst_3050, False, legempty, markersize=1)
    DrawStatSystEmpty(gist_npd0_pass4_3050, gist_npd0_pass4_syst_3050, False, legempty, markersize=1)

    legetarun2 = GetLegend(ncolumns=1, textsize=0.032,
                        xmin=0.165, xmax=0.335,
                        ymin=0.16, ymax=0.30)
    legetarun2.AddEntry(gist_run2_av_syst, 'ALICE Run 2: |#Delta#it{#eta}| > 2.0', ' ')
    legetarun2.AddEntry(gist_run2_av_syst, 'ALICE Run 3: |#Delta#it{#eta}| > 1.3', ' ')

    leg.Draw('same')
    legempty.Draw('same')
    legetarun2.Draw('same')
    legrun2.Draw('same')
    legrun2empty.Draw('same')
    line.Draw('same')
    
    latex.DrawLatexNDC(0.20, 0.86, 'ALICE Preliminary')
    latexdetail.DrawLatexNDC(0.20, 0.8, 'Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV')
    cDv2run3.Update()
    SaveCanvas(cDv2run3, f'{outdir}v2_compilation', suffix)
    input()

#_____________________________________________________________________________________
# D meson vs dueteron
if plot_v2_d_vsdeuteron:
    leg = GetLegend(header='#sqrt{#it{s}_{NN}} = 5.36 TeV, |#it{y}| < 0.8, |#Delta#it{#eta}| > 1.3', xmax=0.4, xmin=0.22, ncolumns=1, ymax=0.68, ymin=0.58, textsize=0.040)
    legempty = leg.Clone('legempty')
    legempty.SetHeader(' ')
    leg.AddEntry(gist_d0_pass4_3040, 'Prompt D^{0}', 'p')
    #leg.AddEntry(gist_dp_pass4_3040, 'D^{+}', 'p')
    #leg.AddEntry(gist_ds_pass4_3040, 'D_{s}^{+}', 'p')
    legrun2 = GetLegend(header='#sqrt{#it{s}_{NN}} = 5.02 TeV, |#it{y}| < 0.5, |#Delta#it{#eta}| > 2.0', xmax=0.4, xmin=0.22, ncolumns=1, ymin=0.7, ymax=0.83, textsize=0.040)
    legrun2empty = legrun2.Clone('legempty')
    legrun2empty.SetHeader(' ')
    legrun2.AddEntry(gh_run2_alice_3040, '#pi^{+} #scale[0.8]{(JHEP 09 (2018) 006)}', 'p')
    legrun2.AddEntry(gist_deu_pass4_3040, 'd #scale[0.8]{(PRC 102 (2020) 055203)}', 'p')
    #legrun2.AddEntry(gist_deu_pass4_3040_276, 'd (2.76 TeV) (EPJC 77 (2017) 658, 2017)', 'p')

    yranges = [-0.0, 0.68]
    cDv2run3, hframes = GetCanvas3sub('cDv2run3', axisnamevareta)
    cDv2run3.cd(1).SetLogx()
    hframes[0].GetYaxis().SetRangeUser(yranges[0], yranges[1])
    hframes[0].GetXaxis().SetRangeUser(0.48, 26)
    hframes[0].GetXaxis().SetMoreLogLabels()

    latexlarge.DrawLatexNDC(0.22, 0.92, 'ALICE Preliminary')
    latexmedium.DrawLatexNDC(0.22, 0.86, 'Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV')
    latexdetail2.DrawLatexNDC(0.82, 0.92, '30#minus40%')
    #ScaleGraph(gist_deu_pass4_3040, 1.25)
    #ScaleGraph(gist_deu_pass4_syst_3040, 1.25)
    #ScaleGraph(gist_d0_pass4_3040, 1/2)
    #ScaleGraph(gist_d0_pass4_syst_3040, 1/2)
    #ScaleGraph(gh_run2_alice_3040, 1/2)
    #ScaleGraph(ghsyst_run2_alice_3040, 1/2)
    DrawStatSystEmpty(gh_run2_alice_3040, ghsyst_run2_alice_3040, False, legrun2empty, markersize=1.2)
    DrawStatSystEmpty(gist_deu_pass4_3040, gist_deu_pass4_syst_3040, False, legrun2empty, markersize=1.)
    #DrawStatSystEmpty(gist_deu_pass4_3040_276, gist_deu_pass4_syst_3040_276, False, legrun2empty, markersize=1.)
    DrawStatSystEmpty(gist_d0_pass4_3040, gist_d0_pass4_syst_3040, False, legempty, markersize=1.2)
    #DrawStatSystEmpty(gist_dp_pass4_3040, gist_dp_pass4_syst_3040, False, legempty)
    #DrawStatSystEmpty(gist_ds_pass4_3040, gist_ds_pass4_syst_3040, False, legempty)
    leg.Draw('same')
    legempty.Draw('same')
    legrun2.Draw('same')
    legrun2empty.Draw('same')
    cDv2run3.Update()
    
    cDv2run3.cd(2).SetLogx()
    hframes[1].GetXaxis().SetLabelSize(0.055)
    hframes[1].GetXaxis().SetTitleSize(0.055)
    hframes[1].GetXaxis().SetLabelOffset(0.0002)
    hframes[1].GetXaxis().SetTitleOffset(1.08)
    hframes[1].GetYaxis().SetRangeUser(yranges[0], yranges[1])
    hframes[1].GetXaxis().SetRangeUser(0.51, 26)
    hframes[1].GetXaxis().SetMoreLogLabels()
    latexdetail2.DrawLatexNDC(0.82, 0.92, '40#minus50%')

    #ScaleGraph(gist_deu_pass4_4050, 1.25)
    #ScaleGraph(gist_deu_pass4_syst_4050, 1.25)
    #ScaleGraph(gist_d0_pass4_4050, 1/2)
    #ScaleGraph(gist_d0_pass4_syst_4050, 1/2)
    #ScaleGraph(gh_run2_alice_4050, 1/2)
    #ScaleGraph(ghsyst_run2_alice_4050, 1/2)
    DrawStatSystEmpty(gh_run2_alice_4050, ghsyst_run2_alice_4050, False, False, markersize=1.2)
    DrawStatSystEmpty(gist_deu_pass4_4050, gist_deu_pass4_syst_4050, False, False, markersize=1.)
    #DrawStatSystEmpty(gist_deu_pass4_4050_276, gist_deu_pass4_syst_4050_276, False, False, markersize=1.2)
    DrawStatSystEmpty(gist_d0_pass4_4050, gist_d0_pass4_syst_4050, False, False, markersize=1.2)
    #DrawStatSystEmpty(gist_dp_pass4_4050, gist_dp_pass4_syst_4050, False, False)
    #DrawStatSystEmpty(gist_ds_pass4_4050, gist_ds_pass4_syst_4050, False, False)
    
    # Dirty trick to create a white box over the label
    x1 = 0.478
    x2 = 0.665
    y1 = -0.055  # Get bottom of the pad
    y2 = -0.009  # Small height to cover only label
    box = ROOT.TBox(x1, y1, x2, y2)
    box.SetFillColor(ROOT.kWhite)  # White box
    box.SetLineColor(ROOT.kWhite)  # Hide border
    box.Draw()
    
    cDv2run3.cd(3).SetLogx()
    hframes[2].GetXaxis().SetLabelSize(0.055)
    hframes[2].GetXaxis().SetTitleSize(0.055)
    hframes[2].GetXaxis().SetLabelOffset(0.0002)
    hframes[2].GetXaxis().SetTitleOffset(1.08)
    hframes[2].GetYaxis().SetRangeUser(yranges[0], yranges[1])
    hframes[2].GetXaxis().SetRangeUser(0.51, 26)
    hframes[2].GetYaxis().SetRangeUser(yranges[0], yranges[1])
    hframes[2].GetXaxis().SetRangeUser(0.48, 17)
    hframes[2].GetXaxis().SetMoreLogLabels()
    latexdetail2.DrawLatexNDC(0.62, 0.92, '#pi^{+}, d: 60#minus70%')
    latexdetail2.DrawLatexNDC(0.62, 0.86, 'D^{0}: 60#minus80%')
    DrawStatSystEmpty(gh_run2_alice_6070, ghsyst_run2_alice_6070, False, False, markersize=1.2)
    DrawStatSystEmpty(gist_deu_pass4_6070, gist_deu_pass4_syst_6070, False, False, markersize=1.)
    DrawStatSystEmpty(gist_d0_pass4_6080, gist_d0_pass4_syst_6080, False, False, markersize=1.2)
    #DrawStatSystEmpty(gist_dp_pass4_6080, gist_dp_pass4_syst_6080, False, False)
    #DrawStatSystEmpty(gist_ds_pass4_6080, gist_ds_pass4_syst_6080, False, False)
    # Dirty trick to create a white box over the label
    cDv2run3.cd(3)
    x1 = 0.478
    x2 = 0.665
    y1 = -0.055  # Get bottom of the pad
    y2 = -0.009  # Small height to cover only label
    box2 = ROOT.TBox(x1, y1, x2, y2)
    box2.SetFillColor(ROOT.kWhite)  # White box
    box2.SetLineColor(ROOT.kWhite)  # Hide border
    box2.Draw()
    SaveCanvas(cDv2run3, f'{outdir}Dmeson_v2_vdeu_run3_3050', suffix)
    input()
    
input('Plotting complete! Press any key to continue.')
