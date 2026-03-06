import os
import yaml
from ROOT import TFile, TLatex, TColor, TGraph

from plot_utils import (
    LoadGraphAndSyst, GetCanvas, GetLegend, DrawLineAt0, SaveCanvas,
    DrawStatSystEmpty, SetGlobalStyle, GetInvMassHistAndFit, GetCanvas2sub,
    get_langevin_prediction, get_root_constant, get_katz_prediction
)

SetGlobalStyle(
    padleftmargin=0.16, padbottommargin=0.14, padtopmargin=0.08,
    opttitle=1, titleoffsety=1.6, labelsize=0.05, titlesize=0.05,
    labeloffset=0.01, titleoffset=1.2, labelfont=42, titlefont=42
)

#______________________________________________________________________________________
# Load configuration
with open('/home/spolitan/alice/hf-vn/figures/config_sqm26_preliminary.yml', 'r') as f:
    config = yaml.safe_load(f)

#______________________________________________________________________________________
# Settings
plot_invmass_fit = config.get('plot_invmass_fit', False)
plot_charm = config.get('plot_charm', False)
plot_lf = config.get('plot_lf', False)
outdir = config.get('outdir', './preliminary_sqm26')
os.makedirs(outdir, exist_ok=True)
suffix = config.get('suffix', 'oo_barlow')


# ______________________________________________________________________________________
# Load data
gist_dp, gist_dp_syst = LoadGraphAndSyst(
    config['input']['Dplus']['v2file'],
    config['input']['Dplus']['graph'],
    config['input']['Dplus']['syst'],
    config['style']['Dplus']['color'],
    config['style']['Dplus']['marker'],
    markersize=config['style']['Dplus']['markersize'],
    alpha=config['style']['Dplus']['alpha']
)
# D0
gist_d0, gist_d0_syst = LoadGraphAndSyst(
    config['input']['Dzero']['v2file'],
    config['input']['Dzero']['graph'],
    config['input']['Dzero']['syst'],
    config['style']['Dzero']['color'],
    config['style']['Dzero']['marker'],
    markersize=config['style']['Dzero']['markersize'],
    alpha=config['style']['Dzero']['alpha']
)
# Ds
gist_ds, gist_ds_syst = LoadGraphAndSyst(
    config['input']['Ds']['v2file'],
    config['input']['Ds']['graph'],
    config['input']['Ds']['syst'],
    config['style']['Ds']['color'],
    config['style']['Ds']['marker'],
    markersize=config['style']['Ds']['markersize'],
    alpha=config['style']['Ds']['alpha']
)

latex = TLatex()
latex.SetTextFont(42)

#_____________________________________________________________________________________
# Simulation fit
if plot_invmass_fit:
    hInvMassD0, fMassTotD0, fMassBkgD0, hV2VsMassD0, fV2TotD0, fV2BkgD0, fMassReflD0 = \
        GetInvMassHistAndFit(config['input']['Dzero']['fitfile'],
                             config['simfit_style']['Dzero']['ptmin'],
                             config['simfit_style']['Dzero']['ptmax'],
                             config['simfit_style']['Dzero']['ptindex'],
                             config['simfit_style']['Dzero']['hasfeflections'])
    hInvMassDp, fMassTotDp, fMassBkgDp, hV2VsMassDp, fV2TotDp, fV2BkgDp = \
        GetInvMassHistAndFit(config['input']['Dplus']['fitfile'],
                             config['simfit_style']['Dplus']['ptmin'],
                             config['simfit_style']['Dplus']['ptmax'],
                             config['simfit_style']['Dplus']['ptindex'], False)
    hInvMassDs, fMassTotDs, fMassBkgDs, hV2VsMassDs, fV2TotDs, fV2BkgDs = \
        GetInvMassHistAndFit(config['input']['Ds']['fitfile'],
                             config['simfit_style']['Ds']['ptmin'],
                             config['simfit_style']['Ds']['ptmax'],
                             config['simfit_style']['Ds']['ptindex'], False)

    xmins_mass = [config['simfit_style']['Dzero']['massmin'],
                  config['simfit_style']['Dplus']['massmin'], config['simfit_style']['Ds']['massmin']]
    xmaxs_mass = [config['simfit_style']['Dzero']['massmax'],
                  config['simfit_style']['Dplus']['massmax'], config['simfit_style']['Ds']['massmax']]
    ymins_mass = [config['simfit_style']['Dzero']['ymin_mass'],
                  config['simfit_style']['Dplus']['ymin_mass'], config['simfit_style']['Ds']['ymin_mass']]
    ymaxs_mass = [config['simfit_style']['Dzero']['ymax_mass'],
                  config['simfit_style']['Dplus']['ymax_mass'], config['simfit_style']['Ds']['ymax_mass']]
    ymins_v2 = [config['simfit_style']['Dzero']['ymin_v2'],
                config['simfit_style']['Dplus']['ymin_v2'], config['simfit_style']['Ds']['ymin_v2']]
    ymaxs_v2 = [config['simfit_style']['Dzero']['ymax_v2'],
                config['simfit_style']['Dplus']['ymax_v2'], config['simfit_style']['Ds']['ymax_v2']]
    axisnamebottoms = [config['simfit_style']['Dzero']['axisnamebottom'],
                       config['simfit_style']['Dplus']['axisnamebottom'], config['simfit_style']['Ds']['axisnamebottom']]
    axisnametops = [f';;Counts per {hInvMassD0.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}',
                    f';;Counts per {hInvMassDp.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}',
                    f';;Counts per {hInvMassDs.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}']

    legD = GetLegend(xmax=0.5, ncolumns=1, ymin=0.35,
                     ymax=0.6, textsize=0.05, xmin=0.25)
    legD.AddEntry(hInvMassDs, 'Data', 'p')
    legD.AddEntry(fV2TotDs, 'Total fit function', 'l')
    legD.AddEntry(fV2BkgDs, 'Combinatorial background', 'l')
    legD.AddEntry(fMassReflD0, 'K#minus#pi reflected', 'f')

    cD0v2run3, _ = GetCanvas2sub(
        'cDv2run3',
        xmins_mass[0], xmaxs_mass[0],
        ymins_mass[0], ymaxs_mass[0],
        ymins_v2[0], ymaxs_v2[0],
        axisnametops[0], axisnamebottoms[0]
    )
    cD0v2run3.cd(1)
    hInvMassD0.GetYaxis().SetMaxDigits(3)
    hInvMassD0.Draw('esame')
    fMassReflD0.Draw('same')
    fMassBkgD0.Draw('same')
    fMassTotD0.Draw('same')
    legD.Draw()

    lat_pad1 = config['simfit_style']['Dzero']['latex']['pad1']
    for text in lat_pad1:
        latex.SetTextSize(text[-1])
        latex.DrawLatexNDC(text[0], text[1], text[2])

    cD0v2run3.cd(2)
    lat_pad2 = config['simfit_style']['Dzero']['latex']['pad2']
    for text in lat_pad2:
        latex.SetTextSize(text[-1])
        latex.DrawLatexNDC(text[0], text[1], text[2])
    hV2VsMassD0.Draw('esame')
    fV2BkgD0.Draw('same')
    fV2TotD0.Draw('same')
    SaveCanvas(cD0v2run3, f'{outdir}D0_invmassfit', suffix)

    cDpv2run3, _ = GetCanvas2sub(
        'cDv2run3',
        xmins_mass[1], xmaxs_mass[1],
        ymins_mass[1], ymaxs_mass[1],
        ymins_v2[1], ymaxs_v2[1],
        axisnametops[1], axisnamebottoms[1]
    )
    legD = GetLegend(xmax=0.5, ncolumns=1, ymin=0.05,
                     ymax=0.24, textsize=0.06, xmin=0.25)
    legD.AddEntry(hInvMassDs, 'Data', 'p')
    legD.AddEntry(fV2TotDs, 'Total fit function', 'l')
    legD.AddEntry(fV2BkgDs, 'Combinatorial background', 'l')
    cDpv2run3.cd(1)
    hInvMassDp.Draw('esame')
    fMassBkgDp.Draw('same')
    fMassTotDp.Draw('same')
    legD.Draw()

    lat_pad1 = config['simfit_style']['Dplus']['latex']['pad1']
    for text in lat_pad1:
        latex.SetTextSize(text[-1])
        latex.DrawLatexNDC(text[0], text[1], text[2])

    cDpv2run3.cd(2)
    lat_pad2 = config['simfit_style']['Dplus']['latex']['pad2']
    for text in lat_pad2:
        latex.SetTextSize(text[-1])
        latex.DrawLatexNDC(text[0], text[1], text[2])
    hV2VsMassDp.GetXaxis().SetNdivisions(510)
    hV2VsMassDp.Draw('esame')
    fV2BkgDp.Draw('same')
    fV2TotDp.Draw('same')
    SaveCanvas(cDpv2run3, f'{outdir}Dp_invmassfit', suffix)

    cDsv2run3, _ = GetCanvas2sub(
        'cDv2run3',
        xmins_mass[2], xmaxs_mass[2],
        ymins_mass[2], ymaxs_mass[2],
        ymins_v2[2], ymaxs_v2[2],
        axisnametops[2], axisnamebottoms[2]
    )
    legD = GetLegend(xmax=0.5, ncolumns=1, ymin=0.4,
                     ymax=0.6, textsize=0.06, xmin=0.25)
    legD.AddEntry(hInvMassDs, 'Data', 'p')
    legD.AddEntry(fV2TotDs, 'Total fit function', 'l')
    legD.AddEntry(fV2BkgDs, 'Combinatorial background', 'l')
    cDsv2run3.cd(1)
    hInvMassDs.Draw('esame')
    fMassBkgDs.Draw('same')
    fMassTotDs.Draw('same')
    legD.Draw()
    lat_pad1 = config['simfit_style']['Ds']['latex']['pad1']
    for text in lat_pad1:
        latex.SetTextSize(text[-1])
        latex.DrawLatexNDC(text[0], text[1], text[2])
    cDsv2run3.cd(2)
    lat_pad2 = config['simfit_style']['Ds']['latex']['pad2']
    for text in lat_pad2:
        latex.SetTextSize(text[-1])
        latex.DrawLatexNDC(text[0], text[1], text[2])
    hV2VsMassDs.Draw('esame')
    fV2BkgDs.Draw('same')
    fV2TotDs.Draw('same')
    SaveCanvas(cDsv2run3, f'{outdir}Ds_invmassfit', suffix)
    input('Simulation fit canvases done. Press Enter to continue.')

#______________________________________________________________________________
# Charm v2 in OO 0-20%
if plot_charm:

    # read the Katz prediction from .txt
    gpred = get_katz_prediction(config['input']['Katz']['v2file'])
    gpred.SetLineWidth(
        get_root_constant(config['style']['Katz']['linewidth']))
    gpred.SetLineColor(get_root_constant(TColor.GetColorTransparent(
        get_root_constant(config['style']['Katz']['color']), config['style']['Katz']['alpha'])))
    gpred.SetLineStyle(
        get_root_constant(config['style']['Katz']['linestyle']))

    # read the Langevin prediction from .dat
    gpred_langevin = get_langevin_prediction(config['input']['Langevin']['v2file'])
    gpred_langevin.SetLineWidth(
        get_root_constant(config['style']['Langevin']['linewidth']))
    gpred_langevin.SetLineColor(TColor.GetColorTransparent(get_root_constant(
        config['style']['Langevin']['color']), config['style']['Langevin']['alpha']))
    gpred_langevin.SetLineStyle(
        get_root_constant(config['style']['Langevin']['linestyle']))

    # read the EPOS4HQ prediction from .dat
    grpred_epos4hq = get_langevin_prediction(config['input']['Epos4HQ']['v2file'])
    grpred_epos4hq.SetLineWidth(
        get_root_constant(config['style']['Epos4HQ']['linewidth']))
    grpred_epos4hq.SetLineColor(TColor.GetColorTransparent(get_root_constant(
        config['style']['Epos4HQ']['color']), config['style']['Epos4HQ']['alpha']))
    grpred_epos4hq.SetLineStyle(
        get_root_constant(config['style']['Epos4HQ']['linestyle']))

    cDv2run3, hframe = GetCanvas(
        'cDv2run3', config['charm_v2_style']['axisname'],
        ymax=config['charm_v2_style']['ymax'], xmin=config['charm_v2_style']['xmin'],
        ymin=config['charm_v2_style']['ymin'], xmax=config['charm_v2_style']['xmax'],
        setlogx=config['charm_v2_style']['setlogx']
    )
    cDv2run3.SetTopMargin(0.05)
    hframe.GetYaxis().SetTitleOffset(1.4)
    hframe.GetYaxis().SetLabelSize(0.04)
    hframe.GetXaxis().SetLabelSize(0.04)
    gpred.GetXaxis().SetRangeUser(
        config['charm_v2_style']['xmin']+0.5, config['charm_v2_style']['xmax']-0.5)

    gpred.Draw('same C')
    gpred_langevin.Draw('same C')
    grpred_epos4hq.Draw('same C')
    DrawStatSystEmpty(gist_ds, gist_ds_syst, False, False)
    DrawStatSystEmpty(gist_dp, gist_dp_syst, False, False)
    DrawStatSystEmpty(gist_d0, gist_d0_syst, False, False)

    leg = GetLegend(header='', xmax=0.6, ncolumns=1, ymin=0.58,
                    textsize=0.033, xmin=0.196, ymax=0.78)
    legempty = GetLegend(header='', xmax=0.9, ncolumns=4, ymin=0.74,
                         textsize=0.033, xmin=0.196, ymax=0.88)
    leg.AddEntry(gist_d0, 'Prompt D^{0}', 'p')
    leg.AddEntry(gist_dp, 'Prompt D^{+}', 'p')
    leg.AddEntry(gist_ds, 'Prompt D_{s}^{+}', 'p')
    legmodels = GetLegend(header='', xmin=0.65, ncolumns=1,
                          ymin=0.58, textsize=0.033, xmax=0.96, ymax=0.78)
    legmodels.AddEntry(gpred, 'Katz et al.', 'l')
    legmodels.AddEntry(gpred_langevin, 'Langevin', 'l')
    legmodels.AddEntry(grpred_epos4hq, 'EPOS4HQ', 'l')
    leg.Draw('same')
    legmodels.Draw('same')
    line = DrawLineAt0(
        config['charm_v2_style']['xmin'], config['charm_v2_style']['xmax']-0.2)

    latexdetail = config['charm_v2_style']['latex']
    for text in latexdetail:
        latex.SetTextSize(text[-1])
        latex.DrawLatexNDC(text[0], text[1], text[2])
    line.Draw('same')
    cDv2run3.Update()
    SaveCanvas(cDv2run3, f'{outdir}charm_v2_3050', suffix)
    input('Charm v2 done. Press Enter to continue.')

#_____________________________________________________________________________________
# Light-flavor vs D0 v2 in OO 0-20%
if plot_lf:
    cDv2run3, hframe = GetCanvas(
        'cDv2run3lf', config['lf_v2_style']['axisname'],
        ymax=config['lf_v2_style']['ymax'], xmin=config['lf_v2_style']['xmin'],
        ymin=config['lf_v2_style']['ymin'], xmax=config['lf_v2_style']['xmax'],
        setlogx=config['lf_v2_style']['setlogx']
    )
    cDv2run3.SetTopMargin(0.05)
    hframe.GetYaxis().SetTitleOffset(1.4)
    hframe.GetYaxis().SetLabelSize(0.04)
    hframe.GetXaxis().SetLabelSize(0.04)

    gist_k0, gist_k0_syst = LoadGraphAndSyst(
        config['input']['K0']['v2file'], config['input']['K0']['graph'],
        config['input']['K0']['syst'], config['style']['K0']['color'],
        config['style']['K0']['marker'], markersize=config['style']['K0']['markersize'],
        alpha=config['style']['K0']['alpha']
    )
    gist_lambda, gist_lambda_syst = LoadGraphAndSyst(
        config['input']['Lambda']['v2file'], config['input']['Lambda']['graph'],
        config['input']['Lambda']['syst'], config['style']['Lambda']['color'],
        config['style']['Lambda']['marker'], markersize=config['style']['Lambda']['markersize'],
        alpha=config['style']['Lambda']['alpha']
    )
    if config['input']['Phi']['v2file']:
        gist_phi, gist_phi_syst = LoadGraphAndSyst(
            config['input']['Phi']['v2file'], config['input']['Phi']['graph'],
            config['input']['Phi']['syst'], config['style']['Phi']['color'],
            config['style']['Phi']['marker'], markersize=config['style']['Phi']['markersize'],
            alpha=config['style']['Phi']['alpha']
        )
    # set error bar on x to match the bin width of D0 points
    lf_bins = [0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0]
    for i in range(gist_k0.GetN()+1):
        delta_x = (lf_bins[i] - lf_bins[i-1]) / 2
        gist_k0.SetPointError(i-1, delta_x, gist_k0.GetErrorY(i-1))
        gist_lambda.SetPointError(i-1, delta_x, gist_lambda.GetErrorY(i-1))
        gist_phi.SetPointError(i-1, delta_x, gist_phi.GetErrorY(i-1)) if config['input']['Phi']['v2file'] else None


    DrawStatSystEmpty(gist_k0, gist_k0_syst, False, False)
    DrawStatSystEmpty(gist_lambda, gist_lambda_syst, False, False)
    DrawStatSystEmpty(gist_phi, gist_phi_syst, False, False) if config['input']['Phi']['v2file'] else None
    DrawStatSystEmpty(gist_d0, gist_d0_syst, False, False)

    leg = GetLegend(header='SP, |#Delta#it{#eta}| > 1.3', xmax=0.6, ncolumns=2,
                    ymin=0.68, textsize=0.033, xmin=0.196, ymax=0.78)
    legempty = GetLegend(header='', xmax=0.9, ncolumns=4, ymin=0.74,
                         textsize=0.033, xmin=0.196, ymax=0.88)
    leg.AddEntry(gist_d0, 'Prompt D^{0}', 'p')
    leglf = GetLegend(header='2PC, |#Delta#it{#eta}| > 1.2', xmax=0.8,
                      ncolumns=1, ymin=0.6, textsize=0.033, xmin=0.7, ymax=0.76)
    leglf.AddEntry(gist_k0, 'K^{0}_{S}', 'p')
    leglf.AddEntry(gist_phi, '#phi', 'p') if config['input']['Phi']['v2file'] else None
    leglf.AddEntry(gist_lambda, '#Lambda', 'p')
    leglf.Draw('same')

    leg.Draw('same')
    leglf.Draw('same')
    line = DrawLineAt0(0., 8.2)
    latexdetail = config['lf_v2_style']['latex']
    for text in latexdetail:
        latex.SetTextSize(text[-1])
        latex.DrawLatexNDC(text[0], text[1], text[2])

    line.Draw('same')
    cDv2run3.Update()
    SaveCanvas(cDv2run3, f'{outdir}charm_v2_3050_lf', suffix)
    input('Light-flavor vs D0 v2 done. Press Enter to continue.')

input('Plotting complete! Press any key to continue.')
