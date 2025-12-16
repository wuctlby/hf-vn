'''
Script for the computation of pT shape weights
run: 
python ./compute_pt_weights.py \
        <cfg>
'''
import os
import sys
import argparse
import yaml
from ROOT import TFile, TCanvas, TLegend  # pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kRed, kAzure # pylint: disable=import-error,no-name-in-module
work_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append('../')
sys.path.append(f"{os.path.dirname(os.path.abspath(__file__))}/../utils")
from utils import logger
from ReadModel import ReadFONLL, ReadTAMU  #pylint: disable=wrong-import-position,import-error
from StyleFormatter import SetObjectStyle     #pylint: disable=wrong-import-position,import-error
sys.path.append(f"{os.path.dirname(os.path.abspath(__file__))}/../src")
from pre_process import get_inputs_sparse, get_input_paths  #pylint: disable=wrong-import-position,import-error

def fill_pt_spectrum(hist, pp_cross_sect, pp_cross_sect_pt_min=None, pp_cross_sect_pt_max=None,
                     RAA=None, RAA_pt_min=None, RAA_pt_max=None):
    """
    Fill histogram h using cross_sect(pt)
    optionally multiplied by RAA(pt) (with boundary saturation).
    """

    for i in range(1, hist.GetNbinsX() + 1):
        pt = hist.GetBinCenter(i)
        if pt < pp_cross_sect_pt_min or pt > pp_cross_sect_pt_max:
            continue  # skip out-of-range bins

        val = pp_cross_sect(pt)
        if RAA is not None:
            if RAA_pt_min <= pt <= RAA_pt_max:
                val *= RAA(pt)
            elif pt > RAA_pt_max:
                val *= RAA(RAA_pt_max)
            else:
                val *= RAA(RAA_pt_min)
        hist.SetBinContent(i, val)

    hist.Sumw2()
    hist.Scale(1.0 / hist.Integral())
    return hist

def compute_pt_weights(cfg):

    # load input cfguration
    with open(cfg, 'r') as ymlCfgFile:
        cfg = yaml.load(ymlCfgFile, yaml.FullLoader)

    cfg_pt_weights = cfg['ptWeights']
    Dspecie = cfg['Dmeson']

    # 'Ds', 'Dplus', 'Dzero', 'Lc'
    if Dspecie not in ['Ds', 'Dplus', 'Dzero', 'Lc']:
        print(f'ERROR: D specie {Dspecie} not supported! Only Ds, Dplus, Dzero, Lc is supported! Exit')
        sys.exit()

    cent = cfg['centrality']
    Bspecie = cfg_pt_weights.get('Bspecie', None)
    rebin = cfg_pt_weights.get('Rebin', 1) # default rebin
    smooth = cfg_pt_weights.get('Smooth', 100) # default smooth
    suffix = cfg_pt_weights.get('Suffix', 'suffix')

    # Retrieve sparse inputs
    #___________________________________________________________________________________________________________________________
    for sparses_cfg in cfg['preprocess']['inputs']:
        if not 'sparses' in sparses_cfg:
            continue
        for sparse_cfg in sparses_cfg['sparses']:
            if sparse_cfg['name'] == 'GenPrompt':
                sparse_gen_prompt = sparse_cfg
                sparse_file_path = get_input_paths(sparses_cfg['files'], "AnalysisResults")[0]
            if sparse_cfg['name'] == 'GenFD':
                sparse_gen_FD = sparse_cfg
                sparse_file_path = get_input_paths(sparses_cfg['files'], "AnalysisResults")[0]

    sparse_file = TFile.Open(sparse_file_path, "read")
    sparseGenPromptD, axesPromptD = get_inputs_sparse(sparse_file, cfg, sparse_gen_prompt)
    sparseGenNonPromptD, axesNonPromptD = get_inputs_sparse(sparse_file, cfg, sparse_gen_FD)
    sparse_file.Close()

    sparseGenPromptD.SetName('sparseGenPromptD')
    sparseGenNonPromptD.SetName('sparseGenNonPromptD')

    hPtGenPromptD = sparseGenPromptD.Projection(axesPromptD['Pt'])
    hPtGenPromptD.SetDirectory(0)
    hPtGenPromptD.SetName('hPtGenPromptD')
    hPtGenPromptD.Sumw2()
    hPtGenPromptD.Rebin(rebin)
    hPtGenPromptD.Scale(1./hPtGenPromptD.Integral())

    if Bspecie:
        if Dspecie == 'Ds' and Bspecie == 'BsBmix':
            # Add LambdaB
            sparseGenDFromLambdaB = sparseGenNonPromptD.Clone('sparseGenDFromLambdaB')
            sparseGenDFromLambdaB.GetAxis(axesNonPromptD['FlagBHad']).SetRange(4, 4)
            hPtGenB = sparseGenDFromLambdaB.Projection(axesNonPromptD['PtBMoth'])

            # Add B+ and B0
            sparseGenDFromBPlusBZero = sparseGenNonPromptD.Clone('sparseGenDFromBPlusBZero')
            sparseGenDFromBPlusBZero.GetAxis(axesNonPromptD['FlagBHad']).SetRange(1, 2)
            hPtGenB.Add(sparseGenDFromBPlusBZero.Projection(axesNonPromptD['PtBMoth']))

            # Add Bs
            sparseGenDFromBs = sparseGenNonPromptD.Clone('sparseGenDFromBs')
            sparseGenDFromBs.GetAxis(axesNonPromptD['FlagBHad']).SetRange(3, 3)
            hPtGenBs = sparseGenDFromBs.Projection(axesNonPromptD['PtBMoth'])
            hPtGenBs.SetDirectory(0)
            hPtGenBs.Scale(1./2 * hPtGenBs.Integral())
            hPtGenB.Scale(1./2)
            hPtGenB.Add(hPtGenBs) # assuming 50% Bs and 50% B, reasonable for non-prompt Ds
        else:
            #TODO: modifications for other B mesons
            hPtGenB = sparseGenNonPromptD.Projection(axesNonPromptD['PtBMoth'])

        hPtGenB.SetName('hPtGenB')
        hPtGenB.SetDirectory(0)
        hPtGenB.Sumw2()
        hPtGenB.Rebin(rebin)
        hPtGenB.Scale(1./hPtGenB.Integral())

    logger('MC input loaded', level='INFO')
    # load models predictions
    #___________________________________________________________________________________________________________________________
    sFONLLD, _, ptMinFONLL, ptMaxFONLL = ReadFONLL(cfg_pt_weights['PtShapeFONLL_D'], True, Dspecie)
    sFONLLB, _, ptMinFONLLB, ptMaxFONLLB = ReadFONLL(cfg_pt_weights['PtShapeFONLL_B'], True, 'B')
    sTAMUD, _, ptMinTAMUD, ptMaxTAMUD = ReadTAMU(cfg_pt_weights['RaaTAMU_D'])
    sTAMUB, _, ptMinTAMUB, ptMaxTAMUB = ReadTAMU(cfg_pt_weights['RaaTAMU_B'])
    if Dspecie == 'Ds' and Bspecie == 'BsBmix' and cfg_pt_weights.get('RaaTAMU_Bs'):
        sTAMUBs, _, ptMinTAMUBs, ptMaxTAMUBs = ReadTAMU(cfg_pt_weights['RaaTAMU_Bs'])

    histoDNames = ['hPtFONLLDcent', 'hPtFONLLDmin', 'hPtFONLLDmax']
    histoBNames = ['hPtFONLLBcent', 'hPtFONLLBmin', 'hPtFONLLBmax']
    modelPred = ['yCent', 'yMin', 'yMax']

    hPtFONLLD, hPtFONLLB, hPtFONLLtimesTAMUD, hPtFONLLtimesTAMUB = [], [], [], []
    hPtWeightsFONLLD, hPtWeightsFONLLB, hPtWeightsFONLLtimesTAMUD, hPtWeightsFONLLtimesTAMUB = [], [], [], []
    logger("Start computing pT weights", level='INFO')

    # D meson weights
    #___________________________________________________________________________________________________________________________
    for histoName, pred in zip(histoDNames, modelPred):
        hPtFONLL = hPtGenPromptD.Clone(histoName)
        hPtFONLLtimesTAMU = hPtGenPromptD.Clone(histoName.replace("FONLL", "FONLLtimesTAMU"))

        fill_pt_spectrum(hPtFONLL, pp_cross_sect = sFONLLD[pred],
                         pp_cross_sect_pt_min = ptMinFONLL, pp_cross_sect_pt_max = ptMaxFONLL)
        hPtFONLLD.append(hPtFONLL)
        fill_pt_spectrum(hPtFONLLtimesTAMU,
                         pp_cross_sect = sFONLLD[pred], pp_cross_sect_pt_min = ptMinFONLL, pp_cross_sect_pt_max = ptMaxFONLL,
                         RAA = sTAMUD["yCent"], RAA_pt_min = ptMinTAMUD, RAA_pt_max = ptMaxTAMUD
        )
        hPtFONLLtimesTAMUD.append(hPtFONLLtimesTAMU)

        hPtWeightsFONLL = hPtFONLL.Clone(histoName.replace("Pt", "PtWeights"))
        hPtWeightsFONLL.Divide(hPtFONLL, hPtGenPromptD)
        hPtWeightsFONLL.Smooth(smooth)
        hPtWeightsFONLLD.append(hPtWeightsFONLL)

        hPtWeightsFONLLtimesTAMU = hPtFONLLtimesTAMU.Clone(hPtFONLLtimesTAMU.GetName().replace("Pt", "PtWeights"))
        hPtWeightsFONLLtimesTAMU.Divide(hPtFONLLtimesTAMU, hPtGenPromptD)
        hPtWeightsFONLLtimesTAMU.Smooth(smooth)
        hPtWeightsFONLLtimesTAMUD.append(hPtWeightsFONLLtimesTAMU)

    # B meson weights
    #___________________________________________________________________________________________________________________________
    if Bspecie:
        for histoName, pred in zip(histoBNames, modelPred):

            hPtFONLL = hPtGenB.Clone(histoName)
            hPtFONLLtimesTAMU = hPtGenB.Clone(histoName.replace("FONLL", "FONLLtimesTAMU"))

            fill_pt_spectrum(hPtFONLL, pp_cross_sect = sFONLLB[pred],
                             pp_cross_sect_pt_min = ptMinFONLL, pp_cross_sect_pt_max = ptMaxFONLL)
            hPtFONLLB.append(hPtFONLL)
            fill_pt_spectrum(hPtFONLLtimesTAMU,
                             pp_cross_sect = sFONLLB[pred], pp_cross_sect_pt_min = ptMinFONLL, pp_cross_sect_pt_max = ptMaxFONLL,
                             RAA = sTAMUB["yCent"] if Bspecie != 'BsBmix' else (lambda pt: (sTAMUB['yCent'](pt) + sTAMUBs['yCent'](pt)) / 2),
                             RAA_pt_min = ptMinTAMUB if Bspecie != 'BsBmix' else min([ptMinTAMUB, ptMinTAMUBs]),
                             RAA_pt_max = ptMaxTAMUB if Bspecie != 'BsBmix' else max([ptMaxTAMUB, ptMaxTAMUBs])
            )
            hPtFONLLtimesTAMUB.append(hPtFONLLtimesTAMU)

            hPtWeightsFONLL = hPtFONLL.Clone(histoName.replace("Pt", "PtWeights"))
            hPtWeightsFONLL.Divide(hPtFONLL, hPtGenB)
            hPtWeightsFONLL.Smooth(smooth)
            hPtWeightsFONLLB.append(hPtWeightsFONLL)

            hPtWeightsFONLLtimesTAMU = hPtFONLLtimesTAMU.Clone(hPtFONLLtimesTAMU.GetName().replace("Pt", "PtWeights"))
            hPtWeightsFONLLtimesTAMU.Divide(hPtFONLLtimesTAMU, hPtGenB)
            hPtWeightsFONLLtimesTAMU.Smooth(smooth)
            hPtWeightsFONLLtimesTAMUB.append(hPtWeightsFONLLtimesTAMU)

    logger("B pT weights calculated", level='INFO')

    # save output
    #___________________________________________________________________________________________________________________________
    outputDir = f'{work_dir}/weights/{Dspecie}/{cent}/'
    os.makedirs(outputDir, exist_ok=True)
    outfile = TFile(f'{outputDir}/pTweight_{Dspecie}_{cent}_{suffix}.root', 'recreate')
    hPtGenPromptD.Write()
    if Bspecie:
        hPtGenB.Write()
    for iHisto, _ in enumerate(hPtFONLLD):
        hPtFONLLD[iHisto].Write()
        hPtWeightsFONLLD[iHisto].Write()
        hPtFONLLtimesTAMUD[iHisto].Write()
        hPtWeightsFONLLtimesTAMUD[iHisto].Write()
        if Bspecie:
            hPtFONLLB[iHisto].Write()
            hPtWeightsFONLLB[iHisto].Write()
            hPtFONLLtimesTAMUB[iHisto].Write()
            hPtWeightsFONLLtimesTAMUB[iHisto].Write()

    # pT shape D
    #___________________________________________________________________________________________________________________________
    canvPtshape = TCanvas('pTshape', 'pTshape', 2000, 900)
    canvPtshape.Divide(2, 1)
    ptD = [0, 36]
    canvPtshape.cd(1).DrawFrame(ptD[0], 0.0000001, ptD[1], 1,
                        f';#it{{p_{{T}}}} (GeV/c);Prompt {cfg["Dmeson"]}')
    canvPtshape.cd(1)
    canvPtshape.cd(1).SetLogy()

    leg = TLegend(0.5, 0.63, 0.7, 0.83)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.04)

    SetObjectStyle(hPtGenPromptD, color=kRed, markersize=0.5)
    leg.AddEntry(hPtGenPromptD, 'Gen Prompt', 'lp')
    SetObjectStyle(hPtFONLLtimesTAMUD[0], color=kBlack, markersize=0.5)
    leg.AddEntry(hPtFONLLtimesTAMUD[0], 'FONLL #times TAMU (R_{AA})', 'lp')
    SetObjectStyle(hPtFONLLD[0], color=kAzure, markersize=0.5)

    hPtFONLLtimesTAMUD[0].Draw('same')
    hPtGenPromptD.Draw('same')
    leg.Draw()

    canvPtshape.cd(2).DrawFrame(0, 0., ptD[1], 5, ';#it{p_{T}} (GeV/c);Ratio')
    canvPtshape.cd(2)

    legR = TLegend(0.5, 0.63, 0.7, 0.83)
    legR.SetFillStyle(0)
    legR.SetBorderSize(0)
    legR.SetTextSize(0.04)

    SetObjectStyle(hPtWeightsFONLLtimesTAMUD[0], color=kBlack, markersize=1)
    legR.AddEntry(hPtWeightsFONLLtimesTAMUD[0], 'FONLL #times TAMU (R_{AA})', 'lp')
    SetObjectStyle(hPtWeightsFONLLD[0], color=kAzure, markersize=1)

    hPtWeightsFONLLtimesTAMUD[0].Draw('same')
    legR.Draw()

    canvPtshape.Write()
    canvPtshape.SaveAs(f'{outputDir}/pTweight_{Dspecie}_{cent}_{suffix}.png')
    canvPtshape.SaveAs(f'{outputDir}/pTweight_{Dspecie}_{cent}_{suffix}.pdf')

    # pT shape B
    #___________________________________________________________________________________________________________________________
    canvPtshapeB = TCanvas('pTshapeB', 'pTshapeB', 2000, 900)
    canvPtshapeB.Divide(2, 1)
    ptB = [0, 72]
    canvPtshapeB.cd(1).DrawFrame(0, 0.0000001, ptB[1], 1, ';#it{p_{T}^{B}} (GeV/c); B hadron')
    canvPtshapeB.cd(1)
    canvPtshapeB.cd(1).SetLogy()

    legB = TLegend(0.5, 0.63, 0.7, 0.83)
    legB.SetFillStyle(0)
    legB.SetBorderSize(0)
    legB.SetTextSize(0.04)

    SetObjectStyle(hPtGenB, color=kRed, markersize=0.5)
    legB.AddEntry(hPtGenB, 'Gen B', 'lp')
    SetObjectStyle(hPtFONLLtimesTAMUB[0], color=kBlack, markersize=0.5)
    legB.AddEntry(hPtFONLLtimesTAMUB[0], 'FONLL #times TAMU (R_{AA})', 'lp')
    SetObjectStyle(hPtFONLLB[0], color=kAzure, markersize=0.5)

    hPtFONLLtimesTAMUB[0].Draw('same')
    hPtGenB.Draw('same')
    legB.Draw()

    canvPtshapeB.cd(2).DrawFrame(0, 0., ptB[1], 10, ';#it{p_{T}}^{B} (GeV/c);Ratio')
    canvPtshapeB.cd(2)

    legBR = TLegend(0.5, 0.63, 0.7, 0.83)
    legBR.SetFillStyle(0)
    legBR.SetBorderSize(0)
    legBR.SetTextSize(0.04)

    SetObjectStyle(hPtWeightsFONLLtimesTAMUB[0], color=kBlack, markersize=1)
    legBR.AddEntry(hPtWeightsFONLLtimesTAMUB[0], 'FONLL #times TAMU (R_{AA})', 'lp')
    SetObjectStyle(hPtWeightsFONLLB[0], color=kAzure, markersize=1)

    hPtWeightsFONLLtimesTAMUB[0].Draw('same')
    legBR.Draw()

    canvPtshapeB.Write()
    canvPtshapeB.SaveAs(f'{outputDir}/pTweightB_{Dspecie}_{cent}_{suffix}.png')
    canvPtshapeB.SaveAs(f'{outputDir}/pTweightB_{Dspecie}_{cent}_{suffix}.pdf')

    outfile.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument("cfg", type=str, help="flow cfg file")
    args = parser.parse_args()

    compute_pt_weights(args.cfg)
