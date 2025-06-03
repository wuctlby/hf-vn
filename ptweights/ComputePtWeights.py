'''
Script for the computation of pT shape weights
run: 
python ./ComputePtweights.py \
        <config> \
        --Bspecie <Bspecie> \
        --suffix <suffix> \
        --fonllD <{fonllD> \
        --fonllB <fonllB> \
        <Raa>
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
from ReadModel import ReadFONLL, ReadTAMU  #pylint: disable=wrong-import-position,import-error
from StyleFormatter import SetObjectStyle     #pylint: disable=wrong-import-position,import-error
from sparse_dicts import get_sparses

def computePtWeights(config):

    # load input configuration
    with open(config, 'r') as ymlCfgFile:
        config = yaml.load(ymlCfgFile, yaml.FullLoader)

    Dspecie = config['Dmeson']

    # 'Ds', 'Dplus', 'Dzero', 'Lc'
    if Dspecie not in ['Ds', 'Dplus', 'Dzero', 'Lc']:
        print(f'ERROR: D specie {Dspecie} not supported! Only Ds, Dplus, Dzero, Lc is supported! Exit')
        sys.exit()
    
    cent = config['Centrality']
    Bspecie = config.get('Bspecie', None)
    rebin = config.get('Rebin', 1) # default rebin
    smooth = config.get('Smooth', 100) # default smooth
    suffix = config.get('Suffix', 'suffix')

    _, _, sparsesGen, axes, _ = get_sparses(config, get_data=False, get_mc=True, debug=False)
    # load thnSparse

    sparseGenD = sparsesGen['GenPrompt'][0]
    sparseGenD.SetName('sparseGenD')
    sparseGenB = sparsesGen['GenFD'][0]
    sparseGenB.SetName('sparseGenB')

    hPtGenD = sparseGenD.Projection(axes['GenPrompt']['Pt'])
    hPtGenD.SetDirectory(0)
    hPtGenD.SetName('hPtGenD')
    hPtGenD.Sumw2()
    hPtGenD.Rebin(rebin)
    hPtGenD.Scale(1./hPtGenD.Integral())

    if Bspecie:
        hPtGenB = sparseGenB.Projection(axes['GenFD']['pt_bmoth'])
        if Dspecie == 'Ds':
            sparseGenBPlusBZero = sparseGenB.Clone('sparseGenBPlusBZero')
            sparseGenBPlusBZero.GetAxis(axes['GenFD']['origin']).SetRange(1, 2)
            sparseGenLambdaBZero = sparseGenB.Clone('sparseGenLambdaBZero')
            sparseGenLambdaBZero.GetAxis(axes['GenFD']['origin']).SetRange(4, 4)
            hPtGenB = sparseGenLambdaBZero.Projection(axes['GenFD']['Pt'])
            hPtGenB.Add(sparseGenBPlusBZero.Projection(axes['GenFD']['Pt']))
        
        #TODO: modifications for other B mesons
        hPtGenB.SetDirectory(0)
        hPtGenB.SetName('hPtGenB')
        hPtGenB.Sumw2()
        hPtGenB.Rebin(rebin)
        hPtGenB.Scale(1./hPtGenB.Integral())
    
    if Bspecie == 'BsBmix':
        sparseGenB.GetAxis(axes['GenFD']['origin']).SetRange(3, 3)
        hPtGenBs = sparseGenB.Projection(axes['GenFD']['pt_bmoth'])
        hPtGenBs.SetDirectory(0)
        hPtGenBs.SetName('hPtGenBs')
        hPtGenBs.Rebin(rebin)
        hPtGenBs.Scale(1./2 * hPtGenBs.Integral())
        hPtGenB.Scale(1./2)
        hPtGenB.Add(hPtGenBs) # assuming 50% Bs and 50% B, reasonable for non-prompt Ds

    print('INFO: MC input loaded')
    # load models predictions
    #___________________________________________________________________________________________________________________________
    sFONLLD, _, ptMinFONLL, ptMaxFONLL = ReadFONLL(config['PtShapeFONLL_D'], True, Dspecie)
    sFONLLB, _, ptMinFONLLB, ptMaxFONLLB = ReadFONLL(config['PtShapeFONLL_B'], True, 'B')
    sTAMU, _, ptMinTAMU, ptMaxTAMU = ReadTAMU(config['RaaTAMU_D'])
    sTAMUB, _, ptMinTAMUB, ptMaxTAMUB = ReadTAMU(config['RaaTAMU_B'])
    if Dspecie == 'Ds' and Bspecie == 'BsBmix' and config.get('RaaTAMU_Bs'):
        sTAMUBs, _, ptMinTAMUBs, ptMaxTAMUBs = ReadTAMU(config['RaaTAMU_Bs'])

    histoDNames = ['hPtFONLLDcent', 'hPtFONLLDmin', 'hPtFONLLDmax']
    histoBNames = ['hPtFONLLBcent', 'hPtFONLLBmin', 'hPtFONLLBmax']
    modelPred = ['yCent', 'yMin', 'yMax']

    hPtFONLLD, hPtFONLLB, hPtFONLLtimesTAMUD, hPtFONLLtimesTAMUB = [], [], [], []
    hPtWeightsFONLLD, hPtWeightsFONLLB, hPtWeightsFONLLtimesTAMUD, hPtWeightsFONLLtimesTAMUB = [], [], [], []
    print('INFO: Start computing pT weights')

    # D meson weights
    #___________________________________________________________________________________________________________________________
    for histoName, pred in zip(histoDNames, modelPred):
        hPtFONLLD.append(hPtGenD.Clone(histoName))
        hPtFONLLtimesTAMUD.append(hPtGenD.Clone(histoName.replace('FONLL', 'FONLLtimesTAMU')))

        for iPt in range(1, hPtFONLLD[-1].GetNbinsX()+1):
            ptCent = hPtFONLLD[-1].GetBinCenter(iPt)
            if ptMinFONLL < ptCent < ptMaxFONLL:
                hPtFONLLD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent))
                if ptMinTAMU <= ptCent <= ptMaxTAMU:
                    hPtFONLLtimesTAMUD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sTAMU['yCent'](ptCent))
                elif ptCent > ptMaxTAMU:
                    hPtFONLLtimesTAMUD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sTAMU['yCent'](ptMaxTAMU))
                else:
                    hPtFONLLtimesTAMUD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sTAMU['yCent'](ptMinTAMU))
            elif ptCent > ptMaxFONLL:
                print(f'WARNING: Results for pT > {ptMaxFONLL} not reliable!')
                continue
            else:
                print(f'WARNING: Results for pT < {ptMinFONLL} not reliable!')
                continue

        hPtFONLLD[-1].Sumw2()
        hPtFONLLD[-1].Scale(1./hPtFONLLD[-1].Integral())
        hPtWeightsFONLLD.append(hPtFONLLD[-1].Clone(histoName.replace('Pt', 'PtWeights')))
        hPtWeightsFONLLD[-1].Divide(hPtFONLLD[-1], hPtGenD)
        hPtWeightsFONLLD[-1].Smooth(smooth)
        hPtFONLLtimesTAMUD[-1].Scale(1./hPtFONLLtimesTAMUD[-1].Integral())
        hPtWeightsFONLLtimesTAMUD.append(
            hPtFONLLtimesTAMUD[-1].Clone(hPtFONLLtimesTAMUD[-1].GetName().replace('Pt', 'PtWeights')))
        hPtWeightsFONLLtimesTAMUD[-1].Divide(hPtFONLLtimesTAMUD[-1], hPtGenD)
        hPtWeightsFONLLtimesTAMUD[-1].Smooth(smooth)
    print('INFO: D pT weights calculated')

    # B meson weights
    #___________________________________________________________________________________________________________________________
    if Bspecie:
        for histoName, pred in zip(histoBNames, modelPred):
            hPtFONLLB.append(hPtGenB.Clone(histoName))
            hPtFONLLtimesTAMUB.append(hPtGenB.Clone(histoName.replace('FONLL', 'FONLLtimesTAMU')))

            for iPt in range(1, hPtFONLLB[-1].GetNbinsX()+1):
                ptCent = hPtFONLLB[-1].GetBinCenter(iPt)
                if ptMinFONLLB < ptCent < ptMaxFONLLB:
                    hPtFONLLB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent))
                    if Bspecie != 'BsBmix':
                        if ptMinTAMUB <= ptCent <= ptMaxTAMUB:
                            hPtFONLLtimesTAMUB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent) * sTAMUB['yCent'](ptCent))
                        elif ptCent > ptMaxTAMUB:
                            hPtFONLLtimesTAMUB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent) * sTAMUB['yCent'](ptMaxTAMUB))
                        else:
                            hPtFONLLtimesTAMUB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent) * sTAMUB['yCent'](ptMinTAMUB))
                    else:
                        ptMaxMix = min([ptMaxTAMUB, ptMaxTAMUBs])
                        ptMinMix = max([ptMinTAMUB, ptMinTAMUBs])
                        if ptMinMix <= ptCent <= ptMaxMix:
                            rAAMix = (sTAMUB['yCent'](ptCent) + sTAMUBs['yCent'](ptCent)) / 2
                        elif ptCent > ptMaxMix:
                            rAAMix = (sTAMUB['yCent'](ptMaxMix) + sTAMUBs['yCent'](ptMaxMix)) / 2
                        else:
                            rAAMix = (sTAMUB['yCent'](ptMinMix) + sTAMUBs['yCent'](ptMinMix)) / 2
                        hPtFONLLtimesTAMUB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent) * rAAMix)
                elif ptCent > ptMaxFONLLB:
                    print(f'WARNING: Results for pT > {ptMaxFONLLB} not reliable! Set weights to 0')
                    continue
                else:
                    print(f'WARNING: Results for pT < {ptMinFONLLB} not reliable! Set weights to 0')
                    continue

            hPtFONLLB[-1].Sumw2()
            hPtFONLLB[-1].Scale(1./hPtFONLLB[-1].Integral())
            hPtWeightsFONLLB.append(hPtFONLLB[-1].Clone(histoName.replace('Pt', 'PtWeights')))
            hPtWeightsFONLLB[-1].Divide(hPtFONLLB[-1], hPtGenB)
            hPtWeightsFONLLB[-1].Smooth(smooth)
            hPtFONLLtimesTAMUB[-1].Scale(1./hPtFONLLtimesTAMUB[-1].Integral())
            hPtWeightsFONLLtimesTAMUB.append(
                hPtFONLLtimesTAMUB[-1].Clone(hPtFONLLtimesTAMUB[-1].GetName().replace('Pt', 'PtWeights')))
            hPtWeightsFONLLtimesTAMUB[-1].Divide(hPtFONLLtimesTAMUB[-1], hPtGenB)
            hPtWeightsFONLLtimesTAMUB[-1].Smooth(smooth)
    print('INFO: B pT weights calculated')

    # save output
    #___________________________________________________________________________________________________________________________
    outputDir = f'{work_dir}/weights/{Dspecie}/{cent}/'
    os.makedirs(outputDir, exist_ok=True)
    outfile = TFile(f'{outputDir}/pTweight_{Dspecie}_{cent}_{suffix}.root', 'recreate')
    hPtGenD.Write()
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
                        f';#it{{p_{{T}}}} (GeV/c);Prompt {config["Dmeson"]}')
    canvPtshape.cd(1)
    canvPtshape.cd(1).SetLogy()

    leg = TLegend(0.5, 0.63, 0.7, 0.83)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.04)

    SetObjectStyle(hPtGenD, color=kRed, markersize=0.5)
    leg.AddEntry(hPtGenD, 'Gen Prompt', 'lp')
    SetObjectStyle(hPtFONLLtimesTAMUD[0], color=kBlack, markersize=0.5)
    leg.AddEntry(hPtFONLLtimesTAMUD[0], 'FONLL #times TAMU (R_{AA})', 'lp')
    SetObjectStyle(hPtFONLLD[0], color=kAzure, markersize=0.5)

    hPtFONLLtimesTAMUD[0].Draw('same')
    hPtGenD.Draw('same')
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
    parser.add_argument("config", type=str, help="flow config file")
    args = parser.parse_args()

    computePtWeights(args.config)