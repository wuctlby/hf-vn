import ROOT
import json
import os
import numpy as np

file_mass_path = '/home/wuct/ALICE/reps/hf-vn-working/correlations/PostProcessing/InvMass/InvMassVsPt_010_.root'
file_v2Delta_path = '/home/wuct/ALICE/reps/hf-vn-working/correlations/PostProcessing/Output_CorrelationFitting_010_negDeta_Root/CorrPhiD0_FinalPlots.root'
json_path = '/home/wuct/ALICE/reps/hf-vn-working/correlations/PostProcessing/config_CorrAnalysis_v2_010_negDeta.json'
v2_delta_hh=0.7

def get_mass_distribution(MassVsPt, pt_min, pt_max):
    binMin = MassVsPt.GetYaxis().FindBin(pt_min*1.0001)
    binMax = MassVsPt.GetYaxis().FindBin(pt_max*0.9999)
    if pt_min != MassVsPt.GetYaxis().GetBinLowEdge(binMin):
        print(f"Warning: pt_min {pt_min} does not match bin edge {MassVsPt.GetYaxis().GetBinLowEdge(binMin)}")
    if pt_max != MassVsPt.GetYaxis().GetBinUpEdge(binMax):
        print(f"Warning: pt_max {pt_max} does not match bin edge {MassVsPt.GetYaxis().GetBinUpEdge(binMax)}")
    MassVsPt.GetYaxis().SetRange(binMin, binMax)
    proj = MassVsPt.ProjectionX(f"mass_proj_{pt_min}_{pt_max}")
    proj.SetDirectory(0)
    return proj
    

with open(json_path, 'r') as f:
    config = json.load(f)
ptBins = config['binsPtCandIntervals']
ptMins = ptBins[0:len(ptBins)-1]
ptMaxs = ptBins[1:len(ptBins)]
ptAssocBins = config['binsPtHadIntervals']
ptAssocMins = ptAssocBins[0:len(ptAssocBins)-1]
ptAssocMaxs = ptAssocBins[1:len(ptAssocBins)]
massBins = config['binsInvMassIntervals']
massMins = massBins[0:len(massBins)-1]
massMaxs = massBins[1:len(massBins)]

fileMass = ROOT.TFile.Open(file_mass_path)
hMassVsPt = fileMass.Get("hMassVsPt")

fileV2Delta = ROOT.TFile.Open(file_v2Delta_path)
hV2Deltas = {}
for obj in fileV2Delta.GetListOfKeys():
    hV2Deltas[obj.GetName()] = fileV2Delta.Get(obj.GetName())

# sort and slice histograms

for indexPtAssoc, (ptAssocMin, ptAssocMax) in enumerate(zip(ptAssocMins, ptAssocMaxs), start=1):
    mass_distributions = []
    v2_mass_distributions = []
    for indexPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs), start=1):
        tmp_mass = hMassVsPt.Clone(f'tMass_PtCand{ptMin*10}to{ptMax*10}_PtAssoc{ptAssocMin*10}to{ptAssocMax*10}')
        mass_dist = get_mass_distribution(tmp_mass, ptMin, ptMax)
        mass_distributions.append(mass_dist)
        v2_mass_dist = ROOT.TH1F(f'v2Mass_PtCand{ptMin*10}to{ptMax*10}_PtAssoc{ptAssocMin*10}to{ptAssocMax*10}', f'v2Mass_PtCand{ptMin*10}to{ptMax*10}_PtAssoc{ptAssocMin*10}to{ptAssocMax*10}',len(massMins), np.array(massBins, dtype='d'))
        v2_mass_dist.SetTitle("Mass vs V2; Mass (GeV/c^2); V2")
        for indexMass, massMin in enumerate(massMins, start=1):
            hV2Delta_name = f"hv2Delta_PtBinAssoc{indexPtAssoc}_InvMassBin{indexMass}"
            hV2Delta = hV2Deltas.get(hV2Delta_name)
            valueV2Delta = hV2Delta.GetBinContent(indexPt)
            errorYV2Delta = hV2Delta.GetBinError(indexPt)
            v2_mass_dist.SetBinContent(indexMass, valueV2Delta)
            v2_mass_dist.SetBinError(indexMass, errorYV2Delta)
        v2_mass_dist.Scale(1/ROOT.TMath.Sqrt(v2_delta_hh))
        v2_mass_distributions.append(v2_mass_dist)
    # save histograms
    os.makedirs("/home/wuct/ALICE/reps/hf-vn-working/correlations/PostProcessing/Output_MassVsV2", exist_ok=True)
    outFileName = f"/home/wuct/ALICE/reps/hf-vn-working/correlations/PostProcessing/Output_MassVsV2/InvMassVsV2_PtAssoc{ptAssocMin*10}to{ptAssocMax*10  }.root"
    outFile = ROOT.TFile.Open(outFileName, "RECREATE")
    for mass_dist, v2_mass_dist in zip(mass_distributions, v2_mass_distributions):
        mass_dist.Write()
        v2_mass_dist.Write()
    outFile.Close()


