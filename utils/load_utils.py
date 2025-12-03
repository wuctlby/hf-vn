import os
from utils import logger
from ROOT import TFile, TH1
from itertools import combinations

TH1.AddDirectory(False)

def load_root_files(inputPath, prefix: str, suffix='.root') -> list[str]:
    """
    Load root files from a specified directory that match the given prefix and suffix.

    Args:
        inputPath (str): Path to the directory containing the root files.
        prefix (str): Prefix of the files to be loaded.
        suffix (str): Suffix of the files to be loaded, default is '.root'.

    Returns:
        list: List of file paths that match the criteria.
    """
    if os.path.exists(inputPath):
        return sorted(
            [f'{os.path.join(inputPath, file)}'
             for file in os.listdir(inputPath) if file.startswith(prefix) and file.endswith(suffix)]
        )
    else:
        logger(f'No folder found in {inputPath}', level='ERROR')
        raise ValueError(f'No folder found in {inputPath}')

def load_reso_histos(an_res_file, wagon_id):
    '''
    Get list of histograms for SP or EP resolution

    Input:
        - an_res_file:
            str, resolution file
        - wagon_id:
            str, wagon ID

    Output:
        - correct_histo_triplets:
            list of TH2D, list of TH2D objects with the SP product or EP cos(deltaphi) values vs centrality
        - correct_histo_labels:
            list of strings, list of detector labels
    '''
    infile_path = f'hf-task-flow-charm-hadrons'
    if wagon_id:
        infile_path = f'{infile_path}_id{wagon_id}'

    infile_path = f'{infile_path}/spReso'  # for sp method only
    prefix = 'hSpReso'  # for sp method only
    infile = TFile(an_res_file, 'READ')
    directory = infile.GetDirectory(infile_path)
    histos = [key.ReadObj() for key in directory.GetListOfKeys()]
    for histo in histos:
        histo.SetDirectory(0)
    pairs = [key.GetName() for key in directory.GetListOfKeys()]

    # generate triplets of pairs (AB, AC, BC)
    triplets = []
    detsA = ['FT0c', 'FT0a', 'FV0a', 'TPCpos', 'FT0m', 'TPCneg']
    triplets = list(combinations(pairs, 3))
    histo_triplets = list(combinations(histos, 3))
    correct_histo_triplets = []
    correct_histo_labels = []
    for i, triplet in enumerate(triplets):
        for detA in detsA:
            detB = triplet[0].replace(prefix, '').replace(detA, '')
            detC = triplet[1].replace(prefix, '').replace(detA, '')
            if (detA in triplet[0] and detA in triplet[1]) and \
               (detB in triplet[0] and detB in triplet[2]) and \
               (detC in triplet[1] and detC in triplet[2]):
                correct_histo_triplets.append(histo_triplets[i])
                correct_histo_labels.append((detA, detB, detC))

    return correct_histo_triplets, correct_histo_labels

def load_eff_histos(effFiles) -> tuple:
    """
    Load efficiency histograms from a file or a list of files.

    Args:
        effFiles (str or list[str]): Path or list of paths to the efficiency files.

    Returns:
        tuple: If input is a single file, returns histograms for prompt and FD efficiencies, and their fractions and corrected fractions.
               If input is a list, returns lists of those histograms for each file.
    """
    def _load_single_eff_histos(effFile: str):
        f = TFile.Open(effFile)
        hEffPrompt      = f.Get('hEffPrompt')
        hEffFD          = f.Get('hEffFD')
        hPromptFrac     = hEffPrompt.Clone('hPromptFrac')
        hFDFrac         = hEffFD.Clone('hFDFrac')
        hPromptFracCorr = hEffPrompt.Clone('hPromptFracCorr')
        hFDFracCorr     = hEffFD.Clone('hFDFracCorr')
        
        hPromptFrac.SetTitle(';#it{p}_{T} (GeV/#it{c}); #it{f}_{prompt}')
        hFDFrac.SetTitle(';#it{p}_{T} (GeV/#it{c}); #it{f}_{FD}')
        hPromptFracCorr.SetTitle(';#it{p}_{T} (GeV/#it{c}); corrected #it{f}_{prompt}')
        hFDFracCorr.SetTitle(';#it{p}_{T} (GeV/#it{c}); corrected #it{f}_{FD}')
        
        f.Close()
        return hEffPrompt, hEffFD, hPromptFrac, hFDFrac, hPromptFracCorr, hFDFracCorr

    if isinstance(effFiles, str):
        return _load_single_eff_histos(effFiles)
    elif isinstance(effFiles, list):
        hEffPrompts, hEffFDs, hPromptFracs, hFDFracs, hPromptFracCorrs, hFDFracCorrs = [], [], [], [], [], []
        for effFile in effFiles:
            hEffPrompt, hEffFD, hPromptFrac, hFDFrac, hPromptFracCorr, hFDFracCorr = _load_single_eff_histos(effFile)
            hEffPrompts.append(hEffPrompt)
            hEffFDs.append(hEffFD)
            hPromptFracs.append(hPromptFrac)
            hFDFracs.append(hFDFrac)
            hPromptFracCorrs.append(hPromptFracCorr)
            hFDFracCorrs.append(hFDFracCorr)
        return (hEffPrompts, hEffFDs, 
                hPromptFracs, hFDFracs,
                hPromptFracCorrs, hFDFracCorrs)
    else:
        raise TypeError("effFiles must be a str or a list of str")

def load_cutVar_histos(cutVarFracFile: str) -> tuple:
    """
    Load histograms from a cut variation file.

    Args:
        cutVarFracFile (str): Path to the cut variation file.

    Returns:
        tuple: Contains histograms for corrected yields and covariance matrices.
    """
    cutVarFracFile   = TFile.Open(cutVarFracFile)
    hCorrYieldPrompt = cutVarFracFile.Get('hCorrYieldPrompt')
    hCorrYieldFD     = cutVarFracFile.Get('hCorrYieldFD')
    hCovPromptPrompt = cutVarFracFile.Get('hCovPromptPrompt')
    hCovPromptFD     = cutVarFracFile.Get('hCovPromptFD')
    hCovFDFD         = cutVarFracFile.Get('hCovFDFD')
    cutVarFracFile.Close()
    
    return (hCorrYieldPrompt, hCorrYieldFD, 
            hCovPromptPrompt, hCovPromptFD, hCovFDFD)