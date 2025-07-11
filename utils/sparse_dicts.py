import sys
from alive_progress import alive_bar
import os
from ROOT import TFile # pyright: ignore # type: ignore
sys.path.append("./")
from utils import logger

def get_sparses_dicts(config):
    
    axes_dict = {}
    ### Data sparse dictionary
    axes_dict['Flow'] = {
        'Mass': 0,
        'Pt': 1,
        'cent': 2,
        'sp': 3,
        'score_bkg': 4,
        'score_FD': 5,
        'occ': 6
    }
    ### MC sparse dictionary
    if config['Dmeson'] == 'Dzero':
        axes_dict['RecoPrompt'] = {
            'score_bkg': 0,
            'score_prompt': 1,
            'score_FD': 2,
            'Mass': 3,
            'Pt': 4,
            'y': 5,
            'cand_type': 6,
            'pt_bmoth': 7,
            'origin': 8,
            'npvcontr': 9,
            'cent': 10,
            'occ': 11,
        }
        axes_dict['RecoFD'] = axes_dict['RecoPrompt']
        axes_dict['RecoRefl'] = axes_dict['RecoPrompt']
        axes_dict['RecoReflPrompt'] = axes_dict['RecoPrompt']
        axes_dict['RecoReflFD'] = axes_dict['RecoPrompt']
        axes_dict['GenPrompt'] = {
            'Pt': 0,
            'pt_bmoth': 1,
            'y': 2,
            'origin': 3,
            'npvcontr': 4,
            'cent': 5,
            'occ': 6
        }
        axes_dict['GenFD'] = {
            'Pt': 0,
            'pt_bmoth': 1,
            'y': 2,
            'origin': 3,
            'npvcontr': 4,
            'cent': 5,
            'occ': 6
        }
    elif config['Dmeson'] == 'Dplus':
        axes_dict['RecoPrompt'] = {
            'Mass': 0,
            'Pt': 1,
            'score_bkg': 2,
            'score_prompt': 3,
            'score_FD': 4,
            'cent': 5,
            'occ': 6,
        }
        axes_dict['RecoFD'] = {
            'Mass': 0,
            'Pt': 1,
            'score_bkg': 2,
            'score_prompt': 3,
            'score_FD': 4,
            'cent': 5,
            'occ': 6,
            'pt_bmoth': 7,
            'flag_bhad': 8,
        }
        axes_dict['GenPrompt'] = {
            'Pt': 0,
            'y': 1,
            'cent': 2,
            'occ': 3
        }
        axes_dict['GenFD'] = {
            'Pt': 0,
            'y': 1,
            'cent': 2,
            'pt_bmoth': 3,
            'flag_bhad': 4,
        }

    elif config['Dmeson'] == 'Ds':
        axes_dict['RecoPrompt'] = {
            'Mass': 0,
            'Pt': 1,
            'cent': 3,  # Check number 2
            'npvcontr': 4,
            'score_bkg': 5,
            'score_prompt': 6,
            'score_FD': 7,
            'occ': 8,
        }
        axes_dict['RecoFD'] = {
            'Mass': 0,
            'Pt': 1,
            'cent': 2,
            'score_bkg': 3,
            'score_prompt': 4,
            'score_FD': 5,
            'pt_bmoth': 6,
            'flag_bhad': 7,
            'occ': 8
        }
        axes_dict['GenPrompt'] = {
            'Pt': 0,
            'y': 1,
            'npvcontr': 2,
            'cent': 3,
            'occ': 4
        }
        axes_dict['GenFD'] = {
            'Pt': 0,
            'y': 1,
            'cent': 2,
            'pt_bmoth': 3,
            'flag_bhad': 4,
            'occ': 5
        }

    return axes_dict

def get_pt_preprocessed_sparses(config, iPt):
    
    logger("Loading preprocessed sparses", level='INFO')
    sparsesFlow, sparsesReco, sparsesGen, axes_dict, resolutions = {}, {}, {}, {}, {}
    pre_cfg = config['preprocess']
    axes_dict['Flow'] = {ax: iax for iax, ax in enumerate(pre_cfg['axes_data']['axis_names'])}
    ptmin = config["ptbins"][iPt]
    ptmax = config["ptbins"][iPt+1]

    if config.get("outdirPrep") and config["outdirPrep"] != "":
        infileprep = TFile(f"{config['outdirPrep']}/preprocess/AnalysisResults_pt_{int(ptmin*10)}_{int(ptmax*10)}.root")
    else:
        infileprep = TFile(f"{config['outdir']}/preprocess/AnalysisResults_pt_{int(ptmin*10)}_{int(ptmax*10)}.root")

    if config["operations"].get("proj_data"):
        for key, _ in pre_cfg["data"].items():
            sparsesFlow[f'Flow_{key}'] = infileprep.Get(f'Data_Flow_{key}/hSparseFlowCharm')
            resolutions[f'Reso_Flow_{key}'] = infileprep.Get(f'Data_Flow_{key}/hResolution')

    if config["operations"].get("proj_mc"):
        subdir = infileprep.Get("MC/Reco")
        for key in subdir.GetListOfKeys():
            obj = key.ReadObj()
            sparsesReco[key.GetName()[1:]] = obj
            axes_dict[key.GetName()[1:]] = {ax: iax for iax, ax in enumerate(pre_cfg['axes_reco']['axis_names'])}

        subdir = infileprep.Get("MC/Gen")
        for key in subdir.GetListOfKeys():
            obj = key.ReadObj()
            sparsesGen[key.GetName()[1:]] = obj
            axes_dict[key.GetName()[1:]] = {ax: iax for iax, ax in enumerate(pre_cfg['axes_gen']['axis_names'])}

    infileprep.Close()

    return sparsesFlow, sparsesReco, sparsesGen, axes_dict, resolutions

def get_sparses(config, get_data=True, get_mc=True, debug=False):
    """Load the sparses and axes infos

    Args:
        config (dict): the flow config dictionary
        get_data (bool, optional): load data sparses. Defaults to True.
        get_mc (bool, optional): load mc sparses. Defaults to True.
        debug (bool, optional): print debug info. Defaults to False.

    Outputs:
        sparsesFlow: thnSparse in the flow task
        sparsesReco: thnSparse of reco level from the D meson task
        sparsesGen: thnSparse of gen level from the D meson task
        axes_dict (dict): dictionary of the axes for each sparse
    """

    sparsesFlow, sparsesReco, sparsesGen, resolutions = {}, {}, {}, {}
    axes_dict = get_sparses_dicts(config)
    pre_cfg = config['preprocess'] if config.get('preprocess') else config
    if get_data:
        logger(f"\t\t[Data] Loading data sparses from: {pre_cfg['data']}")
        infileflow = []
        resolutions = {}
        for name, dataset in pre_cfg['data'].items():
            # Collect all files starting with AnalysisResults_ and ending with .root in the dataset['files'] string
            if isinstance(dataset["files"], str) and not dataset["files"].endswith(".root"):
                list_of_files = [f for f in os.listdir(dataset["files"]) if f.endswith(".root")]
                infileflow = [TFile(os.path.join(dataset["files"], file)) for file in list_of_files]
            elif isinstance(dataset["files"], list):
                if len(dataset["files"]) == 1:
                    if dataset["files"][0].endswith(".root"):
                        infileflow = [TFile(dataset["files"][0])]
                    else:
                        list_of_files = [f for f in os.listdir(dataset["files"][0]) if f.endswith(".root")]
                        infileflow = [TFile(os.path.join(dataset["files"][0], file)) for file in list_of_files]
                elif len(dataset["files"]) > 1:
                    if all(file.endswith(".root") for file in dataset["files"]):
                        infileflow = [TFile(file) for file in dataset["files"]]
                    else:
                        logger("The dataset contains multiple files, but not all of them are root files. Provide a single root file or a list of root files or a directory containing all root files.", level='ERROR')
            else:
                infileflow = [TFile(dataset["files"])] if isinstance(dataset["files"], str) else [TFile(file) for file in dataset["files"]]
            sparsesFlow[f'Flow_{name}'] = []
            with alive_bar(len(infileflow), title=f"[INFO]\t\t[Data] Loading data sparses for {name}") as bar:
                for infile in infileflow:
                    sparsesFlow[f'Flow_{name}'].append(infile.Get('hf-task-flow-charm-hadrons/hSparseFlowCharm'))
                    bar()
            [infile.Close() for infile in infileflow]
            resofile = TFile.Open(dataset["resolution"], 'r')
            det_A = config.get('detA', 'FT0c')
            det_B = config.get('detB', 'FV0a')
            det_C = config.get('detC', 'TPCtot')
            resolutions[f"Reso_Flow_{name}"] = resofile.Get(f'{det_A}_{det_B}_{det_C}/histo_reso_delta_cent')
            resolutions[f"Reso_Flow_{name}"].SetDirectory(0)
            resofile.Close()

    if get_mc:
        print(f"Loading mc sparse from: {pre_cfg['mc']}")
        infiletask = [TFile(pre_cfg['mc'])] if isinstance(pre_cfg['mc'], str) else [TFile(pre_cfg['mc']) for file in pre_cfg['mc']]

        if config['Dmeson'] == 'Dzero':
            sparseD0Path = 'hf-task-d0/hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type'
            sparsesReco['RecoPrompt'] = [file.Get(sparseD0Path) for file in infiletask]
            for ifile in range(len(sparsesReco['RecoPrompt'])):
                sparsesReco['RecoPrompt'][ifile].GetAxis(axes_dict['RecoPrompt']['origin']).SetRange(2, 2)    # make sure it is prompt
                sparsesReco['RecoPrompt'][ifile].GetAxis(axes_dict['RecoPrompt']['cand_type']).SetRange(1, 2) # make sure it is signal

            sparsesReco['RecoFD'] = [file.Get(sparseD0Path) for file in infiletask]
            for ifile in range(len(sparsesReco['RecoFD'])):
                sparsesReco['RecoFD'][ifile].GetAxis(axes_dict['RecoPrompt']['origin']).SetRange(3, 3)       # make sure it is non-prompt
                sparsesReco['RecoFD'][ifile].GetAxis(axes_dict['RecoPrompt']['cand_type']).SetRange(1, 2)    # make sure it is signal

            sparsesReco['RecoRefl'] = [file.Get(sparseD0Path) for file in infiletask]
            for ifile in range(len(sparsesReco['RecoRefl'])):
                sparsesReco['RecoRefl'][ifile].GetAxis(axes_dict['RecoPrompt']['cand_type']).SetRange(3, 4)  # make sure it is reflection

            sparsesReco['RecoReflPrompt'] = [file.Get(sparseD0Path) for file in infiletask]
            for ifile in range(len(sparsesReco['RecoReflPrompt'])):
                sparsesReco['RecoReflPrompt'][ifile].GetAxis(axes_dict['RecoPrompt']['cand_type']).SetRange(3, 4)  # make sure it is reflection
                sparsesReco['RecoReflPrompt'][ifile].GetAxis(axes_dict['RecoPrompt']['origin']).SetRange(2, 2)       # make sure it is prompt   

            sparsesReco['RecoReflFD'] = [file.Get(sparseD0Path) for file in infiletask]
            for ifile in range(len(sparsesReco['RecoReflFD'])):
                sparsesReco['RecoReflFD'][ifile].GetAxis(axes_dict['RecoPrompt']['cand_type']).SetRange(3, 4)    # make sure it is reflection
                sparsesReco['RecoReflFD'][ifile].GetAxis(axes_dict['RecoPrompt']['origin']).SetRange(3, 3)       # make sure it is FD
            #TODO: safety checks for Dmeson reflecton and secondary peak

            sparsesGen['GenPrompt'] = [file.Get('hf-task-d0/hSparseAcc') for file in infiletask]
            for ifile in range(len(sparsesGen['GenPrompt'])):
                sparsesGen['GenPrompt'][ifile].GetAxis(axes_dict['GenPrompt']['origin']).SetRange(2, 2)  # make sure it is prompt

            sparsesGen['GenFD'] = [file.Get('hf-task-d0/hSparseAcc') for file in infiletask]
            for ifile in range(len(sparsesGen['GenFD'])):
                sparsesGen['GenFD'][ifile].GetAxis(axes_dict['GenFD']['origin']).SetRange(3, 3)  # make sure it is non-prompt
        elif config['Dmeson'] == 'Dplus':
            sparsesReco['RecoFD']     = [file.Get('hf-task-dplus/hSparseMassFD') for file in infiletask]
            sparsesReco['RecoPrompt'] = [file.Get('hf-task-dplus/hSparseMassPrompt') for file in infiletask]
            sparsesGen['GenPrompt']   = [file.Get('hf-task-dplus/hSparseMassGenPrompt') for file in infiletask]
            sparsesGen['GenFD']       = [file.Get('hf-task-dplus/hSparseMassGenFD') for file in infiletask]

        elif config['Dmeson'] == 'Ds':
            sparsesReco['RecoPrompt'] = [file.Get('hf-task-ds/MC/Ds/Prompt/hSparseMass') for file in infiletask]
            sparsesReco['RecoFD']     = [file.Get('hf-task-ds/MC/Ds/NonPrompt/hSparseMass') for file in infiletask]
            sparsesGen['GenPrompt']   = [file.Get('hf-task-ds/MC/Ds/Prompt/hSparseGen') for file in infiletask]
            sparsesGen['GenFD']       = [file.Get('hf-task-ds/MC/Ds/NonPrompt/hSparseGen') for file in infiletask]

        [infile.Close() for infile in infiletask]

    logger("Sparses loaded", level='INFO')
    if debug:
        print('\n')
        print('###############################################################')
        for key, value in axes_dict.items():
            logger(f"{key}:", level='DEBUG')
            for sub_key, sub_value in value.items():
                logger(f"    {sub_key}: {sub_value}", level='DEBUG')
        print('###############################################################')
        print('\n')
    return sparsesFlow, sparsesReco, sparsesGen, axes_dict, resolutions
