from asyncio import current_task
import sys
import argparse
import yaml
import numpy as np
import pathlib as PATH
from concurrent.futures import ThreadPoolExecutor, as_completed
from ROOT import TFile, TH1D
from alive_progress import alive_bar
import os
os.environ['ZFIT_DISABLE_TF_WARNINGS'] = '1'
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
import tensorflow as tf
tf.config.threading.set_intra_op_parallelism_threads(4)
tf.config.threading.set_inter_op_parallelism_threads(1)
print(f"Intra-op threads: {tf.config.threading.get_intra_op_parallelism_threads()}")
print(f"Inter-op threads: {tf.config.threading.get_inter_op_parallelism_threads()}")
sys.path.append("../flareflyfitter/")
from raw_yield_fitter import RawYieldFitter
# from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import fitz
from PyPDF2 import PdfWriter, PdfReader, Transformation

def build_raw_yield_fitter(task, fitConfig):
    task_id = task['taskID']
    task_name = task['taskName']
    pt_min = task['ptMin']
    pt_max = task['ptMax']
    pt_str = f"pt_{int(pt_min*10)}_{int(pt_max*10)}"
    i_pt = task['iPtCand'] if 'iPtCand' in task else task['iPt']
    mass_fit_range = fitConfig['MassFitRanges'][i_pt]
    minimizer = fitConfig.get('minimizer', 'flarefly')
    
    sgn_func = {
        'func': fitConfig.get("SgnFunc")[i_pt] if isinstance(fitConfig.get("SgnFunc"), list) else fitConfig.get("SgnFunc"),
        'label': 'sgn',
        'particle': task['Dmeson']
        }
    bkg_func = {
        'func': fitConfig.get("BkgFunc")[i_pt] if isinstance(fitConfig.get("BkgFunc"), list) else fitConfig.get("BkgFunc"),
        'label': 'Comb_bkg'
    }

    # TODO: correlated background
    # if fitConfig.get('corr_bkgs')
    #     corr_bkg_funcs = 

    rebin = fitConfig.get("Rebin")[i_pt] if isinstance(fitConfig.get("Rebin"), list) else fitConfig.get("Rebin", 1)
    fit_pars = fitConfig.get('InitPars', None)

    fitter = RawYieldFitter(task['Dmeson'], pt_min, pt_max, pt_str, minimizer)
    fitter.set_fit_range(mass_fit_range[0], mass_fit_range[1])
    fitter.add_sgn_func(sgn_func["func"], sgn_func["label"], sgn_func["particle"])
    fitter.add_bkg_func(bkg_func["func"], bkg_func["label"])
    fitter.set_name(f'{task_id}_{task_name}')
    fitter.set_rebin(rebin)
    if fit_pars:
        fitter.set_fit_pars(fit_pars, pt_min, pt_max)

    return fitter

def combine_pdfs(pdf_paths):
    from PyPDF2 import PdfMerger, PdfReader
    merger = PdfMerger()
    sum_pdf = pdf_paths[0].parent / f"Summary_Fits.pdf"
    sum_pdf_combined = pdf_paths[0].parent / f"Summary_Fits_Combined.pdf"
    for pdf_path in pdf_paths:
        merger.append(str(pdf_path))
    merger.write(str(sum_pdf))
    merger.close()
    reader = PdfReader(str(sum_pdf))
    cols = 5
    rows = 4
    margins = 0.5
    width = reader.pages[0].mediabox.width * cols
    height = reader.pages[0].mediabox.height * rows
    writer = PdfWriter()
    for i in range(0, len(reader.pages)):
        page = reader.pages[i]
        row = i // cols
        col = i % cols
        x_offset = col * reader.pages[0].mediabox.width
        y_offset = (rows - 1 - row) * reader.pages[0].mediabox.height
        transformation = Transformation().translate(x_offset, y_offset)
        page.add_transformation(transformation)
        writer.add_page(page)
        
    with open(sum_pdf_combined, "wb") as f_out:
        writer.write(f_out)
    docs = fitz.open(sum_pdf_combined)
    for doc in docs:
        doc.select([0])
        page = doc[0]
        mat = fitz.Matrix(2, 2) 
        pix = page.get_pixmap(matrix=mat, alpha=False)
        pix.save(str(sum_pdf.parent / f"Summary_Fits_Combined.png"))
    docs.close()
    # size = reader.pages[0].mediabox
    # width = size.width
    # height = size.height


def interface_raw_yield_fitter(config_path):
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    method = config.get("method", "")
    if method != "DeltaPhiBinning":
        raise ValueError(f"Unsupported method: {method}. Only 'DeltaPhiBinning' is supported in this interface.")

    Dmeson = config["Dmeson"]
    ptBinsCand = config["ptBinsCand"]
    ptBinsHad = config["ptBinsHad"]
    nPtBinsCand = len(ptBinsCand) - 1
    deltaPhiBins = config.get("deltaPhiBins", list(np.linspace(-1.571, 4.712, 24)))  # default 24 bins from -pi/2 to 3pi/2

    outdir = PATH.Path(config["outdir"])
    suffix = config.get("suffix", "")
    outFilePath = outdir / f"AssociatedPairsYields"
    outFilePath.mkdir(parents=True, exist_ok=True)

    inFile = TFile.Open(str(outdir / f"CorrelExtract_{suffix}" / "CorrelationsResults" / "CorrelationsResults.root"))

    # build tasks
    tasks = []
    datas = {}
    hPairsYields_vs_DeltaPhi = {}
    with alive_bar(nPtBinsCand * len(ptBinsHad[:-1]) * len(deltaPhiBins[:-1]), title="Building tasks and loading histograms") as bar:
        for iPtCand, (ptMin, ptMax) in enumerate(zip(ptBinsCand[:-1], ptBinsCand[1:])):
            for iPtHad, (ptHadMin, ptHadMax) in enumerate(zip(ptBinsHad[:-1], ptBinsHad[1:])):
                hTemp = TH1D(f"hPairsYields_vs_DeltaPhi_{iPtCand}_{iPtHad}", "Associated Pair Raw Yields vs #Delta#phi;#Delta#phi (rad); Associated pairs raw yield", len(deltaPhiBins)-1, np.array(deltaPhiBins, dtype='double'))
                hPairsYields_vs_DeltaPhi[(iPtCand, iPtHad)] = hTemp
                for iDeltaPhi, (deltaPhiMin, deltaPhiMax) in enumerate(zip(deltaPhiBins[:-1], deltaPhiBins[1:])):
                    # build task
                    task_id = len(tasks)
                    task = {
                        "taskID": task_id, 
                        "taskName": f"CandPt_{int(ptMin*10)}_{int(ptMax*10)}_HadPt_{int(ptHadMin*10)}_{int(ptHadMax*10)}_dPhi_{int(deltaPhiMin*1000)}_{int(deltaPhiMax*1000)}",
                        "iPtCand": iPtCand, "ptMin": ptMin, "ptMax": ptMax,
                        "iPtHad": iPtHad, "ptHadMin": ptHadMin, "ptHadMax": ptHadMax,
                        "iDeltaPhi": iDeltaPhi, "deltaPhiMin": deltaPhiMin, "deltaPhiMax": deltaPhiMax,
                        "Dmeson": Dmeson,
                        "pdfPath": outFilePath / f"PtCandBin_{int(ptMin*10)}_{int(ptMax*10)}" / f"PtHadBin_{int(ptHadMin*10)}_{int(ptHadMax*10)}" / f"TempFitResult_DeltaPhiBin_{int(deltaPhiMin*1000)}_{int(deltaPhiMax*1000)}.pdf",
                    }
                    tasks.append(task)
                    # load histogram for the task
                    histo_name = f"PtCandBin_{int(ptMin*10)}_{int(ptMax*10)}/PtHadBin_{int(ptHadMin*10)}_{int(ptHadMax*10)}/DeltaPhiBin_{int(deltaPhiMin*1000)}_{int(deltaPhiMax*1000)}/hCorrectedMassPairs_AllPools"
                    bar.text = f"Loading histogram {histo_name} for task {task['taskName']}"
                    datas[task["taskID"]] = inFile.Get(histo_name)
                    datas[task["taskID"]].SetDirectory(0)  # detach from file
                    bar()

    fitters = []
    with alive_bar(len(tasks), title="Building fitters") as bar:
        for task in tasks:
            fitter = build_raw_yield_fitter(task, config["fitConfig"])
            fitters.append(fitter)
            fitter.set_data_to_fit_hist(datas[task["taskID"]])
            fitter.setup()
            bar.text = f"Prepared fitter for task {task['taskName']}"
            bar()

    pdfs = {}
    figs = {}
    with alive_bar(len(fitters), title="Running fits") as bar:
        with ThreadPoolExecutor(max_workers=6) as executor:
            futures = {executor.submit(fitter.fit): fitter for fitter in fitters}
            for future in as_completed(futures):
                fitter = futures[future]
                try:
                    future.result()
                    bar.text = f"Completed fit for task {fitter.fit_name}"
                except Exception as exc:
                    bar.text = f"Fit generated an exception for task {fitter.fit_name}: {exc}"
                bar()
    
    # with alive_bar(len(fitters), title="Running fits") as bar:
    #     for fitter in fitters:
    #         fitter.fit()
    #         bar.text = f"Completed fit for task {fitter.name}"
    #         bar()
    fitters.sort(key=lambda x: (int(x.fit_name.split("_")[0])))
    pdfs = []
    for i, fitter in enumerate(fitters):
        task_id = int(fitter.fit_name.split("_")[0])
        task = tasks[task_id]
        iPtCand = task['iPtCand']
        iPtHad = task['iPtHad']
        iDeltaPhi = task['iDeltaPhi']
        bin_idx = task['iDeltaPhi'] + 1

        fig = fitter.plot_fit(False, True, loc=["lower left", "upper left"])
        fig.savefig(outFilePath / tasks[int(fitter.fit_name.split("_")[0])]['pdfPath'], dpi=300, bbox_inches="tight")
        fit_info, _, _, _, _ = fitter.get_fit_info()
        hPairsYields_vs_DeltaPhi[(iPtCand, iPtHad)].SetBinContent(bin_idx, fit_info['sgn']['ry'])
        hPairsYields_vs_DeltaPhi[(iPtCand, iPtHad)].SetBinError(bin_idx, fit_info['sgn']['ry_unc'])

        # fig = fitter.plot_fit(False, True, loc=["lower left", "upper left"])
        # fig.savefig(pdf_path, dpi=300, bbox_inches="tight")s
        # plt.close(fig)
        # figs.append(fig)
        pdfs.append(task['pdfPath'])
        if iDeltaPhi == len(deltaPhiBins) - 2:
            combine_pdfs(pdfs)
            pdfs = []

    outROOT = TFile.Open(str(outFilePath / f"PairYieldsVsPhi.root"), "RECREATE")
    for (iPtCand, iPtHad), histo in hPairsYields_vs_DeltaPhi.items():
        subdir_pt = f"PtCandBin_{int(ptBinsCand[iPtCand]*10)}_{int(ptBinsCand[iPtCand+1]*10)}"
        subdir_had = f"PtHadBin_{int(ptBinsHad[iPtHad]*10)}_{int(ptBinsHad[iPtHad+1]*10)}"
        outROOT.mkdir(subdir_pt + "/" + subdir_had)
        outROOT.cd(subdir_pt + "/" + subdir_had)
        histo.Write()
    outROOT.Close()

    # outFilePath = outdir / f"AssociatedPairsYields"
    # outFilePath.mkdir(parents=True, exist_ok=True)
    # outROOT = TFile.Open(str(outFilePath / f"PairYieldsVsPhi.root"), "RECREATE")
    # outTasks = []
    # for iCand, (ptMin, ptMax) in enumerate(zip(ptBinsCand[:-1], ptBinsCand[1:])):
    #     subdir_pt = f"PtCandBin_{int(ptMin*10)}_{int(ptMax*10)}"
    #     (outFilePath / subdir_pt).mkdir(parents=True, exist_ok=True)
    #     outROOT.mkdir(subdir_pt)
    #     for iHad, (ptHadMin, ptHadMax) in enumerate(zip(ptBinsHad[:-1], ptBinsHad[1:])):
    #         subdir_had = f"PtHadBin_{int(ptHadMin*10)}_{int(ptHadMax*10)}"
    #         (outFilePath / subdir_pt / subdir_had).mkdir(parents=True, exist_ok=True)
    #         outROOT.GetDirectory(subdir_pt).mkdir(subdir_had)
    #     outTasks.append({
    #             'outdir_pdf': outFilePath / subdir_pt / subdir_had / f"PairYields_vs_DeltaPhi.pdf",
    #             'outdir_png': outFilePath / subdir_pt / subdir_had / f"Summary_Fits.png",
    #             'outROOT_subdir': subdir_pt + "/" + subdir_had,
    #             'pdf_canvas': PdfPages(outFilePath / subdir_pt / subdir_had / f"PairYields_vs_DeltaPhi.pdf")
    #         })

    # for i, task in enumerate(tasks):
    #     fitter = fitters[i]
    #     iPtCand = task['iPtCand']
    #     iPtHad = task['iPtHad']
    #     if (fitter.fit_name != f"{task['taskID']}_{task['taskName']}"):
    #         raise ValueError(f"Fitter name {fitter.fit_name} does not match task name {task['taskName']}")
    #     if task['iDeltaPhi'] == 0:
    #         # outPDF = outTasks[iPtCand * len(ptBinsHad[:-1]) + iPtHad]['pdf']
    #         # pngTemp = []
    #         hTemp = TH1D(f"hPairsYields_vs_DeltaPhi", "Associated Pair Raw Yields vs #Delta#phi;#Delta#phi (rad); Associated pairs raw yield", len(deltaPhiBins)-1, np.array(deltaPhiBins, dtype='double'))
    #     fig = fitter.plot_fit(False, True, loc=["lower left", "upper left"])
    #     current_task = outTasks[iPtCand * len(ptBinsHad[:-1]) + iPtHad]
    #     current_task['pdf_canvas'].savefig(fig)
    #     bin_idx = task['iDeltaPhi'] + 1
    #     fit_info, _, _, _, _ = fitter.get_fit_info()
    #     hTemp.SetBinContent(bin_idx, fit_info['sgn']['ry'])
    #     hTemp.SetBinError(bin_idx, fit_info['sgn']['ry_unc'])
    #     if task['iDeltaPhi'] == len(deltaPhiBins) - 2:
    #             current_task['pdf_canvas'].close()
    #             outROOT.cd(current_task['outROOT_subdir'])
    #             hTemp.Write()
    # outROOT.Close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Raw Yield Fitter Interface")
    parser.add_argument('config', type=str, help='Path to the configuration YAML file')
    args = parser.parse_args()

    config_path = args.config
    interface_raw_yield_fitter(config_path)