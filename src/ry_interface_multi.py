import argparse
import yaml
import numpy as np
import pathlib as PATH
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from ROOT import TFile, TH1D
from alive_progress import alive_bar
import os
os.environ['ZFIT_DISABLE_TF_WARNINGS'] = '1'
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
import tensorflow as tf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing
tf.config.threading.set_intra_op_parallelism_threads(2)
tf.config.threading.set_inter_op_parallelism_threads(1)
import sys
sys.path.append("../flareflyfitter/")
from raw_yield_fitter import RawYieldFitter
import fitz
from PyPDF2 import PdfMerger


def build_raw_yield_fitter(task, fit_config):
    task_id = task['taskID']
    task_name = task['taskName']
    pt_min = task['ptMin']
    pt_max = task['ptMax']
    pt_str = f"pt_{int(pt_min*10)}_{int(pt_max*10)}"
    i_pt = task['iPtCand']
    mass_fit_range = fit_config['MassFitRanges'][i_pt]
    minimizer = fit_config.get('minimizer', 'flarefly')

    sgn_func = {
        'func': fit_config.get('SgnFunc')[i_pt] if isinstance(fit_config.get('SgnFunc'), list) else fit_config.get('SgnFunc'),
        'label': 'sgn',
        'particle': task['Dmeson'],
    }
    bkg_func = {
        'func': fit_config.get('BkgFunc')[i_pt] if isinstance(fit_config.get('BkgFunc'), list) else fit_config.get('BkgFunc'),
        'label': 'Comb_bkg',
    }

    rebin = fit_config.get('Rebin')[i_pt] if isinstance(fit_config.get('Rebin'), list) else fit_config.get('Rebin', 1)
    fit_pars = fit_config.get('InitPars', None)

    fitter = RawYieldFitter(task['Dmeson'], pt_min, pt_max, pt_str, minimizer)
    fitter.set_fit_range(mass_fit_range[0], mass_fit_range[1])
    fitter.add_sgn_func(sgn_func['func'], sgn_func['label'], sgn_func['particle'])
    fitter.add_bkg_func(bkg_func['func'], bkg_func['label'])
    fitter.set_name(f'{task_id}_{task_name}')
    fitter.set_rebin(rebin)
    if fit_pars:
        fitter.set_fit_pars(fit_pars, pt_min, pt_max)

    return fitter


def combine_pdfs(pdf_paths):
    merger = PdfMerger()
    sum_pdf = pdf_paths[0].parent / 'Summary_Fits.pdf'
    for pdf_path in pdf_paths:
        merger.append(str(pdf_path))
    merger.write(str(sum_pdf))
    merger.close()

    for i_total in range(int(np.ceil(len(pdf_paths) / 15))):
        tot, axes = plt.subplots(3, 5, figsize=(20, 12))
        plt.subplots_adjust(wspace=0.01, hspace=0.01, left=0.02, right=0.98, top=0.98, bottom=0.02)
        axes = axes.flatten()
        chunk = pdf_paths[i_total * 15:(i_total + 1) * 15]
        for i, pdf_path in enumerate(chunk):
            doc = fitz.open(str(pdf_path))
            page = doc.load_page(0)
            pix = page.get_pixmap(dpi=300)
            img = np.frombuffer(pix.samples, dtype=np.uint8)
            img = img.reshape(pix.height, pix.width, pix.n)
            axes[i].imshow(img)
            axes[i].axis('off')
            doc.close()

        for j in range(i + 1, 15):
            axes[j].axis('off')

        plt.tight_layout()
        tot_path = sum_pdf.parent / f'Summary_Fits_Combined_{i_total}'
        plt.savefig(tot_path.with_suffix('.pdf'), dpi=300, bbox_inches='tight')
        plt.savefig(tot_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
        plt.close(tot)


def fit_single_task(task, fit_config, root_file_path):
    root_file = TFile.Open(root_file_path)
    hist = root_file.Get(task['histoPath'])
    if not hist:
        raise FileNotFoundError(f"Histogram {task['histoPath']} not found in {root_file_path}")
    hist.SetDirectory(0)
    root_file.Close()

    fitter = build_raw_yield_fitter(task, fit_config)
    fitter.set_data_to_fit_hist(hist)
    fitter.setup()
    fitter.fit()

    fig = fitter.plot_fit(False, True, loc=['lower left', 'upper left'])
    fig.savefig(str(task['pdfPath']), dpi=300, bbox_inches='tight')
    plt.close(fig)

    fit_info, _, _, _, _ = fitter.get_fit_info()
    return {
        'iPtCand': task['iPtCand'],
        'iPtHad': task['iPtHad'],
        'iDeltaPhi': task['iDeltaPhi'],
        'ry': fit_info['sgn']['ry'],
        'ry_unc': fit_info['sgn']['ry_unc'],
        'pdfPath': str(task['pdfPath']),
    }


def interface_raw_yield_fitter_multi(config_path, max_workers=6):
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    method = config.get('method', '')
    if method != 'DeltaPhiBinning':
        raise ValueError(f"Unsupported method: {method}. Only 'DeltaPhiBinning' is supported.")

    Dmeson = config['Dmeson']
    pt_bins_cand = config['ptBinsCand']
    pt_bins_had = config['ptBinsHad']
    delta_phi_bins = config.get('deltaPhiBins', list(np.linspace(-1.5707963705062866, 4.71238911151886, 32 + 1)))

    outdir = PATH.Path(config['outdir'])
    suffix = config.get('suffix', '')
    out_file_path = outdir / 'AssociatedPairsYields'
    out_file_path.mkdir(parents=True, exist_ok=True)

    root_file_path = outdir / f'CorrelExtract_{suffix}' / 'CorrelationsResults' / 'CorrelationsResults.root'
    if not root_file_path.exists():
        raise FileNotFoundError(f'Could not open {root_file_path}')

    tasks = []
    h_pairs_yields = {}

    for i_pt_cand, (pt_min, pt_max) in enumerate(zip(pt_bins_cand[:-1], pt_bins_cand[1:])):
        for i_pt_had, (pt_had_min, pt_had_max) in enumerate(zip(pt_bins_had[:-1], pt_bins_had[1:])):
            h_temp = TH1D(
                f'hPairsYields_vs_DeltaPhi_{i_pt_cand}_{i_pt_had}',
                'Associated Pair Raw Yields vs #Delta#phi;#Delta#phi (rad); Associated pairs raw yield',
                len(delta_phi_bins) - 1,
                np.array(delta_phi_bins, dtype='double'),
            )
            h_pairs_yields[(i_pt_cand, i_pt_had)] = h_temp
            bin_pdflist_dir = out_file_path / f'PtCandBin_{int(pt_min * 10)}_{int(pt_max * 10)}' / f'PtHadBin_{int(pt_had_min * 10)}_{int(pt_had_max * 10)}'
            bin_pdflist_dir.mkdir(parents=True, exist_ok=True)

            for i_delta_phi, (delta_phi_min, delta_phi_max) in enumerate(zip(delta_phi_bins[:-1], delta_phi_bins[1:])):
                task_id = len(tasks)
                pdf_path = bin_pdflist_dir / f'TempFitResult_DeltaPhiBin_{int(delta_phi_min * 1000)}_{int(delta_phi_max * 1000)}.pdf'
                histo_path = (
                    f"PtCandBin_{int(pt_min * 10)}_{int(pt_max * 10)}"
                    f"/PtHadBin_{int(pt_had_min * 10)}_{int(pt_had_max * 10)}"
                    f"/DeltaPhiBin_{int(delta_phi_min * 1000)}_{int(delta_phi_max * 1000)}"
                    '/hCorrectedMassPairs_AllPools'
                )

                tasks.append({
                    'taskID': task_id,
                    'taskName': f"CandPt_{int(pt_min * 10)}_{int(pt_max * 10)}_HadPt_{int(pt_had_min * 10)}_{int(pt_had_max * 10)}_dPhi_{int(delta_phi_min * 1000)}_{int(delta_phi_max * 1000)}",
                    'iPtCand': i_pt_cand,
                    'ptMin': pt_min,
                    'ptMax': pt_max,
                    'iPtHad': i_pt_had,
                    'ptHadMin': pt_had_min,
                    'ptHadMax': pt_had_max,
                    'iDeltaPhi': i_delta_phi,
                    'deltaPhiMin': delta_phi_min,
                    'deltaPhiMax': delta_phi_max,
                    'Dmeson': Dmeson,
                    'pdfPath': pdf_path,
                    'histoPath': histo_path,
                })

    pdf_registry = defaultdict(dict)
    fit_config = config['fitConfig']
    max_workers = min(max_workers, len(tasks)) if tasks else 1

    with alive_bar(len(tasks), title='Fitting (multi-process)') as bar:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(fit_single_task, task, fit_config, str(root_file_path)): task for task in tasks}
            for future in as_completed(futures):
                result = future.result()
                key = (result['iPtCand'], result['iPtHad'])
                hist = h_pairs_yields[key]
                bin_idx = result['iDeltaPhi'] + 1
                hist.SetBinContent(bin_idx, result['ry'])
                hist.SetBinError(bin_idx, result['ry_unc'])
                pdf_registry[key][result['iDeltaPhi']] = PATH.Path(result['pdfPath'])
                bar()

    for key, pdf_map in pdf_registry.items():
        sorted_paths = [pdf_map[phi] for phi in sorted(pdf_map)]
        if sorted_paths:
            combine_pdfs(sorted_paths)

    out_root = TFile.Open(str(out_file_path / 'PairYieldsVsPhi.root'), 'RECREATE')
    for (i_pt_cand, i_pt_had), hist in h_pairs_yields.items():
        subdir_pt = f'PtCandBin_{int(pt_bins_cand[i_pt_cand] * 10)}_{int(pt_bins_cand[i_pt_cand + 1] * 10)}'
        subdir_had = f'PtHadBin_{int(pt_bins_had[i_pt_had] * 10)}_{int(pt_bins_had[i_pt_had + 1] * 10)}'
        out_root.mkdir(f'{subdir_pt}/{subdir_had}')
        out_root.cd(f'{subdir_pt}/{subdir_had}')
        hist.Write('hPairsYields_vs_DeltaPhi')
    out_root.Close()


def main():
    multiprocessing.set_start_method('spawn', force=True)
    parser = argparse.ArgumentParser(description='Multi-process Raw Yield Fitter Interface')
    parser.add_argument('config', type=str, help='Path to the configuration YAML file')
    parser.add_argument('--workers', type=int, default=6, help='Number of worker processes to use')
    args = parser.parse_args()
    interface_raw_yield_fitter_multi(args.config, max_workers=args.workers)


if __name__ == '__main__':
    main()
