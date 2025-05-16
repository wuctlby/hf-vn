# Compute pt-weights
**In O2Physics environment**, after **setting the parameters**, one could **run the following command** to compute the pt-weights:
```bash
bash run_ComputePtweights.sh
```
**Parametes**
- `config`: the flow configuration file, following parameters were used: 
        --`Dmeson`
        --`centrality`, **only `centrality` equalling `k3050` or `k6080` is supported**
        --`eff_filename` (used in `get_sparses`)
- `Bspecie`: the specie of B-mother hadron, **only `'Ball'` or `'BsBmix'` is supported**
- `suffix`: suffix to distinguish the anchored pass MC sample, like `'apass3'` or `'apass4'`

Upon completion, a series of output files following the naming convention `pTweight_{Dspecie}_{cent}_{suffix}` will be generated in the directory `./weights/{Dspecie}/{cent}/`. These files include three formats:
- `.root`
- `.png`
- `.pdf`

# Folder structure
ptweights
- models
  - fonll: cross-section
  - tamu: raa in different centrality classes
    - B
    - Dplus
    - Ds
    - Dzero
- weights: output folder, `{Dspecie}/{cent}`

# TODO
The `computePtWeights` script is currently embedded partly within the main script, the configuration file. Next step is separating it alone.

What is required to separate the script:
1. Function of loading thnsparse from on the `.root` file path

    ```python
    if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Arguments')
        parser.add_argument("MC", type=str, help="MC input file")
        parser.add_argument("--Dspecie", "-D", type=str, default='Dzero', help="D meson specie (Dzero, Dplus, Ds)")
        parser.add_argument("--Bspecie", "-B", type=str, default='Ball', help="B meson specie (Ball, BsBmix)")
        parser.add_argument("--centrality", "-C", type=str, default='k3050', help="centrality")
        parser.add_argument("--suffix", "-s", type=str, default='test', help="Suffix for output files")
        parser.add_argument("--fonllD", type=str, 
                            default='./models/fonll/fonll_pythia_beautyFFee_charmhadrons_5dot5tev_y0dot5.root', help="D meson FONLL file")
        parser.add_argument("--fonllB", type=str, 
                            default='./models/fonll/fonll_pythia_beautyFFee_charmhadrons_5dot5tev_y0dot5.root', help="B meson FONLL file")
        parser.add_argument("--Raa", action='store_true', default=False, help="Use RAA applying to FONLL")
        args = parser.parse_args()

        computePtWeights(
            args.MC,
            args.Dspecie,
            args.Bspecie,
            args.centrality,
            args.suffix,
            args.fonllD,
            args.fonllB,
            args.Raa
        )
    ```