#! /bin/bash

# This script computes the ptweights for the given input config.
#! only k3050 and k6080 are supported
ScriptPath=$(dirname $(readlink -f $0))

# filled by hand
export config="" # path to the flow configuration file
export Bspecie='Ball' # 'Ball' or 'BsBmix'
export suffix='' # 'apass3' pr 'apass4'

# Fixed parameters
export fonllD="${ScriptPath}/models/fonll/fonll_pythia_beautyFFee_charmhadrons_5dot5tev_y0dot5.root"
export fonllB="${ScriptPath}/models/fonll/fonll_pythia_beautyFFee_charmhadrons_5dot5tev_y0dot5.root"
export applyRaa=True
export Raa=$([ "$applyRaa" = "False" ] && echo "" || echo "--Raa")

python ${ScriptPath}/ComputePtweights.py \
    ${config} \
    --Bspecie ${Bspecie} \
    --suffix ${suffix} \
    --fonllD ${fonllD} \
    --fonllB ${fonllB} \
    ${Raa}