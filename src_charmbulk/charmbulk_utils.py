import os
from pathlib import Path
import sys
sys.path.append(str(Path(__file__).parent.parent / "utils"))
from utils import get_paths as _get_paths

def get_paths():
    return _get_paths()

def get_pt_dependent_param(param, nPtBins, isList=False):
    ''' handle parameters that can be either a single value or a list of values for each pt bin.
    - If `isList` is **True**, the parameter is expected to be a `list` of `list` for each pt bin. 
    - If `isList` is **False**, the parameter can be a list of values or a single value to be repeated.
    '''
    def _repeat_last_value(list_param):
        if len(list_param) < nPtBins:
            for i in range(len(list_param), nPtBins):
                list_param.append(list_param[-1])  # repeat last value if not enough values provided
        return list_param
    if isList and isinstance(param, list):
        _repeat_last_value(param)
        parameter = param
    else:
        if isinstance(param, list):
            _repeat_last_value(param)
            parameter = param
        else:
            parameter = [param] * nPtBins
    return parameter