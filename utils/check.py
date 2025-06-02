'''
Utility functions to check the validity of configuration or file
'''
import sys
import os

def check_ptbinned_par(parameter, nPtBins: int):
    """
    Check if the pt-dependent parameters are valid.
    
    Parameters
    -----------
    parameter: list or single value
        The parameter to check. If a single value is provided, it will be replicated for each pt bin.
    nPtBins: int
        The number of pT bins for which the parameter should be defined.

    Returns
    -----------
    parameter: list
        A list of parameters, one for each pT bin.
    """
    if not isinstance(parameter, list):
        parameter = [parameter] * nPtBins
    else:
        if len(parameter) < nPtBins:
            raise ValueError(f"Parameter {parameter} must be defined for each pT bin or be a single value. "
                             f"Expected {nPtBins} values, got {len(parameter)}.")
    return parameter

def check_ptbinned_pars(nPtBins: int, *parameters):
    """
    Check if the pt-dependent parameters are valid.
    
    Parameters
    -----------
    parameters: tuple of lists or single values
        The parameters to check. Each parameter can be a list of values or a single value.
    nPtBins: int
        The number of pT bins for which the parameters should be defined.
    Returns
    -----------
    parameters: tuple of lists
        A tuple containing the validated parameters, each as a list with one value per pT bin.
    """
    return tuple(check_ptbinned_par(par, nPtBins) for par in parameters)

def check_file_exists(filePath):
    """
    Check if a file exists.
    
    Parameters
    -----------
    file_path: str
        The path to the file to check.
    
    Returns
    -----------
    bool
        True if the file exists, False otherwise.
    """
    if not os.path.isfile(filePath):
        print(f"\033[31mERROR: File {filePath} does not exist!\033[0m")
        sys.exit(1)
    return True