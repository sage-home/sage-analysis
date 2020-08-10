import sys
import logging
import os
from typing import Any, Callable, Dict, Optional, Tuple, List

import numpy as np

logger = logging.getLogger(__name__)


def generate_func_dict(
    plot_toggles,
    module_name,
    function_prefix,
    keyword_args={}
) -> Dict[str, Tuple[Callable, Dict[str, Any]]]:
    """
    Generates a dictionary where the keys are the function name and the value is a list
    containing the function itself (0th element) and keyword arguments as a dictionary
    (1st element). All functions in the returned dictionary are expected to have the same
    call signature for non-keyword arguments. Functions are only added when the
    ``plot_toggles`` value is non-zero.

    Functions are required to be named ``<module_name><function_prefix><plot_toggle_key>``
    For example, the default calculation function are kept in the ``model.py`` module and
    are named ``calc_<toggle>``.  E.g., ``sage_analysis.model.calc_SMF()``,
    ``sage_analysis.model.calc_BTF()``, ``sage_analysis.model.calc_sSFR()`` etc.

    Parameters
    ----------

    plot_toggles: dict, [string, int]
        Dictionary specifying the name of each property/plot and whether the values
        will be generated + plotted. A value of 1 denotes plotting, whilst a value of
        0 denotes not plotting.  Entries with a value of 1 will be added to the function
        dictionary.

    module_name: string
        Name of the module where the functions are located. If the functions are located
        in this module, pass an empty string "".

    function_prefix: string
        Prefix that is added to the start of each function.

    keyword_args: dict [string, dict[string, variable]], optional
        Allows the adding of keyword aguments to the functions associated with the
        specified plot toggle. The name of each keyword argument and associated value is
        specified in the inner dictionary.

    Returns
    -------

    func_dict: dict [string, tuple(function, dict[string, variable])]
        The key of this dictionary is the name of the function.  The value is a list with
        the 0th element being the function and the 1st element being a dictionary of
        additional keyword arguments to be passed to the function. The inner dictionary is
        keyed by the keyword argument names with the value specifying the keyword argument
        value.
    """

    # Check if the specified module is present.
    try:
        module = sys.modules[module_name]
    except KeyError:
        raise KeyError(
            f"Module ``{module_name}`` has not been imported.\nPerhaps you need to create an empty ``__init__.py`` "
            f"file to ensure your package can be imported.\nAlso, ensure ``import {module_name}`` is at the top of "
            f"your script, before ``generate_func_dict`` is called."
        )

    # Only populate those methods that have been marked in the `plot_toggles` dictionary.
    func_dict = {}
    for toggle, value in plot_toggles.items():
        if value:

            func_name = "{0}{1}".format(function_prefix, toggle)

            # Be careful.  Maybe the func for a specified `plot_toggle` value wasn't
            # added to the module.
            try:
                func = getattr(module, func_name)
            except AttributeError:
                raise AttributeError(
                    "Tried to get the func named ``{func_name}`` corresponding to ``plot_toggle`` value ``{toggle}``. "
                    f"However, no func named ``{func_name}`` could be found in ``{module_name}`` module."
                )

            # We may have specified some keyword arguments for this plot toggle. Check.
            try:
                key_args = keyword_args[toggle]
            except KeyError:
                # No extra arguments for this.
                key_args = {}

            func_dict[toggle] = (func, key_args)

    return func_dict


def select_random_indices(
    inds: np.ndarray,
    global_num_inds_available: int,
    global_num_inds_requested: int,
    seed: Optional[int] = None,
) -> np.ndarray:
    """
    Select a random subset of indices if the total number of indices (across all files) is known.  This function is
    used if selecting (e.g.,) 100 galaxies from a sample of 10,000.

    However, if the total number of indices is **NOT** known, then this function is not valid.  For example, if one
    wanted to select 100 spiral galaxies, we may not know how many spiral galaxies are present across all files. In
    such scenarios,
    :py:meth:`~sage_analysis.model.Model.select_random_indices_assumed_equal_distribution` should be used.

    Parameters
    ----------
    vals : :obj:`~numpy.ndarray` of values
        Values that the random subset is selected from.

    global_num_inds_available : int
        The total number of indices available across all files.

    global_num_inds_requested : int
        The total number of indices requested across all files.

    seed : int, optional
        If specified, seeds the random number generator with the specified seed.

    Returns
    -------
    random_inds : :obj:`~numpy.ndarray` of values
        Values chosen.
    """

    if seed is not None:
        np.random.seed(seed)

    # First find out the fraction of value that we need to select.
    num_inds_to_choose = int(len(inds) / global_num_inds_available * global_num_inds_requested)

    # Do we have more values than we need?
    if len(inds) > num_inds_to_choose:
        # Randomly select them.
        random_inds = np.random.choice(inds, size=num_inds_to_choose)
    else:
        # Otherwise, we will just use all the indices we were passed.
        random_inds = inds

    return random_inds


def read_generic_sage_params(sage_file_path: str) -> Dict[str, Any]:
    """
    Reads the **SAGE** parameter file values. This function is used for the default ``sage_binary`` and ``sage_hdf5``
    formats. If you have a custom format, you will need to write a ``read_sage_params`` function in your own data
    class.

    Parameters
    ----------
    sage_file_path: string
        Path to the **SAGE** parameter file.

    Returns
    -------
    model_dict: dict [str, var]
        Dictionary containing the parameter names and their values.

    Errors
    ------
    FileNotFoundError
        Raised if the specified **SAGE** parameter file is not found.
    """

    # Fields that we will be reading from the ini file.
    SAGE_fields = [
        "FileNameGalaxies",
        "OutputDir",
        "FirstFile",
        "LastFile",
        "OutputFormat",
        "NumSimulationTreeFiles",
        "FileWithSnapList",
        "Hubble_h",
        "BoxSize",
        "PartMass"
    ]
    SAGE_dict = {}

    # Ignore lines starting with one of these.
    comment_characters = [";", "%", "-"]

    try:
        with open(sage_file_path, "r") as SAGE_file:
            data = SAGE_file.readlines()

            # Each line in the parameter file is of the form...
            # parameter_name       parameter_value.
            for line in range(len(data)):

                # Remove surrounding whitespace from the line.
                stripped = data[line].strip()

                # May have been an empty line.
                try:
                    first_char = stripped[0]
                except IndexError:
                    continue

                # Check for comment.
                if first_char in comment_characters:
                    continue

                # Split into [name, value] list.
                split = stripped.split()

                # Then check if the field is one we care about.
                if split[0] in SAGE_fields:

                    SAGE_dict[split[0]] = split[1]

    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find SAGE ini file {sage_file_path}")

    # Now we have all the fields, rebuild the dictionary to be exactly what we need for
    # initialising the model.
    model_dict = {}

    model_dict["_label"] = SAGE_dict["FileNameGalaxies"]

    try:
        model_dict["_output_format"] = SAGE_dict["OutputFormat"]
    except KeyError:
        pass

    model_dict["_parameter_dirpath"] = os.path.dirname(sage_file_path)

    # ``FileWithSnapList`` may either be an absolute or relative path (wrt to ``_parameter_dirpath``).
    try:
        fname_absolute = f"{model_dict['_parameter_dirpath']}/{SAGE_dict['FileWithSnapList']}"
        alist = np.loadtxt(fname_absolute)
    except IOError:
        fname_relative = f"{SAGE_dict['FileWithSnapList']}"
        logger.debug(f"Could not find snapshot file {fname_absolute}. Trying as {fname_relative} instead.")
        alist = np.loadtxt(f"{SAGE_dict['FileWithSnapList']}")

    redshifts = 1.0 / alist - 1.0
    model_dict["_redshifts"] = redshifts
    model_dict["_snapshot"] = len(alist) - 1  # By default, plot the final snapshot.

    base_sage_output_path_absolute = f"{model_dict['_parameter_dirpath']}/{SAGE_dict['OutputDir']}/{SAGE_dict['FileNameGalaxies']}"  # noqa: E501
    model_dict["_base_sage_output_path_absolute"] = base_sage_output_path_absolute

    base_sage_output_path_relative = f"{SAGE_dict['OutputDir']}/{SAGE_dict['FileNameGalaxies']}"  # noqa: E501
    model_dict["_base_sage_output_path_relative"] = base_sage_output_path_relative

    model_dict["_output_dir"] = SAGE_dict['OutputDir']
    model_dict["_hubble_h"] = float(SAGE_dict["Hubble_h"])
    model_dict["_box_size"] = float(SAGE_dict["BoxSize"])
    model_dict["_num_sim_tree_files"] = int(SAGE_dict["NumSimulationTreeFiles"])

    return model_dict

def find_closest_indices(values: List[float], target_values: List[float]) -> List[int]:
    """
    Finds the indices in ``values`` that result in values closest to ``target_values``.
    """

    closest_indices = [(np.abs(values - target_value)).argmin() for target_value in target_values]
    return closest_indices
