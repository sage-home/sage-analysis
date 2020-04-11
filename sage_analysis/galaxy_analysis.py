"""

TODO
----
Split this up into smaller functions. Its ghastly as is.
Doc showing exactly how to generate non-default calc/plot functions.
Doc showing exactly how to generate non-default extra properties.
"""

import logging
import os
from typing import List, Dict, Union, Optional, Any, Tuple, Callable

import sage_analysis.example_calcs
import sage_analysis.example_plots

from sage_analysis.utils import generate_func_dict, read_generic_sage_params
from sage_analysis.model import Model
from sage_analysis.sage_binary import SageBinaryData
try:
    from sage_analysis.sage_hdf5 import SageHdf5Data
except ImportError:
    print("h5py not found.  If you're reading in HDF5 output from SAGE, please install this package.")

import matplotlib
import numpy as np


logger = logging.getLogger(__name__)


class GalaxyAnalysis:

    def __init__(
        self,
        sage_parameter_fnames: List[str],
        snapshots_to_plot: Optional[List[int]] = None,
        IMFs: Optional[List[str]] = None,
        labels: Optional[List[str]] = None,
        sage_output_formats: Optional[List[str]] = None,
        first_files_to_analyse: Optional[List[int]] = None,
        last_files_to_analyse: Optional[List[int]] = None,
        num_sage_output_files: Optional[List[int]] = None,
        output_format_data_classes: Optional[Dict[str, Any]] = None,
        random_seeds: Optional[List[int]] = None,
        plot_toggles: Optional[Dict[str, bool]] = None,
        calculation_functions: Optional[Dict[str, Tuple[Callable, Dict[str, Any]]]] = None,
        plot_functions: Optional[Dict[str, Tuple[Callable, Dict[str, Any]]]] = None,
        galaxy_properties_to_analyse: Optional[Dict[str, Dict[str, Union[str, List[str]]]]] = None,
        plot_output_format: str = "png",
        plot_output_path: str = "./plots",
    ):
        """
        snapshots_to_plot : list of ints, optional
            If not specified, uses the lowest redshift snapshot as specified by the redshift file read from ``sage_file``.

        first_files_to_analyse, last_files_to_analyse : list of ints, optional
            The output **SAGE** files to be analysed.  This is an inclusive range, with the output files analyzed ranging
            from ``[first_file_to_analyse, last_file_to_analyse]`` for each model.

            If the corresponding entry in ``sage_output_format`` is ``sage_binary`` (whether passed explicitly or read from
            ``sage_file``), these two variables **MUST** be specified.  Otherwise, if not specified, will analyse **ALL**
            output HDF5 files.

        calculation_functions : dict [string, tuple(function, dict[string, variable])], optional
            A dictionary of functions that are used to compute the properties of galaxies being analysed.  Here, the string
            is the name of the function (e.g., ``"calc_SMF"``), the value is a tuple containing the function itself (e.g.,
            ``calc_SMF()``), and another dictionary which specifies any optional keyword arguments to that function with
            keys as the name of variable (e.g., ``"calc_sub_populations"``) and values as the variable value (e.g.,
            ``True``).

            The functions in this dictionary are called for all files analyzed and **MUST** have a signature ``func(Model,
            gals, optional_keyword_arguments)``. This dict can be generated using
            :py:func:`~sage_analysis.utils.generate_func_dict`.

            If not specified, will use the functions found in :py:module:`~sage_analysis.example_calcs`, filtered to ensure
            that only those functions necessary to plot the plots specified by ``plot_toggles`` are run.

        plot_functions : dict [string, tuple(function, dict[string, variable])], optional
            A dictionary of functions that are used to plot the properties of galaxies being analysed.  Here, the string is
            the name of the function (e.g., ``"plot_SMF"``), the value is a tuple containing the function itself (e.g.,
            ``plot_SMF()``), and another dictionary which specifies any optional keyword arguments to that function with
            keys as the name of variable (e.g., ``"plot_sub_populations"``) and values as the variable value (e.g.,
            ``True``).

            The functions in this dictionary are called for all files analyzed and **MUST** have a signature ``func(Models,
            plot_output_path, plot_output_format, optional_keyword_arguments)``. This dict can be generated using
            :py:func:`~sage_analysis.utils.generate_func_dict`.

            If not specified, will use the functions found in :py:module:`~sage_analysis.example_plots`, filtered to ensure
            that only those functions necessary to plot the plots specified by ``plot_toggles`` are run.
        """

        self._plot_output_path = plot_output_path
        if not os.path.exists(plot_output_path):
            os.makedirs(plot_output_path)

        self._plot_output_format = plot_output_format

        num_models = len(sage_parameter_fnames)
        self._num_models = num_models

        if snapshots_to_plot is None:
            snapshots_to_plot = [None] * num_models
        self._snapshots_to_plot = snapshots_to_plot

        if IMFs is None:
            IMFs = ["Chabrier"] * num_models
        self.IMFs = IMFs

        if labels is None:
            labels = [None] * num_models
        self._labels = labels

        if sage_output_formats is None:
            sage_output_formats = [None] * num_models
        self._sage_output_formats = sage_output_formats

        if first_files_to_analyse is None:
            first_files_to_analyse = [None] * num_models
        self._first_files_to_analyse = first_files_to_analyse

        if last_files_to_analyse is None:
            last_files_to_analyse = [None] * num_models
        self._last_files_to_analyse = last_files_to_analyse

        if num_sage_output_files is None:
            num_sage_output_files = [None] * num_models
        self._num_sage_output_files = num_sage_output_files

        if output_format_data_classes is None:
            output_format_data_classes = {"sage_binary": SageBinaryData, "sage_hdf5": SageHdf5Data}
        self._output_format_data_classes = output_format_data_classes

        if random_seeds is None:
            random_seeds = [None] * num_models
        self._random_seeds = random_seeds

        parameters = [
            sage_parameter_fnames,
            snapshots_to_plot,
            IMFs,
            labels,
            sage_output_formats,
            first_files_to_analyse,
            last_files_to_analyse,
            num_sage_output_files,
            random_seeds,
        ]

        # All the parameters should have the same lengths.
        for parameter_vals in parameters:
            if len(parameter_vals) != self._num_models:
                raise ValueError(
                    f"The number of parameter files is {self._num_models}. Ensure that all parameters passed to "
                    f"``analyse_sage_output`` are lists of this length."
                )

        if plot_toggles is None:
            plot_toggles = {
                "SMF" : True,  # Stellar mass function.
                "BMF" : True,  # Baryonic mass function.
                "GMF" : True,  # Gas mass function (cold gas).
                "BTF" : True,  # Baryonic Tully-Fisher.
                "sSFR" : True,  # Specific star formation rate.
                "gas_fraction" : True,  # Fraction of galaxy that is cold gas.
                "metallicity" : True,  # Metallicity scatter plot.
                "bh_bulge" : True,  # Black hole-bulge relationship.
                "quiescent" : True,  # Fraction of galaxies that are quiescent.
                "bulge_fraction" : True,  # Fraction of galaxies that are bulge/disc dominated.
                "baryon_fraction" : True,  # Fraction of baryons in galaxy/reservoir.
                "reservoirs" : True,  # Mass in each reservoir.
                "spatial" : True   # Spatial distribution of galaxies.
            }
        self._plot_toggles = plot_toggles

        # ``parameters`` is a matrix of parameters with each "column" specifying the parameters for a single model. Hence
        # we want to iteratre through column-wise and use these to build the ``Model`` class instance. Here, the ``map``
        # function does this tranpose into a column-wise iterable.
        models = [
            Model(*model_parameters, plot_toggles) for model_parameters in map(list, zip(*parameters))
        ]
        self._models = models

        # Then populate the `calculation_methods` dictionary. This dictionary will control
        # which properties each model will calculate.  The dictionary is populated using
        # the plot_toggles defined above.
        # Our functions are inside the `example_calcs.py` module and are named "calc_<toggle>". If
        # your functions are in a different module or different function prefix, change it
        # here.
        # ALL FUNCTIONS MUST HAVE A FUNCTION SIGNATURE `func(Model, gals, optional_kwargs=...)`.
        if calculation_functions is None:
            calculation_functions = generate_func_dict(
                plot_toggles, module_name="sage_analysis.example_calcs", function_prefix="calc_"
            )

        if plot_functions is None:
            plot_functions = generate_func_dict(
                plot_toggles, module_name="sage_analysis.example_plots", function_prefix="plot_"
            )
        self._plot_functions = plot_functions

        if galaxy_properties_to_analyse is None:
            galaxy_properties_to_analyse = {
                "stellar_mass_bins": {
                    "type": "binned",
                    "bin_low": 8.0,
                    "bin_high": 12.0,
                    "bin_width": 0.1,
                    "property_names": [
                        "SMF", "red_SMF", "blue_SMF", "BMF", "GMF",
                        "centrals_MF", "satellites_MF", "quiescent_galaxy_counts",
                        "quiescent_centrals_counts", "quiescent_satellites_counts",
                        "fraction_bulge_sum", "fraction_bulge_var",
                        "fraction_disk_sum", "fraction_disk_var"
                    ],
                },
                "halo_mass_bins": {
                    "type": "binned",
                    "bin_low": 10.0,
                    "bin_high": 14.0,
                    "bin_width": 0.1,
                    "property_names": ["fof_HMF"] + [f"halo_{component}_fraction_sum"
                        for component in ["baryon", "stars", "cold", "hot", "ejected", "ICS", "bh"]
                    ],
                },
                "scatter_properties": {
                    "type": "scatter",
                    "property_names": [
                        "BTF_mass", "BTF_vel", "sSFR_mass", "sSFR_sSFR",
                        "gas_frac_mass", "gas_frac", "metallicity_mass",
                        "metallicity", "bh_mass", "bulge_mass", "reservoir_mvir",
                        "reservoir_stars", "reservoir_cold", "reservoir_hot",
                        "reservoir_ejected", "reservoir_ICS", "x_pos",
                        "y_pos", "z_pos"
                    ],
                },
            }
        self._galaxy_properties_to_analyse = galaxy_properties_to_analyse

        for model in self._models:

            # Read the parameter files and update any of the missing parameters.
            self._read_sage_file(model)

            # Also initialise all of the properties that are required.
            for name, galaxy_properties in self._galaxy_properties_to_analyse.items():
                self._initialise_properties(name, model, galaxy_properties)

            # Finally, set some attributes that are more suited to a per-model basis.
            model._calculation_functions = calculation_functions


    def _read_sage_file(self, model: Model) -> None:

        # If the format wasn't defined, then attempt to read a default parameter file to determine format.
        if model.sage_output_format is None:
            logger.info(
                f"No SAGE output format specified. Attempting to read ``{model.sage_file}`` and using the format "
                f"specified inside."
            )
            sage_dict = read_generic_sage_params(model.sage_file)

            model.sage_output_format = sage_dict["_output_format"]
            logger.info(f"Using ``{model.sage_output_format}`` output format.")

        # Each SAGE output has a specific class written to read in the data.
        model.data_class = self._output_format_data_classes[model.sage_output_format](model, model.sage_file)

        # The data class has read the SAGE ini file.  Update the model with the parameters
        # read and those specified by the user. We will also log some of these.
        for key, value in model.data_class.sage_model_dict.items():

            # Check if the attribute has already been set to a non-default value.
            try:
                attr_value = getattr(model, key)
            except AttributeError:
                pass
            else:
                if attr_value is not None:
                    continue
            setattr(model, key, value)

            default_messages = {
                "_snapshot": f"Snapshot to analyse not specified; using the final snapshot of the simulation ({value})",
                "_label": f"Label not specified; using the FileNameGalaxies from parameter file ({value})",
                "_first_file_to_analyse": f"First file to analyse not specified; using {value}",
                "_last_file_to_analyse": f"Last file to analyse not specified; using num cores SAGE ran with minus 1 ({value})",
            }

            try:
                logger.info(default_messages[key])
            except KeyError:
                pass

        model.volume = model.data_class.determine_volume_analysed(model)

        # Some properties require the stellar mass function to normalize their values. For
        # these, the SMF plot toggle is explicitly required.
        try:
            if model.plot_toggles["SMF"]:
                model.calc_SMF = True
            else:
                model.calc_SMF = False
        except KeyError:  # Maybe we've removed "SMF" from plot_toggles...
                model.calc_SMF = False

    def _initialise_properties(
        self,
        name: str,
        model: Model,
        galaxy_properties: Dict[str, Dict[str, Union[str, List[str]]]],
    ) -> None:

            # Properties can be binned (e.g., how many galaxies with mass between 10^8.0
            # and 10^8.1), scatter plotted (e.g., for 1000 galaxies plot the specific star
            # formation rate versus stellar mass) or a single number (e.g., the sum
            # of the star formation rate at a snapshot). Properties can be accessed using
            # `Model.properties["property_name"]`; e.g., `Model.properties["SMF"]`.

            # First let's do the properties binned on stellar mass. The bins themselves can be
            # accessed using `Model.bins["bin_name"]`; e.g., `Model.bins["stellar_mass_bins"]

            allowed_property_types = ["binned", "scatter", "single"]
            if galaxy_properties["type"] not in allowed_property_types:
                raise ValueError(
                    f"Requested to analyse a galaxy property with unkown type.  The galaxy properties were "
                    f"{galaxy_properties} and the only accepted types are {allowed_property_types}."
                )

            if galaxy_properties["type"] == "binned":
                model.init_binned_properties(
                    galaxy_properties["bin_low"],
                    galaxy_properties["bin_high"],
                    galaxy_properties["bin_width"],
                    name,
                    galaxy_properties["property_names"],
                )
            elif galaxy_properties["type"] == "scatter":
                model.init_scatter_properties(galaxy_properties["property_names"])
            elif galaxy_properties["type"] == "single":
                model.init_single_properties(galaxy_properties["property_names"])

            logger.info(f"Initialized galaxy properties {galaxy_properties}")

    def analyse_galaxies(self):

        for model in self._models:

            # To be more memory concious, we calculate the required properties on a
            # file-by-file basis. This ensures we do not keep ALL the galaxy data in memory.
            model.calc_properties_all_files(model._calculation_functions, model._snapshot, debug=False)

    def generate_plots(self) -> List[matplotlib.figure.Figure]:

        # Now do the plotting.
        figs = []
        for func, kwargs in self._plot_functions.values():
            fig = func(self._models, self._plot_output_path, self._plot_output_format, **kwargs)

            if type(fig) == list:
                figs.extend(fig)
            else:
                figs.append(fig)

        return figs
