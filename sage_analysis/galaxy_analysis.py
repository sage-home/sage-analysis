"""

TODO
----
Split this up into smaller functions. Its ghastly as is.
Doc showing exactly how to generate non-default calc/plot functions.
Doc showing exactly how to generate non-default extra properties.
"plot_output_path" should be passed to "generate_plots".
Add more nuance to checking when/if SMF needs to be calculated.
Add "*args" to data class functions to allow them to pass a variable number of arguments.
Add extra tests to ensure that things are being logged when there are no galaxies present.
"""

import logging
import os
from typing import List, Dict, Union, Optional, Any, Tuple, Callable

import sage_analysis.example_calcs
import sage_analysis.example_plots

from sage_analysis.default_analysis_arguments import default_plot_toggles, default_galaxy_properties_to_analyse
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
        IMFs: Optional[List[str]] = None,
        labels: Optional[List[str]] = None,
        sage_output_formats: Optional[List[str]] = None,
        first_files_to_analyse: Optional[List[int]] = None,
        last_files_to_analyse: Optional[List[int]] = None,
        num_sage_output_files: Optional[List[int]] = None,
        output_format_data_classes: Optional[Dict[str, Any]] = None,
        random_seeds: Optional[List[int]] = None,
        history_redshifts: Optional[Dict[str, Union[List[float], str]]] = None,
        plot_toggles: Optional[Dict[str, bool]] = None,
        history_plot_toggles_keys: Optional[List[str]] = None,
        calculation_functions: Optional[Dict[str, Tuple[Callable, Dict[str, Any]]]] = None,
        plot_functions: Optional[Dict[str, Tuple[Callable, Dict[str, Any]]]] = None,
        galaxy_properties_to_analyse: Optional[Dict[str, Dict[str, Union[str, List[str]]]]] = None,
        plots_that_need_smf: Optional[List[str]] = None,
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
            is the name of the plot toggle (e.g., ``"SMF"``), the value is a tuple containing the function itself (e.g.,
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

        plots_that_need_smf : list of strings, optional
            The plot toggles that require the stellar mass function to be properly computed and analysed. For example,
            plotting the quiescent fraction of galaxies requires knowledge of the total number of galaxies. The strings
            here must **EXACTLY** match the keys in ``plot_toggles``.

            If not specified, uses a default value of ``["SMF", "quiescent", "bulge_fraction", "SMF_history"]``.
        """

        self._plot_output_path = plot_output_path
        if not os.path.exists(plot_output_path):
            os.makedirs(plot_output_path)

        self._plot_output_format = plot_output_format

        num_models = len(sage_parameter_fnames)
        self._num_models = num_models

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

        if history_redshifts is None:
            history_redshifts = {
                "SMF_history": "All",
                "SMD_history": "All",
                "SFRD_history": "All",
            }
        self._history_redshifts = history_redshifts

        if plots_that_need_smf is None:
            plots_that_need_smf = ["SMF", "quiescent", "bulge_fraction", "SMF_history"]

        individual_model_parameters = [
            sage_parameter_fnames,
            IMFs,
            labels,
            sage_output_formats,
            first_files_to_analyse,
            last_files_to_analyse,
            num_sage_output_files,
            random_seeds,
        ]

        global_model_parameters = [
            history_redshifts,
            plots_that_need_smf,
        ]

        if plot_toggles is None:
            plot_toggles = default_plot_toggles
        self._plot_toggles = plot_toggles

        # ``parameters`` is a matrix of parameters with each "column" specifying the parameters for a single model. Hence
        # we want to iteratre through column-wise and use these to build the ``Model`` class instance. Here, the ``map``
        # function does this tranpose into a column-wise iterable.
        models = [
            Model(*model_parameters, *global_model_parameters, plot_toggles) for model_parameters in map(list, zip(*individual_model_parameters))
        ]
        self._models = models

        # Determine if the stellar mass function needs to be computed for each model. Important we do this before
        # checking ``calculation_functions``.
        for model in models:
            if self._does_smf_need_computing(model):
                model.plot_toggles["SMF"] = True
                plot_toggles["SMF"] = True

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

        # Because we have adjusted the SMF, check if it's a calculation.
        for model in models:
            # Condition checks for condition where SMF is being compputed but not in ``calculation_functions``.
            if model._plot_toggles.get("SMF", None) and not calculation_functions.get("SMF", None):
                raise ValueError(
                    "The stellar mass function (``SMF``) is being computed, either because it was set manually "
                    "(through ``plot_toggles``) or because another plot requires it (``plots_that_need_smf``).\n"
                    "However, ``calc_SMF`` was not found in ``calculation_functions``. Ensure that it is added."
                )

        print(calculation_functions)

        if plot_functions is None:
            plot_functions = generate_func_dict(
                plot_toggles, module_name="sage_analysis.example_plots", function_prefix="plot_"
            )
        self._plot_functions = plot_functions

        if galaxy_properties_to_analyse is None:
            galaxy_properties_to_analyse = default_galaxy_properties_to_analyse
        self._galaxy_properties_to_analyse = galaxy_properties_to_analyse

        for model in self._models:

            # Read the parameter files and update any of the missing parameters.
            self._read_sage_file(model)

            # Also initialise all of the properties that are required.
            for name, galaxy_properties in self._galaxy_properties_to_analyse.items():
                for snapshot, _ in enumerate(model._redshifts):
                    self._initialise_properties(name, model, galaxy_properties, snapshot)

            # Finally, set some attributes that are more suited to a per-model basis.
            model._calculation_functions = calculation_functions

            # Go through the calculation functions and pull out those that are actually being analysed over a number of
            # snapshots. Ensure these aren't in ``_calculation_functions`` because otherwise we'd double count.
            history_calculation_functions = {}
            for func_name in self._history_redshifts.keys():
                try:
                    calculation_function = calculation_functions[func_name]
                except KeyError:
                    continue
                history_calculation_functions[func_name] = calculation_function
                del model._calculation_functions[func_name]

            model._history_calculation_functions = history_calculation_functions

            self._set_history_snapshots(model)

            # Based on the snapshots required to analyse over history, we may have to loop over different snapshots.
            history_snaps_to_loop = []
            for property_name in self._history_redshifts.keys():

                # If the plot toggles have been changed, then there's no guarantee that the keys for
                # ``_history_redshifts`` and ``_plot_toggles`` match.
                try:
                    plot_toggle_value = self._plot_toggles[property_name]
                except KeyError:
                    continue

                # Furthermore, maybe plotting has been disabled for this property.
                if not plot_toggle_value:
                    continue

                snaps = getattr(model, f"_history_{property_name}_snaps")
                history_snaps_to_loop.extend(snaps)

            model._history_snaps_to_loop = np.unique(history_snaps_to_loop)
            logger.info(f"Looping through snapshots {model._history_snaps_to_loop}")


    def _set_history_snapshots(self, model: Model) -> None:

        # Maybe there were no history redshifts specified.
        if self._history_redshifts is None:
            self._history_snapshots = None
            return

        # Convert these redshifts into snapshots.
        for property_name, property_redshifts in self._history_redshifts.items():

            if property_redshifts == "All":
                redshifts = model._redshifts
            else:
                redshifts = property_redshifts

            attrname = f"_history_{property_name}_redshifts"
            setattr(model, attrname, redshifts)

            # Find the snapshots that are closest to the requested redshifts.
            property_snaps = [(np.abs(model._redshifts - redshift)).argmin() for redshift in redshifts]

            attrname = f"_history_{property_name}_snaps"
            setattr(model, attrname, property_snaps)
        return

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

    def _initialise_properties(
        self,
        name: str,
        model: Model,
        galaxy_properties: Dict[str, Dict[str, Union[str, List[str]]]],
        snapshot: int,
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
                    snapshot,
                )
            elif galaxy_properties["type"] == "scatter":
                model.init_scatter_properties(galaxy_properties["property_names"], snapshot)
            elif galaxy_properties["type"] == "single":
                model.init_single_properties(galaxy_properties["property_names"], snapshot)

            logger.info(f"Initialized galaxy properties {galaxy_properties} for Snapshot {snapshot}")


    def _does_smf_need_computing(self, model: Model) -> bool:

        # Maybe the SMF has already been toggled on.
        try:
            toggle = model._plot_toggles["SMF"]
        except KeyError:
            # Perhaps "SMF" isn't present in the plot toggles.
            pass
        else:
            # "SMF" is present and is turned on -> SMF does need to be computed.
            if toggle:
                return True

        # Determine those plots that are being computed.
        plots = [toggle for toggle, value in model._plot_toggles.items() if value]

        # Then, check if any of these plots need the SMF.
        if any([plot in model._plots_that_need_smf for plot in plots]):
            logger.info(f"One of your plots require the SMF to be calculated. Turning the SMF plot toggle on.")
            return True

        # Otherwise, they don't need the SMF.
        return False


    def analyse_galaxies(self, snapshot: Optional[int] = None):

        for model in self._models:

            if snapshot is None:
                snapshot = len(model._redshifts) - 1

            model.calc_properties_all_files(model._calculation_functions, snapshot, debug=False)

            # Now handle calculate properties that are tracked over redshift (if any).
            if model._history_snaps_to_loop is None:
                continue

            for snap in model._history_snaps_to_loop:
                model.calc_properties_all_files(
                    model._history_calculation_functions, snap, debug=False, close_file=False
                )

    def generate_plots(self, snapshot: Optional[int] = None) -> List[matplotlib.figure.Figure]:

        # Snapshot needs to be a list to allow plotting different snapshots for different models.
        if snapshot is None:
            snapshot = len(self._models[0]._redshifts) - 1

        # Now do the plotting.
        figs = []

        for func, kwargs in self._plot_functions.values():
            fig = func(self._models, snapshot, self._plot_output_path, self._plot_output_format, **kwargs)

            if type(fig) == list:
                figs.extend(fig)
            else:
                figs.append(fig)

        return figs
