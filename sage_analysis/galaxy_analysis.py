"""

TODO
----
Split this up into smaller functions. Its ghastly as is.
Doc showing exactly how to generate non-default calc/plot functions.
Doc showing exactly how to generate non-default extra properties.
Add "*args" to data class functions to allow them to pass a variable number of arguments.
Add extra tests to ensure that things are being logged when there are no galaxies present.
Test to plot different snapshots for multiple models.
"""

import logging
import os
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import matplotlib
matplotlib.use('Agg')
import numpy as np

import sage_analysis.example_calcs
import sage_analysis.example_plots  # noqa: F401
from sage_analysis.default_analysis_arguments import (
    default_calculation_functions,
    default_plot_functions,
    default_galaxy_properties_to_analyze,
    default_plot_toggles
)
from sage_analysis.model import Model
from sage_analysis.plot_helper import PlotHelper
from sage_analysis.sage_binary import SageBinaryData
from sage_analysis.utils import generate_func_dict, read_generic_sage_params, find_closest_indices

try:
    from sage_analysis.sage_hdf5 import SageHdf5Data
except ImportError:
    print("h5py not found.  If you're reading in HDF5 output from SAGE, please install this package.")

logger = logging.getLogger(__name__)


class GalaxyAnalysis:
    """
    Handles the ingestion, analysis, and plotting of **SAGE** galaxy outputs.
    """

    def __init__(
        self,
        sage_parameter_fnames: List[str],
        plot_toggles: Optional[Dict[str, bool]] = None,
        sage_output_formats: Optional[List[str]] = None,
        labels: Optional[List[str]] = None,
        first_files_to_analyze: Optional[List[int]] = None,
        last_files_to_analyze: Optional[List[int]] = None,
        num_sage_output_files: Optional[List[int]] = None,
        output_format_data_classes_dict: Optional[Dict[str, Any]] = None,
        random_seeds: Optional[List[int]] = None,
        history_redshifts: Optional[Dict[str, Union[List[float], str]]] = None,
        calculation_functions: Optional[Dict[str, Tuple[Callable, Dict[str, Any]]]] = None,
        plot_functions: Optional[Dict[str, Tuple[Callable, Dict[str, Any]]]] = None,
        galaxy_properties_to_analyze: Optional[Dict[str, Dict[str, Union[str, List[str]]]]] = None,
        plots_that_need_smf: Optional[List[str]] = None,
        IMFs: Optional[List[str]] = None,
    ):
        """
        Parameters
        ----------
        sage_parameter_fnames : list of strings
            The name of the **SAGE** parameter files that are to be analyzed. These are the ``.ini`` files used to
            generate the galaxy files. The length of this variable is equal to the number of models to be analyzed.

        plot_toggles : dict [str, bool], optional
            Specifies which properties should be analyzed and plotted.

            If not specified, uses

            .. highlight:: python
            .. code-block:: python

                default_plot_toggles = {
                    "SMF" : True,
                    "BMF" : True,
                    "GMF" : True,
                    "BTF" : True,
                    "sSFR" : True,
                    "gas_fraction" : True,
                    "metallicity" : True,
                    "bh_bulge" : True,
                    "quiescent" : True,
                    "bulge_fraction" : True,
                    "baryon_fraction" : True,
                    "reservoirs" : True,
                    "spatial" : True,
                    "SMF_history": False,
                    "SFRD_history": False,
                    "SMD_history": False,
                }

        sage_output_formats : list of strings, optional
            The output formats of each **SAGE** model being analyzed.  Each value here **MUST** have a corresponding
            entry in ``output_format_data_classes_dict``. The length of this variable is equal to the number of models
            to be analyzed.

            If not specified, will use the ``OutputFormat`` entry from the respective **SAGE** parameter file.

        labels : list of strings, optional
            The labels to be used in the legend for each model.  The length of this variable is equal to the number of
            models to be analyzed.

            If not specified, will use the ``FileNameGalaxies`` entry from the respective **SAGE** parameter file.

        first_files_to_analyze, last_files_to_analyze : list of ints, optional-ish
            The output **SAGE** files to be analyzed.  This is an inclusive range, with the output files analyzed
            ranging from ``[first_file_to_analyze, last_file_to_analyze]`` for each model. The length of this variable
            is equal to the number of models to be analyzed.

            If the corresponding entry in ``sage_output_format`` is ``sage_binary`` (whether passed explicitly or read
            from ``sage_file``), these two variables **MUST** be specified.  Otherwise, if not specified, will analyze
            **ALL** output HDF5 files.

        num_sage_output_files : list of ints, optional-ish
            Specifies the number of output files that were generated by running **SAGE**. This will generally be equal
            to the number of processors used to run **SAGE** and can be different to the range specified by
            ``[first_file_to_analyze, last_file_to_analyze]``.

            If the corresponding entry in ``sage_output_format`` is ``sage_binary`` (whether passed explicitly or read
            from ``sage_file``), this **MUST** be specified.  Otherwise, this variable is **NOT** used.

        output_format_data_classes_dict : dict [string, class], optional
            A dictionary that maps the output format name to the corresponding data class. Each value in
            ``sage_output_formats`` **MUST** have an entry in this dictionary.

            If not specified, will use a default value
            ``output_format_data_classes_dict = {"sage_binary":`` :py:class:`~sage_analysis.sage_binary.SageBinaryData`
            ``,
            "sage_hdf5":`` :py:class:`~sage_analysis.sage_hdf5.SageHdf5Data` ``}``.

        random_seeds : list of ints, optional
            The values to seed the random number generator for each model.  If the value is ``None``, then the
            generator is seeded using the ``np.random.seed()`` method. The length of this variable is equal to the
            number of models to be analyzed.

            If not specified, uses ``None`` for each model (i.e., no predetermined seed).

        history_redshifts : dict [string, string or list of floats], optional
            Specifies which redshifts should be analyzed for properties and plots that are tracked over time. The keys
            here **MUST** have the same name as in ``plot_toggles``.

            If the value of the entry is ``"All"``, then all snapshots will be analyzed. Otherwise, will search for the
            closest snapshots to the requested redshifts.

            If not specified, uses

            .. highlight:: python
            .. code-block:: python

                history_redshifts = {
                    "SMF_history": "All",
                    "SMD_history": "All",
                    "SFRD_history": "All",
                }

        calculation_functions : dict [string, tuple(function, dict[string, variable])], optional
            A dictionary of functions that are used to compute the properties of galaxies being analyzed.  Here, the
            string is the name of the plot toggle (e.g., ``"SMF"``), the value is a tuple containing the function
            itself (e.g., ``calc_SMF()``), and another dictionary which specifies any optional keyword arguments to
            that function with keys as the name of variable (e.g., ``"calc_sub_populations"``) and values as the
            variable value (e.g., ``True``).

            The functions in this dictionary are called for all files analyzed and **MUST** have a signature
            ``func(model, gals, snapshot, optional_keyword_arguments)``. This dict can be generated using
            :py:func:`~sage_analysis.utils.generate_func_dict`.

            If not specified, will use the functions found in :py:mod:`~sage_analysis.example_calcs`, filtered to
            ensure that only those functions necessary to plot the plots specified by ``plot_toggles`` are run.

        plot_functions : dict [string, tuple(function, dict[string, variable])], optional
            A dictionary of functions that are used to plot the properties of galaxies being analyzed.  Here, the
            string is the name of the function (e.g., ``"plot_SMF"``), the value is a tuple containing the function
            itself (e.g., ``plot_SMF()``), and another dictionary which specifies any optional keyword arguments to
            that function with keys as the name of variable (e.g., ``"plot_sub_populations"``) and values as the
            variable value (e.g., ``True``).

            The functions in this dictionary are called for all files analyzed and **MUST** have a signature
            ``func(models, snapshots, plot_helper, optional_keyword_arguments)``. This dict
            can be generated using :py:func:`~sage_analysis.utils.generate_func_dict`.

            If not specified, will use the functions found in :py:mod:`~sage_analysis.example_plots`, filtered to
            ensure that only those functions necessary to plot the plots specified by ``plot_toggles`` are run.

        galaxy_properties_to_analyze : dict [string, dict[str, float or str or list of strings]], optional
            The galaxy properties that are used when running ``calculation_functions``. The properties initialized here
            will be accessible through ``model.properties["property_name"]``.

            This variable is a nested dictionary with the outer dictionary specifying the name of the bins (if the
            properties will be binned), or a unique name otherwise.

            The inner dictionary has a number of fields that depend upon the type of property.  We support properties
            being either binned against a property (e.g., the stellar or halo mass functions are binned on stellar/halo
            mass), plotted as x-vs-y scatter plots (e.g., specific star formation rate vs stellar mass for 1000
            galaxies), or as a single value (e.g., the stellar mass density).

            For binned against a property, the key/value pairs are:  ``"type": "binned"``, ``bin_low: The lower bound
            of the bin (float)``, ``bin_high: The upper bound of the bin (float)``, ``bin_width: The width of the bin
            (float)``, ``property_names: A list of strings denoting the properties to be initialised``. The bin values
            are all initialized as 0.0.

            For properties to be plotted as x-vs-y scatter plots, the key/value pairs are:  ``"type": "scatter"``,
            ``property_names: A list of strings denoting the properties to be initialised``. All properties are
            initialized as empty lists.

            For properties that are single values, the key/value pairs are:  ``"type": "single"``,
            ``property_names: A list of strings denoting the properties to be initialised``. All properties are
            initialized with a value of 0.0.

            If not specified, uses

            .. highlight:: python
            .. code-block:: python

                default_galaxy_properties_to_analyze = {
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
                            "fraction_disk_sum", "fraction_disk_var", "SMF_history",
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
                    "single_properties": {
                        "type": "single",
                        "property_names": ["SMD_history", "SFRD_history"],
                    },
                }

        plots_that_need_smf : list of strings, optional
            The plot toggles that require the stellar mass function to be properly computed and analyzed. For example,
            plotting the quiescent fraction of galaxies requires knowledge of the total number of galaxies. The strings
            here must **EXACTLY** match the keys in ``plot_toggles``.

            If not specified, uses a default value of ``["SMF", "quiescent", "bulge_fraction", "SMF_history"]``.

        IMFs : list of strings, optional, ``{"Chabrier", "Salpeter"}``
            The initial mass functions used during the analysis of the galaxies.  This is used to shift the
            observational data points. The length of this variable is equal to the number of models to be analyzed.

            If not specified, uses a ``"Chabrier"`` IMF for each model.
        """

        num_models = len(sage_parameter_fnames)
        self._num_models = num_models

        if labels is None:
            labels = [None] * num_models

        if sage_output_formats is None:
            sage_output_formats = [None] * num_models

        if first_files_to_analyze is None:
            first_files_to_analyze = [None] * num_models

        if last_files_to_analyze is None:
            last_files_to_analyze = [None] * num_models

        if num_sage_output_files is None:
            num_sage_output_files = [None] * num_models

        if output_format_data_classes_dict is None:
            output_format_data_classes_dict = {"sage_binary": SageBinaryData, "sage_hdf5": SageHdf5Data}
        self._output_format_data_classes_dict = output_format_data_classes_dict

        if random_seeds is None:
            random_seeds = [None] * num_models

        if IMFs is None:
            IMFs = ["Chabrier"] * num_models

        # These are parameters that are model dependant.
        individual_model_parameters = [
            sage_parameter_fnames,
            sage_output_formats,
            labels,
            first_files_to_analyze,
            last_files_to_analyze,
            num_sage_output_files,
            random_seeds,
            IMFs,
        ]

        if plot_toggles is None:
            plot_toggles = default_plot_toggles
        self._plot_toggles = plot_toggles

        if history_redshifts is None:
            history_redshifts = {
                "SMF_history": "All",
                "SMD_history": "All",
                "SFRD_history": "All",
            }
        self._history_redshifts = history_redshifts

        if plots_that_need_smf is None:
            plots_that_need_smf = ["SMF", "quiescent", "bulge_fraction", "SMF_history"]

        global_model_parameters = [
            plot_toggles,
            plots_that_need_smf,
        ]

        # ``parameters`` is a matrix of parameters with each "column" specifying the parameters for a single model.
        # Hence we want to iterate through column-wise and use these to build the ``Model`` class instance. Here, the
        # ``map`` function does this transpose into a column-wise iterable.
        models = [
            Model(*model_parameters, *global_model_parameters)
            for model_parameters in map(list, zip(*individual_model_parameters))
        ]
        self._models = models

        # Determine if the stellar mass function needs to be computed for each model. Important we do this before
        # checking ``calculation_functions``.
        for model in models:
            if self._does_smf_need_computing(model):
                model.plot_toggles["SMF"] = True
                plot_toggles["SMF"] = True

        # Then populate the `calculation_methods` dictionary. This dictionary will control which properties each model
        # will calculate.  The dictionary is populated using the plot_toggles defined above.
        if calculation_functions is None:
            calculation_functions = generate_func_dict(
                plot_toggles, module_name="sage_analysis.example_calcs", function_prefix="calc_"
            )

        # Because we have adjusted the SMF, check if it's in the calculation dict.
        for model in models:
            # Condition checks for condition where SMF is being compputed but not in ``calculation_functions``.
            if model._plot_toggles.get("SMF", None) and not calculation_functions.get("SMF", None):
                raise ValueError(
                    "The stellar mass function (``SMF``) is being computed, either because it was set manually "
                    "(through ``plot_toggles``) or because another plot requires it (``plots_that_need_smf``).\n"
                    "However, ``calc_SMF`` was not found in ``calculation_functions``. Ensure that it is added."
                )

        if plot_functions is None:
            plot_functions = generate_func_dict(
                plot_toggles, module_name="sage_analysis.example_plots", function_prefix="plot_"
            )
        self._plot_functions = plot_functions

        if galaxy_properties_to_analyze is None:
            galaxy_properties_to_analyze = default_galaxy_properties_to_analyze

        for model in self._models:

            # Read the parameter files and update any of the missing parameters.
            self._read_sage_file(model)

            # Also initialise all of the properties that are required.
            for name, galaxy_properties in galaxy_properties_to_analyze.items():
                for snapshot, _ in enumerate(model._redshifts):
                    self._initialise_properties(name, model, galaxy_properties, snapshot)

            model._calculation_functions = calculation_functions

            # Go through the calculation functions and pull out those that are actually being analyzed over a number of
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

            # Determine what snapshots we need to loop over to compute the properties that are tracked over time.
            history_snaps_to_loop = self._determine_history_snapshots(model)
            model._history_snaps_to_loop = history_snaps_to_loop

        logger.info(f"Looping through snapshots {model._history_snaps_to_loop}")

    @property
    def num_models(self) -> int:
        """
        int : The number of models being analyzed.
        """
        return self._num_models

    @property
    def output_format_data_classes_dict(self) -> Dict[str, Any]:
        """
        dict [str, class] :  A dictionary that maps the output format name to the corresponding data class.
        """
        return self._output_format_data_classes_dict

    @property
    def history_redshifts(self) -> Dict[str, Union[str, List[float]]]:
        """
        dict [string, string or list of floats] : Specifies which redshifts should be analyzed for properties and
        plots that are tracked over time. The keys here **MUST** correspond to the keys in :py:attr:`~plot_toggles`. If
        the value of the entry is ``"All"``, then all snapshots will be analyzed. Otherwise, will search for the
        closest snapshots to the requested redshifts.
        """
        return self._history_redshifts

    @property
    def models(self) -> List[Model]:
        """
        list of :py:class:`~sage_analysis.model.Model` class instances : The :py:class:`~sage_analysis.model.Model` s
        being analyzed.
        """
        return self._models

    @property
    def plot_toggles(self) -> Dict[str, bool]:
        """
        dict [str, bool] : Specifies which properties should be analyzed and plotted.
        """
        return self._plot_toggles

    @property
    def plot_functions(self) -> Dict[str, Tuple[Callable, Dict[str, Any]]]:
        """
        dict [str, tuple(function, dict [str, any])] : A dictionary of functions that are used to plot the properties
        of galaxies being analyzed.  Here, the outer key is the name of the corresponding plot toggle (e.g.,
        ``"SMF"``), the value is a tuple containing the function itself (e.g., ``plot_SMF()``), and another dictionary
        which specifies any optional keyword arguments to that function with keys as the name of variable (e.g.,
        ``"plot_sub_populations"``) and values as the variable value (e.g., ``True``).

        The functions in this dictionary are called for all files analyzed and **MUST** have a signature ``func(Models,
        snapshot, plot_helper, plot_output_format, optional_keyword_arguments)``. This dict can be generated using
        :py:func:`~sage_analysis.utils.generate_func_dict`.
        """
        return self._plot_functions

    def _determine_history_snapshots(self, model: Model) -> Optional[List[int]]:
        """
        Determines which snapshots need to be iterated over to track properties over time. For each
        :py:class:`~sage_analysis.model.Model`, the ``_history_<property>_redshifts`` and
        ``_history_<property>_snapshots`` attributes are updated.

        Parameters
        ----------
        model : :py:class:`~sage_analysis.model.Model`
            The :py:class:`~sage_analysis.model.Model` instance to be updated.

        Returns
        -------
        snapshots_to_loop : list of ints
            The snapshots that need to be analyzed for this model to ensure that the requested redshifts are analyzed
            for the history properties.
        """

        # Maybe there were no history redshifts specified.
        if self._history_redshifts is None or self._history_redshifts == {}:
            return None

        # Convert these redshifts into snapshots.
        for property_name, property_redshifts in self._history_redshifts.items():

            # "All" denotes that we want all of the redshifts, otherwise use the values that were specified.
            if property_redshifts == "All":
                redshifts = model._redshifts
            else:
                redshifts = property_redshifts

            attrname = f"_history_{property_name}_redshifts"
            setattr(model, attrname, redshifts)

            # Find the snapshots that are closest to the requested redshifts.
            property_snaps = find_closest_indices(model._redshifts, redshifts)

            attrname = f"_history_{property_name}_snapshots"
            setattr(model, attrname, property_snaps)

        # Based on the snapshots required to analyze over history, we may have to loop over different snapshots.
        snapshots_to_loop: List[int] = []
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

            # Otherwise, we are following this history of this property and hence will need to loop over its snapshots.
            snaps = getattr(model, f"_history_{property_name}_snapshots")
            snapshots_to_loop.extend(snaps)

        return list(np.unique(snapshots_to_loop))

    def _read_sage_file(self, model: Model) -> None:
        """
        Reads a **SAGE** parameter file to determine all parameters such as cosmology, redshift list, etc. In
        particular, also initializes the :py:attr:`~sage_analysis.model.Model.data_class` for each model. This
        attribute is unique depending upon the value of :py:attr:`~sage_analysis.model.Model.sage_output_format` and
        the corresponding entry in :py:attr:`~output_format_data_classes_dict`.

        Parameters
        ----------
        model : :py:class:`~sage_analysis.model.Model`
            The :py:class:`~sage_analysis.model.Model` instance to be updated.
        """

        # If the format wasn't defined, then attempt to read a default parameter file to determine format.
        if model._sage_output_format is None:
            logger.info(
                f"No SAGE output format specified. Attempting to read ``{model.sage_file}`` and using the format "
                f"specified inside."
            )
            sage_dict = read_generic_sage_params(model.sage_file)

            model._sage_output_format = sage_dict["_output_format"]
            logger.info(f"Using ``{model._sage_output_format}`` output format.")

        # Each SAGE output has a specific class written to read in the data.
        model.data_class = self._output_format_data_classes_dict[model._sage_output_format](model, model.sage_file)

        # The data class has read the SAGE ini file.  Update the model with the parameters read and those specified by
        # the user. We will also log some of these.
        for key, value in model.data_class.sage_model_dict.items():

            # Check if the attribute has already been set to a non-default value.
            try:
                attr_value = getattr(model, key)
            except AttributeError:
                pass
            else:
                if attr_value is not None:
                    continue

            # At this point, the user has not specified a non-default value. Use the one read from the ini file.
            setattr(model, key, value)

            default_messages = {
                "_snapshot": f"Snapshot to analyze not specified; using final snapshot of the simulation ({value})",
                "_label": f"Label not specified; using the FileNameGalaxies from parameter file ({value})",
                "_first_file_to_analyze": f"First file to analyze not specified; using {value}",
                "_last_file_to_analyze":
                    f"Last file analyze not specified; using 1 - num cores SAGE ran with ({value})",
            }

            try:
                logger.info(default_messages[key])
            except KeyError:
                pass

        # Finally, compute the volume based on the number of files being analyzed.
        model._volume = model.data_class.determine_volume_analyzed(model)

    def _initialise_properties(
        self,
        name: str,
        model: Model,
        galaxy_properties: Dict[str, Union[str, List[str]]],
        snapshot: int,
    ) -> None:
        """
        Initialises galaxy properties that will be analyzed.

        Parameters
        ----------
        name : string
            The name of the bins if the properties will be binned or a unique identifying name otherwise.

        model : :py:class:`~sage_analysis.model.Model`
            The :py:class:`~sage_analysis.model.Model` instance to be updated.

        galaxy_properties : dict[str, float or str or list of strings]]
            The galaxy properties that will be initialized. We defer to ``galaxy_properties_to_analyze`` in the
            :py:method:`~__init__` method for a full description of this variable.

        snapshot : int
            The snapshot the properties are being updated for.
        """

        # Only currently support a few property types.
        allowed_property_types = ["binned", "scatter", "single"]
        if galaxy_properties["type"] not in allowed_property_types:
            raise ValueError(
                f"Requested to analyze a galaxy property with unkown type.  The galaxy properties were "
                f"{galaxy_properties} and the only accepted types are {allowed_property_types}."
            )

        # TODO: Should this be a dict to allow the user to specify their own property type?
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
        """
        Determines whether the stellar mass function needs to be calculated based on the values of
        :py:attr:`~plot_toggles` :py:attr:`~sage_analysis.model.Model.plots_that_need_smf`.

        Parameters
        ----------
        model : :py:class:`~sage_analysis.model.Model`
            The :py:class:`~sage_analysis.model.Model` instance we're checking.

        Returns
        -------
        bool
            A boolean indicating whether the stellar mass function needs to be computed or not.
        """

        # Maybe the SMF has already been toggled on.
        toggle = model._plot_toggles.get("SMF", None)
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

    def _determine_snapshots_to_use(
        self, snapshots: Optional[List[List[int]]], redshifts: Optional[List[List[int]]]
    ) -> List[List[int]]:
        """
        Determine which snapshots should be analyzed/plotted based on the input from the user.

        Parameters
        ---------
        snapshots : nested list of ints or string, optional
            The snapshots to analyze for each model. If both this variable and ``redshifts`` are not specified, uses
            the highest snapshot (i.e., lowest redshift) as dictated by the
            :py:attr:`~sage_analysis.model.Model.redshifts` attribute from the parameter file read for each model.

            If an entry if ``"All"``, then all snapshots for that model will be analyzed.

            The length of the outer list **MUST** be equal to :py:attr:`~num_models`.

            Warnings
            --------
            Only **ONE** of ``snapshots`` and ``redshifts`` can be specified.

        redshifts : nested list of ints, optional
            The redshift to analyze for each model. If both this variable and ``snapshots`` are not specified, uses
            the highest snapshot (i.e., lowest redshift) as dictated by the
            :py:attr:`~sage_analysis.model.Model.redshifts` attribute from the parameter file read for each model.

            The snapshots selected for analysis will be those that result in the redshifts closest to those requested.
            If an entry if ``"All"``, then all snapshots for that model will be analyzed.

            The length of the outer list **MUST** be equal to :py:attr:`~num_models`.

            Warnings
            --------
            Only **ONE** of ``snapshots`` and ``redshifts`` can be specified.

        Returns
        -------
        snapshots_for_models : nested list of ints
            The snapshots to be analyzed for each model.

        Errors
        ------
        ValueError
            Thrown if **BOTH** ``snapshots`` and ``redshifts`` are specified.
        """
        # The user cannot have non-None values for both snapshots and redshifts.
        if snapshots is not None and redshifts is not None:
            raise ValueError("Both the ``snapshots`` and ``redshift`` arguments CANNOT be non-none. Only specify one.")

        if snapshots is None and redshifts is None:
            # User hasn't explicitly specified which snapshots or redshifts they want -> use the lowest redshift ones.
            snapshots_for_models = [[len(model._redshifts) - 1] for model in self._models]
        elif snapshots == "All" or redshifts == "All":
            # User wants all snapshots (or equivalently redshifts).
            snapshots_for_models = [list(np.arange(len(model._redshifts)) - 1) for model in self._models]
        elif redshifts is not None:
            # User has specified which redshifts they want; convert to the corresponding snapshots.
            snapshots_for_models = [find_closest_indices(model._redshifts, redshifts) for model in self._models]
        elif snapshots is not None:
            # Otherwise the user has specified exactly what snapshots they want.
            snapshots_for_models = snapshots

        return snapshots_for_models

    def analyze_galaxies(
        self,
        snapshots: Optional[List[List[Union[int, str]]]] = None,
        redshifts: Optional[List[List[Union[float, str]]]] = None,
        analyze_history_snapshots: bool = True,
    ) -> None:
        """
        Analyses the galaxies of the initialized :py:attr:`~models`. These attributes will be updated directly, with
        the properties accessible via ``GalaxyAnalysis.models[<model_num>].properties[<snapshot>][<property_name>]``.

        Also, all snapshots required to track the properties over time (as specified by
        :py:attr:`~sage_analysis.model.Model._history_snaps_to_loop`) will be analyzed, unless
        ``analyze_history_snapshots`` is ``False``.

        Parameters
        ----------
        snapshots : nested list of ints or string, optional
            The snapshots to analyze for each model. If both this variable and ``redshifts`` are not specified, uses
            the highest snapshot (i.e., lowest redshift) as dictated by the
            :py:attr:`~sage_analysis.model.Model.redshifts` attribute from the parameter file read for each model.

            If an entry if ``"All"``, then all snapshots for that model will be analyzed.

            The length of the outer list **MUST** be equal to :py:attr:`~num_models`.

            Notes
            -----
            If ``analyze_history_snapshots`` is ``True``, then the snapshots iterated over will be the unique
            combination of the snapshots required for history snapshots and those specified by this variable.

            Warnings
            --------
            Only **ONE** of ``snapshots`` and ``redshifts`` can be specified.

        redshifts : nested list of ints, optional
            The redshift to analyze for each model. If both this variable and ``snapshots`` are not specified, uses
            the highest snapshot (i.e., lowest redshift) as dictated by the
            :py:attr:`~sage_analysis.model.Model.redshifts` attribute from the parameter file read for each model.

            The snapshots selected for analysis will be those that result in the redshifts closest to those requested.
            If an entry if ``"All"``, then all snapshots for that model will be analyzed.

            The length of the outer list **MUST** be equal to :py:attr:`~num_models`.

            Notes
            -----
            If ``analyze_history_snapshots`` is ``True``, then the snapshots iterated over will be the unique
            combination of the snapshots required for history snapshots and those specified by this variable.

            Warnings
            --------
            Only **ONE** of ``snapshots`` and ``redshifts`` can be specified.

        analyze_history_snapshots : bool, optional
            Specifies whether the snapshots required to analyze the properties tracked over time (e.g., stellar mass or
            star formation rate density) should be iterated over.  If not specified, then only ``snapshot`` will be
            analyzed.

        Notes
        -----
        If you wish to analyze different properties to when you initialized an instance of :py:class:`~GalaxyAnalysis`,
        you **MUST** re-initialize another instance.  Otherwise, the properties will be non-zeroed and not initialized
        correctly.

        Errors
        ------
        ValueError
            Thrown if **BOTH** ``snapshots`` and ``redshifts`` are specified.
        """

        if self._plot_toggles == {}:
            logger.debug(f"No plot toggles specified.")
            return

        baseline_snapshots_models = self._determine_snapshots_to_use(snapshots, redshifts)

        for model, baseline_snapshots in zip(self._models, baseline_snapshots_models):

            logger.info(f"Analyzing baseline snapshots {baseline_snapshots}")

            for snap in baseline_snapshots:
                # First compute all of the "normal" properties that aren't tracked over time.
                model.calc_properties_all_files(
                    model._calculation_functions, snap, debug=False, close_file=False
                )

                # Then check if this is a snapshot we're analyzing properties over time.
                if model._history_snaps_to_loop is None or not analyze_history_snapshots:
                    continue

                # Can't combine this condition with the line above because it would throw an error if ``None``.
                if snap not in model._history_snaps_to_loop:
                    continue

                model.calc_properties_all_files(
                    model._history_calculation_functions, snap, debug=False, close_file=False
                )

            # Finally, determine if there are any remaining snapshots that need to be analyzed for the history
            # properties.
            if model._history_snaps_to_loop is None or not analyze_history_snapshots:
                continue

            history_snaps = list(set(model._history_snaps_to_loop).difference(set(baseline_snapshots)))
            logger.info(f"Also analyzing snapshots {history_snaps} for the properties over redshift.")

            for snap in history_snaps:
                model.calc_properties_all_files(
                    model._history_calculation_functions, snap, debug=False, close_file=False
                )
            model.data_class.close_file(model)

    def generate_plots(
        self,
        snapshots: Optional[List[List[Union[int, str]]]] = None,
        redshifts: Optional[List[List[Union[float, str]]]] = None,
        plot_helper: Optional[PlotHelper] = None,
    ) -> Optional[List[matplotlib.figure.Figure]]:
        """
        Generates the plots for the :py:attr:`~models` being analyzed. The plots to be created are defined by the
        values of :py:attr:`~plot_toggles` specified when an instance of :py:class:`~GalaxyAnalysis` was initialized.
        If you wish to analyze different properties or create different plots, you **MUST** initialize another instance
        of :py:class:`~GalaxyAnalysis` with the new values for :py:attr:`~plot_toggles` (ensuring that values of
        ``calcuations_functions`` and ``plot_functions`` are updated if using non-default values for ``plot_toggles``).

        This method should be run after analysing the galaxies using :py:method:`~analyze_galaxies`.

        Parameters
        ----------
        snapshots : nested list of ints or string, optional
            The snapshots to plot for each model. If both this variable and ``redshifts`` are not specified, uses
            the highest snapshot (i.e., lowest redshift) as dictated by the
            :py:attr:`~sage_analysis.model.Model.redshifts` attribute from the parameter file read for each model.

            If an entry if ``"All"``, then all snapshots for that model will be analyzed.

            The length of the outer list **MUST** be equal to :py:attr:`~num_models`.

            For properties that aren't analyzed over redshift, the snapshots for each model will be plotted on each
            figure.  For example, if we are plotting a single model, setting this variable to ``[[63, 50]]`` will
            give results for snapshot 63 and 50 on each figure. For some plots (e.g., those properties that are scatter
            plotted), this is undesirable and one should instead iterate over single snapshot values instead.

            Notes
            -----
            If ``analyze_history_snapshots`` is ``True``, then the snapshots iterated over will be the unique
            combination of the snapshots required for history snapshots and those specified by this variable.

            Warnings
            --------
            Only **ONE** of ``snapshots`` and ``redshifts`` can be specified.

        redshifts : nested list of ints, optional
            The redshift to plot for each model. If both this variable and ``snapshots`` are not specified, uses
            the highest snapshot (i.e., lowest redshift) as dictated by the
            :py:attr:`~sage_analysis.model.Model.redshifts` attribute from the parameter file read for each model.

            The snapshots selected for analysis will be those that result in the redshifts closest to those requested.
            If an entry if ``"All"``, then all snapshots for that model will be analyzed.

            The length of the outer list **MUST** be equal to :py:attr:`~num_models`.

            Warnings
            --------
            Only **ONE** of ``snapshots`` and ``redshifts`` can be specified.

        plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`, optional
            A helper class that contains attributes and methods to assist with plotting. In particular, the path where
            the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
            how to initialize this class and its use.

            If not specified, then will initialize a default instance of
            :py:class:`~sage_analysis.plot_helper.PlotHelper`. Refer to the
            :py:class:`~sage_analysis.plot_helper.PlotHelper` documentation for a list of default attributes.

        Returns
        -------
        None
            Returned if :py:attr:`~plot_toggles` is an empty dictionary.

        figs
            The figures generated by the :py:attr:`~plot_functions` functions.
        """

        if self._plot_toggles == {}:
            logger.debug(f"No plot toggles specified.")
            return None

        if plot_helper is None:
            plot_helper = PlotHelper()

        snapshots = self._determine_snapshots_to_use(snapshots, redshifts)

        # Now do the plotting.
        figs: List[matplotlib.figure.Figure] = []

        for func, kwargs in self._plot_functions.values():
            fig = func(
                    self._models,
                    snapshots,
                    plot_helper,
                    **kwargs
                )

            if type(fig) == list:
                figs.extend(fig)
            else:
                figs.append(fig)

        return figs
