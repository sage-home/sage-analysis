import os
from typing import List, Dict, Union, Optional

import sage_analysis.example_calcs
import sage_analysis.example_plots

from sage_analysis.utils import generate_func_dict
from sage_analysis.model import Model
from sage_analysis.sage_binary import SageBinaryData
try:
    from sage_analysis.sage_hdf5 import SageHdf5Data
except ImportError:
    print("h5py not found.  If you're reading in HDF5 output from SAGE, please install this package.")

import numpy as np

import matplotlib

# Sometimes we divide a galaxy that has zero mass (e.g., no cold gas). Ignore these
# warnings as they spam stdout. Also remember the old settings.
old_error_settings = np.seterr()
np.seterr(all="ignore")


def analyse_sage_output(
    par_fnames: List[str],
    snapshots_to_plot: Optional[List[int]] = None,
    IMFs: Optional[List[str]] = None,
    labels: Optional[List[str]] = None,
    sage_output_formats: Optional[List[str]] = None,
    first_files_to_plot: Optional[List[int]] = None,
    last_files_to_plot: Optional[List[int]] = None,
    num_output_files_to_plot: Optional[List[int]] = None,
    plot_toggles: Optional[Dict[str, int]] = None,
    plot_output_format: str = "png",
    plot_output_path: str = "./plots_new",
) -> List[matplotlib.figure.Figure]:

    if snapshots_to_plot is None:
        snapshots_to_plot = [63]

    if IMFs is None:
        IMFs = ["Chabrier"]

    if labels is None:
        labels = ["Mini-Millennium"]

    if sage_output_formats is None:
        sage_output_formats = ["sage_binary"]

    if first_files_to_plot is None:
        first_files_to_plot = [0]

    if last_files_to_plot is None:
        last_files_to_plot = [0]

    if num_output_files_to_plot is None:
        num_output_files_to_plot = [1]

    parameters = [
        par_fnames,
        snapshots_to_plot,
        IMFs,
        labels,
        sage_output_formats,
        first_files_to_plot,
        last_files_to_plot,
        num_output_files_to_plot
    ]

    # ``parameters`` is a matrix of parameters with each "column" specifying the parameters for a single model. Hence
    # we want to iteratre through column-wise and use these to build the ``Model`` class instance.
    parameters = np.array(parameters)
    models_to_plot = [
        Model(*model_parameters) for model_parameters in parameters.T
    ]

    if plot_toggles is None:
        plot_toggles = {
            "SMF" : 1,  # Stellar mass function.
            "BMF" : 1,  # Baryonic mass function.
            "GMF" : 1,  # Gas mass function (cold gas).
            "BTF" : 1,  # Baryonic Tully-Fisher.
            "sSFR" : 1,  # Specific star formation rate.
            "gas_fraction" : 1,  # Fraction of galaxy that is cold gas.
            "metallicity" : 1,  # Metallicity scatter plot.
            "bh_bulge" : 1,  # Black hole-bulge relationship.
            "quiescent" : 1,  # Fraction of galaxies that are quiescent.
            "bulge_fraction" : 1,  # Fraction of galaxies that are bulge/disc dominated.
            "baryon_fraction" : 1,  # Fraction of baryons in galaxy/reservoir.
            "reservoirs" : 1,  # Mass in each reservoir.
            "spatial" : 1   # Spatial distribution of galaxies.
        }


    # Generate directory for output plots.
    if not os.path.exists(plot_output_path):
        os.makedirs(plot_output_path)

    # Go through each model and calculate all the required properties.
    for model in models_to_plot:

        model.plot_output_format = plot_output_format

        # Each SAGE output has a specific class written to read in the data.
        if model.sage_output_format == "sage_binary":
            model.data_class = SageBinaryData(
                model, model.num_output_files, model.sage_file, model.snapshot
            )
        elif model.sage_output_format == "sage_hdf5":
            model.data_class = SageHdf5Data(model, model.sage_file)

        # The data class has read the SAGE ini file.  Update the model with the parameters
        # read and those specified by the user.
        model.update_attributes({}, plot_toggles)

        #my_model.data_class.set_cosmology(my_model)

        # Some properties require the stellar mass function to normalize their values. For
        # these, the SMF plot toggle is explicitly required.
        try:
            if plot_toggles["SMF"]:
                model.calc_SMF = True
            else:
                model.calc_SMF = False
        except KeyError:  # Maybe we've removed "SMF" from plot_toggles...
                model.calc_SMF = False

        # Then populate the `calculation_methods` dictionary. This dictionary will control
        # which properties each model will calculate.  The dictionary is populated using
        # the plot_toggles defined above.
        # Our functions are inside the `example_calcs.py` module and are named "calc_<toggle>". If
        # your functions are in a different module or different function prefix, change it
        # here.
        # ALL FUNCTIONS MUST HAVE A FUNCTION SIGNATURE `func(Model, gals, optional_kwargs=...)`.
        calculation_functions = generate_func_dict(plot_toggles, module_name="sage_analysis.example_calcs",
                                                   function_prefix="calc_")

        # Finally, before we calculate the properties, we need to decide how each property
        # is stored. Properties can be binned (e.g., how many galaxies with mass between 10^8.0
        # and 10^8.1), scatter plotted (e.g., for 1000 galaxies plot the specific star
        # formation rate versus stellar mass) or a single number (e.g., the sum
        # of the star formation rate at a snapshot). Properties can be accessed using
        # `Model.properties["property_name"]`; e.g., `Model.properties["SMF"]`.

        # First let's do the properties binned on stellar mass. The bins themselves can be
        # accessed using `Model.bins["bin_name"]`; e.g., `Model.bins["stellar_mass_bins"]
        stellar_properties = ["SMF", "red_SMF", "blue_SMF", "BMF", "GMF",
                              "centrals_MF", "satellites_MF", "quiescent_galaxy_counts",
                              "quiescent_centrals_counts", "quiescent_satellites_counts",
                              "fraction_bulge_sum", "fraction_bulge_var",
                              "fraction_disk_sum", "fraction_disk_var"]
        model.init_binned_properties(8.0, 12.0, 0.1, "stellar_mass_bins",
                                        stellar_properties)

        # Properties binned on halo mass.
        halo_properties = ["fof_HMF"]
        component_properties = ["halo_{0}_fraction_sum".format(component) for component in
                               ["baryon", "stars", "cold", "hot", "ejected", "ICS", "bh"]]
        model.init_binned_properties(10.0, 14.0, 0.1, "halo_mass_bins",
                                        halo_properties+component_properties)

        # Now properties that will be extended as lists.
        scatter_properties = ["BTF_mass", "BTF_vel", "sSFR_mass", "sSFR_sSFR",
                              "gas_frac_mass", "gas_frac", "metallicity_mass",
                              "metallicity", "bh_mass", "bulge_mass", "reservoir_mvir",
                              "reservoir_stars", "reservoir_cold", "reservoir_hot",
                              "reservoir_ejected", "reservoir_ICS", "x_pos",
                              "y_pos", "z_pos"]
        model.init_scatter_properties(scatter_properties)

        # Finally those properties that are stored as a single number.
        single_properties = []
        model.init_single_properties(single_properties)

        # To be more memory concious, we calculate the required properties on a
        # file-by-file basis. This ensures we do not keep ALL the galaxy data in memory.
        model.calc_properties_all_files(calculation_functions, debug=False)

    # Similar to the calculation functions, all of the plotting functions are in the
    # `plots.py` module and are labelled `plot_<toggle>`.
    plot_functions = generate_func_dict(plot_toggles, module_name="sage_analysis.example_plots",
                                        function_prefix="plot_")

    # Now do the plotting.
    figs = []
    for func_name in plot_functions.keys():
        func = plot_functions[func_name][0]
        keyword_args = plot_functions[func_name][1]

        fig = func(models_to_plot, plot_output_path, plot_output_format, **keyword_args)

        if type(fig) == list:
            figs.append(*fig)
        else:
            figs.append(fig)

    # Set the error settings to the previous ones so we don't annoy the user.
    np.seterr(divide=old_error_settings["divide"], over=old_error_settings["over"],
              under=old_error_settings["under"], invalid=old_error_settings["invalid"])

    return figs
