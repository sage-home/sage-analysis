from typing import List
import warnings
import logging

import numpy as np
import pytest
import matplotlib.pyplot as plt
import os

from matplotlib.testing.compare import compare_images

from sage_analysis.utils import generate_func_dict
from sage_analysis.galaxy_analysis import analyse_sage_output, _generate_plots
from sage_analysis.model import Model

logger = logging.getLogger(__name__)

def calc_num_particles_in_halo(model: Model, gals):

    non_zero_stellar = np.where(gals["StellarMass"][:] > 0.0)[0]
    halo_len = gals["Len"][:][non_zero_stellar]

    halos_per_bin, _ = np.histogram(halo_len, bins=model.bins["halo_len_bins"])
    model.properties["halo_len"] += halos_per_bin

def plot_num_particles_in_halo(models: List[Model], plot_output_path: str, plot_output_format: str = "png"):

    # This should make a matplotlib plot but it's adequate enough to just test that we pass through this function.
    logger.debug(models[0].properties["halo_len"])

def my_compare_images(baseline_image_path: str, generated_image_path: str, cleanup_plots: bool = True) -> None:

    # First ensure that we actually have baseline plots to compare against.
    baseline_image_names = [os.path.join(baseline_image_path, name) for name in os.listdir(baseline_image_path)]
    baseline_image_names = [image for image in baseline_image_names if ".DS_Store" not in image]

    if len(baseline_image_names) == 0:
        raise ValueError(f"No baseline images found.")

    # Finally, compare and ensure that the images are as expected.
    generated_image_names = [os.path.join(generated_image_path, name) for name in os.listdir(generated_image_path)]
    generated_image_names = [image for image in generated_image_names if ".DS_Store" not in image]

    for baseline, generated in zip(baseline_image_names, generated_image_names):
        comparison = compare_images(baseline, generated, tol=1)
        assert comparison is None

        # Cleanup.
        if cleanup_plots:
            os.remove(generated)

@pytest.mark.parametrize("sage_output_formats", [(["sage_binary"]), (["sage_hdf5"])])
def test_sage_output_format(sage_output_formats):

    fname = ["test_data/mini-millennium.par"]
    labels = ["Mini-Millennium"]
    first_files_to_analyse = [0]
    last_files_to_analyse = [0]
    num_sage_output_files = [1]
    generated_image_path = "test_data/generated_plots/"
    random_seeds = [666]

    analyse_sage_output(
        fname,
        sage_output_formats=sage_output_formats,
        num_sage_output_files=num_sage_output_files,
        first_files_to_analyse=first_files_to_analyse,
        last_files_to_analyse=last_files_to_analyse,
        plot_output_path=generated_image_path,
        labels=labels,
        random_seeds=random_seeds,
    )

    baseline_image_path = "test_data/baseline_plots/"
    my_compare_images(baseline_image_path, generated_image_path)


def test_binary_sage_num_output_file_error():
    """
    When analysing binary SAGE output, the user MUST specify the number of output files SAGE generated. Otherwise, an
    error is thrown.
    """

    fname = ["test_data/mini-millennium.par"]
    sage_output_formats = ["sage_binary"]

    with pytest.raises(RuntimeError):
        analyse_sage_output(fname, sage_output_formats=sage_output_formats)


def test_hdf5_sage_num_output_file_message(caplog):
    """
    When analysing HDF5 SAGE output, the number of output files that SAGE produced (i.e., the number of cores SAGE ran
    on) is read from the HDF5 itself.

    In the scenario where the user specifies a value for the number of output files, we generate a "DEBUG" message to
    inform that this is unecessary. Furthermore, if the passed value is differently to what is inside the HDF5 file, we
    generate a "UserWarning".
    """

    caplog.set_level(logging.DEBUG)

    fname = ["test_data/mini-millennium.par"]
    sage_output_formats = ["sage_hdf5"]
    num_sage_output_files = [100]

    with pytest.warns(UserWarning) as record:
        analyse_sage_output(
            fname,
            sage_output_formats=sage_output_formats,
            num_sage_output_files=num_sage_output_files,
            generate_plots=False,
        )

    # Check that the messages appeared in the log.
    assert "It is not required to specify the number of SAGE output files when analysing HDF5 output." in caplog.text
    assert record[0].message.args[0] ==  f"The number of SAGE output files according to the master HDF5 file is 1. " \
        f"However, ``analyse_sage_output`` was called with {num_sage_output_files[0]}. Using the number of files " \
        "from the HDF5 file as the correct value."


def test_use_parameter_file_format(caplog) -> None:
    """
    If no SAGE output format is specified, then it should use the output format specified in the SAGE parameter file.
    In this instance, it should log some info.
    """

    caplog.set_level(logging.INFO)

    fname = ["test_data/mini-millennium.par"]
    num_sage_output_files = [1]
    first_files_to_analyse = [0]
    last_files_to_analyse = [0]

    models = analyse_sage_output(
        fname,
        num_sage_output_files=num_sage_output_files,
        first_files_to_analyse=first_files_to_analyse,
        last_files_to_analyse=last_files_to_analyse,
        generate_plots=False,
    )

    assert models[0].sage_output_format == "sage_binary"
    assert "No SAGE output format specified. Attempting to read ``test_data/mini-millennium.par`` and using the " \
        "format specified inside." in caplog.text
    assert "Using ``sage_binary`` output format."


def test_multiple_models() -> None:
    """
    If the user specifies multiple models, they should be analysed and plotted no problemo.
    """

    fname = ["test_data/mini-millennium.par", "test_data/mini-millennium.par"]
    sage_output_formats = ["sage_binary", "sage_hdf5"]
    labels = ["Binary", "HDF5"]
    num_sage_output_files = [1, None]
    first_files_to_analyse = [0, None]
    last_files_to_analyse = [0, None]

    models = analyse_sage_output(
        fname,
        sage_output_formats=sage_output_formats,
        num_sage_output_files=num_sage_output_files,
        first_files_to_analyse=first_files_to_analyse,
        last_files_to_analyse=last_files_to_analyse,
        labels=labels,
        generate_plots=False,
    )

    assert len(models) == 2

    # To ensure we have visual difference between the stellar mass function, let's add some offset to one model and
    # then plot the output.
    models[0].properties["SMF"] *= 1.5

    generated_image_path = "test_data/generated_plots/"
    plot_functions = generate_func_dict(
        models[0].plot_toggles, module_name="sage_analysis.example_plots", function_prefix="plot_"
    )
    _generate_plots(models, generated_image_path, "png", plot_functions)

    # Save these and test equality.

# Need to add an option to add a tag to outputs?
# Test to ensure that wrong number of parameters yield error.


def test_defaults(caplog):
    """
    If parameters aren't specified, they should be read from the parameter file.
    """

    caplog.set_level(logging.INFO)

    fname = ["test_data/mini-millennium.par"]
    sage_output_formats = ["sage_hdf5"]
    random_seeds = [666]

    models = analyse_sage_output(
        fname,
        sage_output_formats=sage_output_formats,
        random_seeds=random_seeds,
        generate_plots=False,
    )

    assert models[0].snapshot == 63
    assert "Snapshot to analyse not specified; using the final snapshot of the simulation (63)" in caplog.text

    assert models[0].label == "correct-mini-millennium-output"
    assert "Label not specified; using the FileNameGalaxies from parameter file (correct-mini-millennium-output)" \
        in caplog.text

    assert models[0].first_file_to_analyse == 0
    assert "First file to analyse not specified; using 0" in caplog.text

    assert models[0].last_file_to_analyse == 0
    assert "Last file to analyse not specified; using num cores SAGE ran with minus 1 (0)"

    # The only thing that should be different to the baseline plots is the label in the legend. Hence if we update the
    # label and generate the plots, they should be identical.
    models[0].label = "Mini-Millennium"

    generated_image_path = "test_data/generated_plots/"
    plot_functions = generate_func_dict(
        models[0].plot_toggles, module_name="sage_analysis.example_plots", function_prefix="plot_"
    )
    _generate_plots(models, generated_image_path, "png", plot_functions)

    baseline_image_path = "test_data/baseline_plots/"
    my_compare_images(baseline_image_path, generated_image_path)


def test_additional_property(caplog):
    """
    A new property can be analysed and accounted for.
    """

    caplog.set_level(logging.INFO)

    fname = ["test_data/mini-millennium.par"]
    sage_output_formats = ["sage_hdf5"]

    baseline_plot_toggles = {"SMF": True}
    baseline_calculation_functions = generate_func_dict(baseline_plot_toggles, "sage_analysis.example_calcs", "calc_")
    # baseline_plot_functions = generate_func_dict(baseline_plot_toggles, "sage_analysis.example_plots", "plot_")

    extra_plot_toggles = {"num_particles_in_halo": True}
    extra_calculation_functions = generate_func_dict(extra_plot_toggles, "sage_analysis.tests.test_sage_analysis", "calc_")
    # extra_plot_functions = generate_func_dict(extra_plot_toggles, "sage_analysis.tests.test_sage_analysis", "plot_")

    baseline_plot_toggles.update(extra_plot_toggles)
    baseline_calculation_functions.update(extra_calculation_functions)
    # baseline_plot_functions.update(extra_plot_functions)

    galaxy_properties_to_analyse = {
        "stellar_mass_bins": {
                    "type": "binned",
                    "bin_low": 8.0,
                    "bin_high": 12.0,
                    "bin_width": 0.1,
                    "property_names": ["SMF", "red_SMF", "blue_SMF"],
        },
        "halo_len_bins": {
            "type": "binned",
            "bin_low": 1.0,
            "bin_high": 1000.0,
            "bin_width": 1,
            "property_names": ["halo_len"],
        }
    }

    models = analyse_sage_output(
        fname,
        sage_output_formats=sage_output_formats,
        plot_toggles=baseline_plot_toggles,
        calculation_functions=baseline_calculation_functions,
        # plot_functions=baseline_plot_functions,
        galaxy_properties_to_analyse=galaxy_properties_to_analyse,
        generate_plots=False,
    )

    assert f"Initialized galaxy properties {galaxy_properties_to_analyse['stellar_mass_bins']}" in caplog.text
    assert f"Initialized galaxy properties {galaxy_properties_to_analyse['halo_len_bins']}" in caplog.text


def test_incorrect_additional_property():
    """
    The only accepted types of the properties are one of ``["binned", "scatter", "single"]``. Otherwise, it should
    throw a ``ValueError``.
    """

    fname = ["test_data/mini-millennium.par"]
    sage_output_formats = ["sage_hdf5"]

    galaxy_properties_to_analyse = {
        "incorrect_property": {
                    "type": "this is an incorrect property type",
        },
    }

    with pytest.raises(ValueError):
        models = analyse_sage_output(
            fname,
            sage_output_formats=sage_output_formats,
            galaxy_properties_to_analyse=galaxy_properties_to_analyse,
        )
