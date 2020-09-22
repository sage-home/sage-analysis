import logging
import os
import warnings
from pathlib import Path
from typing import Dict, List

import hypothesis.strategies as st
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from hypothesis import given, settings
from matplotlib.testing.compare import compare_images

from sage_analysis.default_analysis_arguments import (
    default_calculation_functions, default_plot_functions,
    default_plot_toggles)
from sage_analysis.galaxy_analysis import GalaxyAnalysis
from sage_analysis.model import Model
from sage_analysis.plot_helper import PlotHelper
from sage_analysis.sage_binary import SageBinaryData
from sage_analysis.sage_dark import SageDarkData
from sage_analysis.sage_dusty import SageDustyData
from sage_analysis.sage_hdf5 import SageHdf5Data
from sage_analysis.utils import generate_func_dict

logger = logging.getLogger(__name__)

test_path = Path(__file__).parent
baseline_image_path = f"{test_path}/test_data/baseline_plots/"
generated_image_path = f"{test_path}/test_data/generated_plots/"
plot_helper = PlotHelper(output_path=generated_image_path)


def calc_num_particles_in_halo(model: Model, gals, snapshot: int):

    non_zero_stellar = np.where(gals["StellarMass"][:] > 0.0)[0]
    halo_len = gals["Len"][:][non_zero_stellar]

    halos_per_bin, _ = np.histogram(halo_len, bins=model.bins["halo_len_bins"])
    model.properties[f"snapshot_{snapshot}"]["halo_len"] += halos_per_bin

def plot_num_particles_in_halo(
    models: List[Model], snapshot: int, plot_helper: PlotHelper
) -> matplotlib.figure.Figure:

    # This should make a matplotlib plot but it's adequate enough to just test that we pass through this function.
    logger.info(f"Passed through ``plot_num_particles_in_halo``.")


def my_compare_images(baseline_image_path: str, generated_image_path: str, cleanup_plots: bool = True) -> None:

    image_names = os.listdir(generated_image_path)

    for name in image_names:
        if ".DS_Store" in name:
            continue

        baseline_image = os.path.join(baseline_image_path, name)
        generated_image = os.path.join(generated_image_path, name)

        comparison = compare_images(baseline_image, generated_image, tol=1)
        assert comparison is None

        # Cleanup.
        if cleanup_plots:
            os.remove(generated_image)

@pytest.mark.parametrize("sage_output_formats", [(["sage_binary"]), (["sage_hdf5"])])
def test_sage_output_format(sage_output_formats):

    parameter_fnames = [f"{test_path}/test_data/mini-millennium.par"]
    labels = ["Mini-Millennium"]
    first_files_to_analyze = [0]
    last_files_to_analyze = [0]
    num_sage_output_files = [1]
    random_seeds = [666]

    galaxy_analysis = GalaxyAnalysis(
        parameter_fnames,
        sage_output_formats=sage_output_formats,
        num_sage_output_files=num_sage_output_files,
        first_files_to_analyze=first_files_to_analyze,
        last_files_to_analyze=last_files_to_analyze,
        labels=labels,
        random_seeds=random_seeds,
    )

    galaxy_analysis.analyze_galaxies()
    galaxy_analysis.generate_plots(plot_helper=plot_helper)

    my_compare_images(baseline_image_path, generated_image_path, False)


def test_binary_sage_num_output_file_error():
    """
    When analysing binary SAGE output, the user MUST specify the number of output files SAGE generated. Otherwise, an
    error is thrown.
    """

    parameter_fnames = [f"{test_path}/test_data/mini-millennium.par"]
    sage_output_formats = ["sage_binary"]

    with pytest.raises(RuntimeError):
        galaxy_analysis = GalaxyAnalysis(parameter_fnames, sage_output_formats=sage_output_formats)


def test_hdf5_sage_num_output_file_message(caplog):
    """
    When analysing HDF5 SAGE output, the number of output files that SAGE produced (i.e., the number of cores SAGE ran
    on) is read from the HDF5 itself.

    In the scenario where the user specifies a value for the number of output files, we generate a "DEBUG" message to
    inform that this is unecessary. Furthermore, if the passed value is differently to what is inside the HDF5 file, we
    generate a "UserWarning".
    """

    caplog.set_level(logging.DEBUG)

    parameter_fnames = [f"{test_path}/test_data/mini-millennium.par"]
    sage_output_formats = ["sage_hdf5"]
    num_sage_output_files = [100]

    with pytest.warns(UserWarning) as record:
        galaxy_analysis = GalaxyAnalysis(
            parameter_fnames,
            sage_output_formats=sage_output_formats,
            num_sage_output_files=num_sage_output_files,
        )
        galaxy_analysis.analyze_galaxies()

    # Check that the messages appeared in the log.
    assert "It is not required to specify the number of SAGE output files when analysing HDF5 output." in caplog.text
    assert record[0].message.args[0] ==  f"The number of SAGE output files according to the master HDF5 file is 1. " \
        f"However, ``analyze_sage_output`` was called with {num_sage_output_files[0]}. Using the number of files " \
        "from the HDF5 file as the correct value."


def test_use_parameter_file_format(caplog) -> None:
    """
    If no SAGE output format is specified, then it should use the output format specified in the SAGE parameter file.
    In this instance, it should log some info.
    """

    caplog.set_level(logging.INFO)

    parameter_fnames = [f"{test_path}/test_data/mini-millennium.par"]
    num_sage_output_files = [1]
    first_files_to_analyze = [0]
    last_files_to_analyze = [0]

    galaxy_analysis = GalaxyAnalysis(
        parameter_fnames,
        num_sage_output_files=num_sage_output_files,
        first_files_to_analyze=first_files_to_analyze,
        last_files_to_analyze=last_files_to_analyze,
    )
    galaxy_analysis.analyze_galaxies()

    assert galaxy_analysis._models[0].sage_output_format == "sage_binary"
    assert f"No SAGE output format specified. Attempting to read ``{parameter_fnames[0]}`` and using the " \
        "format specified inside." in caplog.text
    assert "Using ``sage_binary`` output format."


def test_multiple_models() -> None:
    """
    If the user specifies multiple models, they should be analyzed and plotted no problemo.
    """

    parameter_fnames = [f"{test_path}/test_data/mini-millennium.par", f"{test_path}/test_data/mini-millennium.par"]
    sage_output_formats = ["sage_binary", "sage_hdf5"]
    labels = ["Binary", "HDF5"]
    num_sage_output_files = [1, None]
    first_files_to_analyze = [0, None]
    last_files_to_analyze = [0, None]

    galaxy_analysis = GalaxyAnalysis(
        parameter_fnames,
        sage_output_formats=sage_output_formats,
        num_sage_output_files=num_sage_output_files,
        first_files_to_analyze=first_files_to_analyze,
        last_files_to_analyze=last_files_to_analyze,
        labels=labels,
    )
    galaxy_analysis.analyze_galaxies()

    assert len(galaxy_analysis._models) == 2

    # To ensure we have visual difference between the stellar mass function, let's add some offset to one model and
    # then plot the output.
    galaxy_analysis._models[0].properties["snapshot_63"]["SMF"] *= 1.5
    galaxy_analysis.generate_plots()

    # Save these and test equality.
# Test to ensure that wrong number of parameters yield error.

def test_defaults(caplog):
    """
    If parameters aren't specified, they should be read from the parameter file.
    """

    caplog.set_level(logging.INFO)

    parameter_fnames = [f"{test_path}/test_data/mini-millennium.par"]
    sage_output_formats = ["sage_hdf5"]
    random_seeds = [666]

    galaxy_analysis = GalaxyAnalysis(
        parameter_fnames,
        sage_output_formats=sage_output_formats,
        random_seeds=random_seeds,
    )

    assert galaxy_analysis._models[0].snapshot == 63
    assert "Snapshot to analyze not specified; using final snapshot of the simulation (63)" in caplog.text

    assert galaxy_analysis._models[0].label == "correct-mini-millennium-output"
    assert "Label not specified; using the FileNameGalaxies from parameter file (correct-mini-millennium-output)" \
        in caplog.text

    assert galaxy_analysis._models[0].first_file_to_analyze == 0
    assert "First file to analyze not specified; using 0" in caplog.text

    assert galaxy_analysis._models[0].last_file_to_analyze == 0
    assert "Last file to analyze not specified; using num cores SAGE ran with minus 1 (0)"

    # The only thing that should be different to the baseline plots is the label in the legend. Hence if we update the
    # label and generate the plots, they should be identical.
    galaxy_analysis.analyze_galaxies()
    galaxy_analysis._models[0]._label = "Mini-Millennium"

    generated_image_path = "test_data/generated_plots/"
    galaxy_analysis.generate_plots(plot_helper=plot_helper)

    my_compare_images(baseline_image_path, generated_image_path)

def test_additional_property(caplog):
    """
    A new property can be analyzed and accounted for.
    """

    caplog.set_level(logging.INFO)

    parameter_fnames = [f"{test_path}/test_data/mini-millennium.par"]
    sage_output_formats = ["sage_hdf5"]

    baseline_plot_toggles = {"SMF": True}
    baseline_calculation_functions = generate_func_dict(baseline_plot_toggles, "sage_analysis.example_calcs", "calc_")
    baseline_plot_functions = generate_func_dict(baseline_plot_toggles, "sage_analysis.example_plots", "plot_")

    extra_plot_toggles = {"num_particles_in_halo": True}
    extra_calculation_functions = generate_func_dict(extra_plot_toggles, "sage_analysis.tests.test_galaxy_analysis", "calc_")
    extra_plot_functions = generate_func_dict(extra_plot_toggles, "sage_analysis.tests.test_galaxy_analysis", "plot_")

    baseline_plot_toggles.update(extra_plot_toggles)
    baseline_calculation_functions.update(extra_calculation_functions)
    baseline_plot_functions.update(extra_plot_functions)

    galaxy_properties_to_analyze = {
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

    galaxy_analysis = GalaxyAnalysis(
        parameter_fnames,
        sage_output_formats=sage_output_formats,
        plot_toggles=baseline_plot_toggles,
        calculation_functions=baseline_calculation_functions,
        plot_functions=baseline_plot_functions,
        galaxy_properties_to_analyze=galaxy_properties_to_analyze,
    )

    galaxy_analysis.analyze_galaxies()
    assert f"Initialized galaxy properties {galaxy_properties_to_analyze['stellar_mass_bins']} for Snapshot 63" in caplog.text
    assert f"Initialized galaxy properties {galaxy_properties_to_analyze['halo_len_bins']} for Snapshot 63" in caplog.text
    galaxy_analysis.generate_plots(plot_helper=plot_helper)
    assert "Passed through ``plot_num_particles_in_halo``." in caplog.text
    os.remove(f"{generated_image_path}/1.StellarMassFunction.png")

def test_incorrect_additional_property():
    """
    The only accepted types of the properties are one of ``["binned", "scatter", "single"]``. Otherwise, it should
    throw a ``ValueError``.
    """

    parameter_fnames = [f"{test_path}/test_data/mini-millennium.par"]
    sage_output_formats = ["sage_hdf5"]

    galaxy_properties_to_analyze = {
        "incorrect_property": {
                    "type": "this is an incorrect property type",
        },
    }

    with pytest.raises(ValueError):
        galaxy_analysis = GalaxyAnalysis(
            parameter_fnames,
            sage_output_formats=sage_output_formats,
            galaxy_properties_to_analyze=galaxy_properties_to_analyze,
        )

@pytest.mark.hypothesis
@given(st.dictionaries(st.sampled_from(list(default_plot_toggles.keys())), st.sampled_from([True, False])))
@settings(deadline=None)
def test_random_combinations(plot_toggles: Dict[str, bool]) -> None:
    """
    Ensure that any random combination of plot toggles run as expected using default arugments.

    Note: This runs 100 different combinations of plot toggles, each time opening/closing files and generating plots.
    Consequently, it takes ~2 minutes to run. This test can be skipped by running ``pytest -m "not hypothesis"``.
    """

    parameter_fnames = [f"{test_path}/test_data/mini-millennium.par"]
    sage_output_formats = ["sage_hdf5"]
    labels = ["Mini-Millennium"]
    random_seeds = [666]
    generated_image_path = "test_data/generated_plots/"

    galaxy_analysis = GalaxyAnalysis(
        parameter_fnames,
        sage_output_formats=sage_output_formats,
        random_seeds=random_seeds,
        labels=labels,
        plot_toggles=plot_toggles,
    )

    galaxy_analysis.analyze_galaxies()
    figs = galaxy_analysis.generate_plots(plot_helper=plot_helper)

    # When ``plot_toggles = {}``, no figures are generated. This is still a valid scenario to test however there isn't
    # anything we want to be comparing.
    if figs is None:
        return

    my_compare_images(baseline_image_path, generated_image_path)


def test_multiple_snapshots():
    """
    The user can request to analyse + plot multiple snapshots at the same time.
    """

    parameter_fnames = [f"{test_path}/test_data/mini-millennium.par"]
    sage_output_formats = ["sage_hdf5"]
    labels = ["Mini-Millennium"]
    random_seeds = [666]
    generated_image_path = "test_data/generated_plots/"

    galaxy_analysis = GalaxyAnalysis(
        parameter_fnames,
        sage_output_formats=sage_output_formats,
        random_seeds=random_seeds,
        labels=labels,
    )

    snapshots = [[63, 37, 32]]
    galaxy_analysis.analyze_galaxies(snapshots=snapshots)

    # galaxy_analysis.generate_plots(plot_helper=plot_helper, snapshots=snapshots)

    # Check these plots to ensure they're all good.


def test_keyword_arguments():

    parameter_fnames = [f"{test_path}/test_data/mini-millennium.par"]
    sage_output_formats = ["sage_hdf5"]
    labels = ["Mini-Millennium"]
    random_seeds = [666]
    generated_image_path = "test_data/generated_plots/keyword_args"

    plot_toggles = {"SMF": True, "quiescent": True, "baryon_fraction": True}

    kwargs = {"SMF": {"calc_sub_populations": True}}
    calculation_functions = generate_func_dict(plot_toggles, "sage_analysis.example_calcs", "calc_", kwargs)

    kwargs = {
        "SMF": {"plot_sub_populations": True},
        "quiescent": {"plot_sub_populations": True},
        "baryon_fraction": {"plot_sub_populations": True}
    }
    plot_functions = generate_func_dict(plot_toggles, "sage_analysis.example_plots", "plot_", kwargs)


    galaxy_analysis = GalaxyAnalysis(
        parameter_fnames,
        sage_output_formats=sage_output_formats,
        random_seeds=random_seeds,
        labels=labels,
        plot_toggles=plot_toggles,
        calculation_functions=calculation_functions,
        plot_functions=plot_functions
    )

    galaxy_analysis.analyze_galaxies()
    # galaxy_analysis.generate_plots(plot_helper=plot_helper)


def test_both_snapshots_and_redshift():
    """
    User requests to read both snapshots and redshifts. Should throw an error.
    """

    parameter_fnames = [f"{test_path}/test_data/mini-millennium.par"]
    sage_output_formats = ["sage_hdf5"]

    galaxy_analysis = GalaxyAnalysis(
        parameter_fnames,
        sage_output_formats=sage_output_formats,
    )

    snapshots = "All"
    redshifts = "All"

    with pytest.raises(ValueError):
        galaxy_analysis.analyze_galaxies(snapshots=snapshots, redshifts=redshifts)

    with pytest.raises(ValueError):
        galaxy_analysis.generate_plots(snapshots=snapshots, redshifts=redshifts)

def test_highest_snapshot():
    """
    User explicitly asked for Snapshot 63 (the highest snapshot).  This should correspond to the default scenario.
    """

    parameter_fnames = [f"{test_path}/test_data/mini-millennium.par"]
    sage_output_formats = ["sage_hdf5"]
    labels = ["Mini-Millennium"]
    random_seeds = [666]
    generated_image_path = "test_data/generated_plots/"

    galaxy_analysis = GalaxyAnalysis(
        parameter_fnames,
        sage_output_formats=sage_output_formats,
        random_seeds=random_seeds,
        labels=labels,
    )

    snapshots = [[63]]
    galaxy_analysis.analyze_galaxies(snapshots=snapshots)

    galaxy_analysis.generate_plots(plot_helper=plot_helper, snapshots=snapshots)

    my_compare_images(baseline_image_path, generated_image_path)



def test_redshift_zero():
    """
    User explicitly asked for redshift 0.  This should correspond to the default scenario.
    """

    parameter_fnames = [f"{test_path}/test_data/mini-millennium.par"]
    sage_output_formats = ["sage_hdf5"]
    labels = ["Mini-Millennium"]
    random_seeds = [666]
    generated_image_path = "test_data/generated_plots/"

    galaxy_analysis = GalaxyAnalysis(
        parameter_fnames,
        sage_output_formats=sage_output_formats,
        random_seeds=random_seeds,
        labels=labels,
    )

    redshifts = [[0.0]]
    galaxy_analysis.analyze_galaxies(redshifts=redshifts)

    galaxy_analysis.generate_plots(plot_helper=plot_helper, redshifts=redshifts)

    my_compare_images(baseline_image_path, generated_image_path)


def test_other_sages():

    parameter_fnames = [
        f"{test_path}/test_data/mini-millennium.par",
        "/Users/jacobseiler/Desktop/DarkSage/input/millennium.par",
        "/Users/jacobseiler/Desktop/dusty-sage/src/auxdata/trees/mini-millennium/mini-millennium.par",
    ]
    sage_output_formats = ["sage_hdf5", "sage_dark", "sage_dusty"]
    labels = ["Base SAGE", "Dark SAGE", "Dusty SAGE"]
    random_seeds = [666, 666, 666]
    generated_image_path = "test_data/generated_plots/other_sages/"
    output_format_data_classes_dict = {
        "sage_binary": SageBinaryData,
        "sage_hdf5": SageHdf5Data,
        "sage_dark": SageDarkData,
        "sage_dusty": SageDustyData,
    }
    plot_toggles = {"SMF": True}

    first_files_to_analyze = [0, 0, 0]
    last_files_to_analyze = [0, 7, 7]
    num_sage_output_files = [1, 8, 8]

    galaxy_analysis = GalaxyAnalysis(
        parameter_fnames,
        sage_output_formats=sage_output_formats,
        random_seeds=random_seeds,
        plot_toggles=plot_toggles,
        labels=labels,
        output_format_data_classes_dict=output_format_data_classes_dict,
        num_sage_output_files=num_sage_output_files,
        first_files_to_analyze=first_files_to_analyze,
        last_files_to_analyze=last_files_to_analyze,
    )

    galaxy_analysis.analyze_galaxies()

    galaxy_analysis.generate_plots(plot_helper=plot_helper)

    # my_compare_images(baseline_image_path, generated_image_path)
