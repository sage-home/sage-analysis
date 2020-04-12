from typing import List
import warnings
import logging

import numpy as np
import pytest
import matplotlib.pyplot as plt
import os

from matplotlib.testing.compare import compare_images

from sage_analysis.utils import generate_func_dict
from sage_analysis.default_analysis_arguments import default_plot_toggles
from sage_analysis.galaxy_analysis import GalaxyAnalysis
from sage_analysis.model import Model
from sage_analysis.tests.test_galaxy_analysis import my_compare_images

logger = logging.getLogger(__name__)
baseline_image_path = "test_data/baseline_plots/"
generated_image_path = "test_data/generated_plots/"

@pytest.mark.parametrize("sage_output_formats", [(["sage_binary"]), (["sage_hdf5"])])
def test_history(sage_output_formats: List[str]) -> None:
    """
    Ensure generated history plots are identical for both output formats.
    """

    parameter_fnames = ["test_data/mini-millennium.par"]
    labels = ["Mini-Millennium"]
    first_files_to_analyse = [0]
    last_files_to_analyse = [0]
    num_sage_output_files = [1]
    plot_toggles = {"SMF_history": True, "SFRD_history": True, "SMD_history": True}

    galaxy_analysis = GalaxyAnalysis(
        parameter_fnames,
        sage_output_formats=sage_output_formats,
        num_sage_output_files=num_sage_output_files,
        first_files_to_analyse=first_files_to_analyse,
        last_files_to_analyse=last_files_to_analyse,
        plot_toggles=plot_toggles,
        plot_output_path=generated_image_path,
        labels=labels,
    )

    galaxy_analysis.analyse_galaxies()
    galaxy_analysis.generate_plots()

    my_compare_images(baseline_image_path, generated_image_path)


@pytest.mark.parametrize("sage_output_formats", [(["sage_binary"]), (["sage_hdf5"])])
def test_history_and_baseline(sage_output_formats: List[str]) -> None:
    """
    Ensure that if ALL plot toggles are turned on, then all plots can be generated in one go.
    """

    parameter_fnames = ["test_data/mini-millennium.par"]
    labels = ["Mini-Millennium"]
    first_files_to_analyse = [0]
    last_files_to_analyse = [0]
    num_sage_output_files = [1]
    random_seeds = [666]

    plot_toggles = default_plot_toggles
    history_properties = ["SMF_history", "SFRD_history", "SMD_history"]

    for key in history_properties:
        plot_toggles[key] = True

    galaxy_analysis = GalaxyAnalysis(
        parameter_fnames,
        sage_output_formats=sage_output_formats,
        num_sage_output_files=num_sage_output_files,
        first_files_to_analyse=first_files_to_analyse,
        last_files_to_analyse=last_files_to_analyse,
        random_seeds=random_seeds,
        plot_toggles=plot_toggles,
        plot_output_path=generated_image_path,
        labels=labels,
    )

    galaxy_analysis.analyse_galaxies()
    galaxy_analysis.generate_plots()

    my_compare_images(baseline_image_path, generated_image_path)
