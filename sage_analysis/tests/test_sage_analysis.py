import pytest
import matplotlib.pyplot as plt
import os

from matplotlib.testing.compare import compare_images
from sage_analysis.galaxy_analysis import analyse_sage_output

"""
def test_binary():

    fname = ["test_data/mini-millennium.par"]
    analyse_sage_output(fname)
"""


def test_sage_output_format():

    # First ensure that we actually have baseline plots to compare against.
    baseline_image_path = "test_data/baseline_plots/"
    baseline_image_names = [os.path.join(baseline_image_path, name) for name in os.listdir(baseline_image_path)]

    if len(baseline_image_names) == 0:
        raise ValueError(f"No baseline images found.")

    # Now generate the plots.
    fname = ["test_data/mini-millennium.par"]
    sage_output_formats = ["sage_hdf5"]
    generated_image_path = "test_data/generated_plots/"

    analyse_sage_output(fname, sage_output_formats=sage_output_formats, plot_output_path=generated_image_path)

    # Finally, compare and ensure that the images are as expected.
    generated_image_names = [os.path.join(generated_image_path, name) for name in os.listdir(generated_image_path)]
    for baseline, generated in zip(baseline_image_names, generated_image_names):
        compare_images(baseline, generated, tol=10)

        # Cleanup.
        os.remove(generated)
