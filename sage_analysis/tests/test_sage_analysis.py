import warnings
import logging

import pytest
import matplotlib.pyplot as plt
import os

from matplotlib.testing.compare import compare_images
from sage_analysis.galaxy_analysis import analyse_sage_output


@pytest.mark.parametrize("sage_output_formats", [(["sage_binary"]), (["sage_hdf5"])])
def test_sage_output_format(sage_output_formats):

    # First ensure that we actually have baseline plots to compare against.
    baseline_image_path = "test_data/baseline_plots/"
    baseline_image_names = [os.path.join(baseline_image_path, name) for name in os.listdir(baseline_image_path)]

    if len(baseline_image_names) == 0:
        raise ValueError(f"No baseline images found.")

    # Now generate the plots.
    fname = ["test_data/mini-millennium.par"]
    num_sage_output_files = [1]
    generated_image_path = "test_data/generated_plots/"

    analyse_sage_output(
        fname,
        sage_output_formats=sage_output_formats,
        num_sage_output_files=num_sage_output_files,
        plot_output_path=generated_image_path
    )

    # Finally, compare and ensure that the images are as expected.
    generated_image_names = [os.path.join(generated_image_path, name) for name in os.listdir(generated_image_path)]
    for baseline, generated in zip(baseline_image_names, generated_image_names):
        compare_images(baseline, generated, tol=10)

        # Cleanup.
        os.remove(generated)



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
    generated_image_path = "test_data/generated_plots/"

    with pytest.warns(UserWarning) as record:
        analyse_sage_output(
            fname,
            sage_output_formats=sage_output_formats,
            num_sage_output_files=num_sage_output_files,
            plot_output_path=generated_image_path
        )

    # Check that the messages appeared in the log.
    assert "It is not required to specify the number of SAGE output files when analysing HDF5 output." in caplog.text
    assert record[0].message.args[0] ==  f"The number of SAGE output files according to the master HDF5 file is 1. " \
        f"However, ``analyse_sage_output`` was called with {num_sage_output_files[0]}. Using the number of files " \
        "from the HDF5 file as the correct value."

# Write a test where if the output format is NOT specified, then it uses the output format specified in the sage file.
