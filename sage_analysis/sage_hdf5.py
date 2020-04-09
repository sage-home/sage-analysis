#!/usr/bin/env python
"""
This module defines the ``SageHdf5Data`` class. This class interfaces with the
:py:class:`~sage_analysis.model.Model` class to read in binary data written by **SAGE**.
The value of :py:attr:`~sage_analysis.model.Model.sage_output_format` is generally
``sage_hdf5`` if it is to be read with this class.

We refer to :doc:`../user/custom_data_classes` for more information about adding your
own Data Class to ingest data.

Author: Jacob Seiler.
"""

import warnings

from sage_analysis.model import Model
from sage_analysis.utils import read_generic_sage_params
from typing import Dict, Any, Optional
import logging
import os

import numpy as np
import h5py


logger = logging.getLogger(__name__)

## DO NOT TOUCH ##
sage_data_version = "1.00"
## DO NOT TOUCH ##

class SageHdf5Data():
    """
    Class intended to inteface with the :py:class:`~sage_analysis.model.Model` class to
    ingest the data written by **SAGE**. It includes methods for reading the output
    galaxies, setting cosmology etc. It is specifically written for when
    :py:attr:`~sage_analysis.model.Model.sage_output_format` is ``sage_hdf5``.
    """

    def __init__(self, model: Model, sage_file_to_read: str) -> None:
        """
        Instantiates the Data Class for reading in **SAGE** HDF5 data. In particular,
        opens up the file and ensures the data version matches the expected value.

        Parameters
        ----------

        model: :py:class:`~sage_analysis.model.Model` instance
            The model that this data class is associated with; this class will read the
            data for this model.
        """

        logger.info("Reading using SAGE HDF5 output format.")

        # Use the SAGE parameter file to generate a bunch of attributes.
        sage_dict = self._read_sage_params(sage_file_to_read)
        self.sage_model_dict = sage_dict
        logger.info(f"The read SAGE parameters are {sage_dict}")

        # The output data will be named via the parameter file with the ``.hdf5`` extension.
        model_path = f"{sage_dict['_base_sage_output_path']}.hdf5"
        model.hdf5_file = h5py.File(model_path, "r")

        # Due to how attributes are created in C, they will need to be decoded to get cast to a string.
        model.sage_version = model.hdf5_file["Header"]["Misc"].attrs["sage_version"].decode("ascii")
        model.sage_data_version = model.hdf5_file["Header"]["Misc"].attrs["sage_data_version"].decode("ascii")

        # Check that this module is current for the SAGE data version.
        if model.sage_data_version != sage_data_version:
            msg = "The 'sage_hdf5.py' module was written to be compatible with " \
                  "sage_data_version {0}.  Your version of SAGE HDF5 has data version " \
                  "{1}. Please update your version of SAGE or submit an issue at " \
                  "https://github.com/sage-home/sage-model/issues".format(sage_data_version, \
                  model.sage_data_version)
            raise ValueError(msg)

        # Finally, perform some checks to ensure that the data in ``model`` is compatible with our assumptions
        # regarding the HDF5 file.
        self._check_model_compatibility(model, sage_dict)

        model._num_output_files = model.hdf5_file["Header"]["Misc"].attrs["num_cores"]

        # The cores to analyze may have a default value.  To allow for this, explicitly set these as if they were read
        # from the parameter file.
        self.sage_model_dict["_first_file_to_analyse"] = 0
        self.sage_model_dict["_last_file_to_analyse"] = model._num_output_files - 1

    def _check_model_compatibility(self, model: Model, sage_dict: Optional[Dict[str, Any]]) -> None:
        """
        Ensures that the attributes in the :py:class:`~sage_analysis.model.Model` instance are compatible with the
        variables read from the **SAGE** parameter file (if read at all).

        Parameters
        ----------
        model : :py:class:`~sage_analysis.model.Model` instance
            The model that this data class is associated with.

        sage_dict : optional, dict[str, Any]
            A dictionary containing all of the fields read from the **SAGE** parameter file.

        Warnings
        --------
        UserWarning
            Raised if the user initialized ``Model`` with a value of
            :py:attr:`~sage_analysis.model.Model.num_sage_output_files` that is different to the value specified in the
            HDF5 file.
        """

        if model._num_sage_output_files is not None:
            logger.debug("It is not required to specify the number of SAGE output files when analysing HDF5 output.")

            # Check to ensure that there isn't a mismatch in the number specified and the number in the file.
            hdf5_num_files = model.hdf5_file["Header"]["Misc"].attrs["num_cores"]
            if model._num_sage_output_files != hdf5_num_files:
                warnings.warn(
                    f"The number of SAGE output files according to the master HDF5 file is {hdf5_num_files}."
                    f" However, ``analyse_sage_output`` was called with {model._num_sage_output_files}. "
                    f"Using the number of files from the HDF5 file as the correct value."
                )

    def determine_volume_analysed(self, model: Model) -> float:
        """
        Determines the volume analysed. This can be smaller than the total simulation box.

        Parameters
        ----------
        model : :py:class:`~sage_analysis.model.Model` instance
            The model that this data class is associated with.

        Returns
        -------
        volume : float
            The numeric volume being processed during this run of the code in (Mpc/h)^3.
        """

        # Each output file may have processed a different fraction of the total volume. Hence to determine the total
        # volume processed, loop through all of the files that we're analysing and use their volume fractions.
        total_volume_frac_processed = 0.0

        for core_idx in range(model._first_file_to_analyse, model._last_file_to_analyse + 1):

            core_key = "Core_{0}".format(core_idx)
            frac_processed = model.hdf5_file[core_key]["Header"]["Runtime"].attrs["frac_volume_processed"]
            total_volume_frac_processed += frac_processed
            logger.info(
                f"{core_key} processed {frac_processed} fraction of the volume. Total is {total_volume_frac_processed}"
            )

        volume = pow(model._box_size, 3) * total_volume_frac_processed
        logger.info(
            f"Total fraction of volume processed is {total_volume_frac_processed}.  Box size is "
            f"{model._box_size}. The size of the volume processed is hence {volume} (Mpc/h)^3."
        )

        return volume

    def _read_sage_params(self, sage_file_path: str) -> Dict[str, Any]:
        """
        Read the **SAGE** parameter file.

        Parameters
        ----------
        sage_file_path: string
            Path to the **SAGE** parameter file.

        Returns
        -------
        model_dict: dict [str, var]
            Dictionary containing the parameter names and their values.
        """

        model_dict = read_generic_sage_params(sage_file_path)

        return model_dict

    def determine_num_gals(self, model):
        """
        Determines the number of galaxies in all cores for this model at the specified
        :py:attr:`~sage_analysis.model.Model.snapshot`.

        Parameters
        ----------

        model: :py:class:`~sage_analysis.model.Model` class
            The :py:class:`~sage_analysis.model.Model` we're reading data for.
        """

        ngals = 0
        snap_key = "Snap_{0}".format(model.snapshot)

        for core_idx in range(model.first_file_to_analyse, model.last_file_to_analyse + 1):

            core_key = "Core_{0}".format(core_idx)
            ngals += model.hdf5_file[core_key][snap_key].attrs["num_gals"]

        model.num_gals_all_files = ngals


    def read_gals(self, model, core_num, pbar=None, plot_galaxies=False, debug=False):
        """
        Reads the galaxies of a single core at the specified
        :py:attr:`~sage_analysis.model.Model.snapshot`.

        Parameters
        ----------

        model: :py:class:`~sage_analysis.model.Model` class
            The :py:class:`~sage_analysis.model.Model` we're reading data for.

        core_num: Integer
            The core group we're reading.

        pbar: ``tqdm`` class instance, optional
            Bar showing the progress of galaxy reading.  If ``None``, progress bar will
            not show.

        plot_galaxies : Boolean, optional
            If set, plots and saves the 3D distribution of galaxies for this file.

        debug : Boolean, optional
            If set, prints out extra useful debug information.

        Returns
        -------

        gals : ``h5py`` group
            The galaxies for this file.

        Notes
        -----

        ``tqdm`` does not play nicely with printing to stdout. Hence we disable
        the ``tqdm`` progress bar if ``debug=True``.
        """

        core_key = "Core_{0}".format(core_num)
        snap_key = "Snap_{0}".format(model.snapshot)

        num_gals_read = model.hdf5_file[core_key][snap_key].attrs["num_gals"]

        # If there aren't any galaxies, exit here.
        if num_gals_read == 0:
            return None

        gals = model.hdf5_file[core_key][snap_key]

        # If we're using the `tqdm` package, update the progress bar.
        if pbar is not None:
            pbar.set_postfix(file=core_key, refresh=False)
            pbar.update(num_gals_read)

        if debug:
            print("")
            print("Core {0}, Snapshot {1} contained {2} galaxies".format(core_num,
                                                                         model.snapshot,
                                                                         num_gals_read))

            w = np.where(gals["StellarMass"][:] > 1.0)[0]
            print("{0} of these galaxies have mass greater than 10^10Msun/h".format(len(w)))

        if plot_galaxies:
            # Show the distribution of galaxies in 3D.
            from sage_analysis.plots import plot_spatial_3d

            pos = np.empty((len(gals["Posx"]), 3), dtype=np.float32)

            dim_name = ["x", "y", "z"]
            for (dim_num, dim_name) in enumerate(dim_name):
                key = "Pos{0}".format(dim_num)
                pos[:, dim_num] = gals[key][:]

            output_file = "./galaxies_{0}.{1}".format(core_num, model.plot_output_format)
            plot_spatial_3d(pos, output_file, model.box_size)

        return gals


    def update_snapshot_and_data_path(self, model: Model, snapshot: int):
        """
        Updates the :py:attr:`~sage_analysis.Model.snapshot` attribute to ``snapshot``.  As the HDF5 file contains all
        snapshot information, we do **not** need to update the path to the output data.
        """
        model._snapshot = snapshot


    def close_file(self, model):
        """
        Closes the open HDF5 file.
        """

        model.hdf5_file.close()
