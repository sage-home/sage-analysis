#!/usr/bin/env python
"""
This module defines the ``SageBinaryData`` class. This class interfaces with the
:py:class:`~sage_analysis.model.Model` class to read in binary data written by **SAGE**.
The value of :py:attr:`~sage_analysis.model.Model.sage_output_format` is generally
``sage_binary`` if it is to be read with this class.

If you wish to ingest data from your own flavour of SAGE, please open a Github issue, I plan to add this documentation
in future :)

Author: Jacob Seiler.
"""
import logging
import os
from typing import Any, Dict, Optional

import numpy as np
from tqdm import tqdm

from sage_analysis.data_class import DataClass
from sage_analysis.model import Model
from sage_analysis.utils import read_generic_sage_params

logger = logging.getLogger(__name__)


class SageBinaryData(DataClass):
    """
    Class intended to inteface with the :py:class:`~sage_analysis.model.Model` class to
    ingest the data written by **SAGE**. It includes methods for reading the output
    galaxies, setting cosmology etc. It is specifically written for when
    :py:attr:`~sage_analysis.model.Model.sage_output_format` is ``sage_binary``.
    """

    def __init__(self, model: Model, sage_file_to_read: str) -> None:
        """
        Instantiates the Data Class for reading in **SAGE** binary data. In particular,
        generates the ``numpy`` structured array to read the output galaxies.

        model: :py:class:`~sage_analysis.model.Model` instance
            The model that this data class is associated with; this class will read the
            data for this model.
        """

        logger.info("Reading using SAGE binary output format.")

        self._get_galaxy_struct()

        # Use the SAGE parameter file to generate a bunch of attributes.
        sage_dict = self.read_sage_params(sage_file_to_read)
        self.sage_model_dict = sage_dict
        logger.info(f"The read SAGE parameters are {sage_dict}")

    def determine_volume_analyzed(self, model: Model) -> float:
        """
        Determines the volume analyzed. This can be smaller than the total simulation box.

        Parameters
        ----------
        model : :py:class:`~sage_analysis.model.Model` instance
            The model that this data class is associated with.

        Returns
        -------
        volume : float
            The numeric volume being processed during this run of the code in (Mpc/h)^3.
        """

        # To properly scale properties that use the simulation volume (e.g., SMF), we need
        # to know how much of the volume this model is analysing.  SAGE is formulated such
        # that every processor writes out a single file.  However, each model here can
        # analyze fewer files than were simulated by SAGE.

        # For example, SAGE could have run on 4 processors, and hence 4 files would be
        # produced.  To allow the quick inspection of results, we may only be running our
        # analysis on one file. Hence, we should scale some of the properties by a factor
        # of 4.

        # Importantly: There is no way for the binary output of SAGE to know this factor!
        # Hence, if the user is running in binary mode, they MUST specify the total number
        # of files that SAGE output (i.e., the number of processors they ran SAGE with).
        frac_volume_analyzed = (model._last_file_to_analyze - model._first_file_to_analyze + 1) / \
            model._num_sage_output_files
        volume = pow(model._box_size, 3) * frac_volume_analyzed

        logger.info(
            f"The files read is [{model._first_file_to_analyze}, {model._last_file_to_analyze}] with a total number "
            f"of {model._num_sage_output_files}; resulting a volume fraction analyzed of {frac_volume_analyzed}.\nThe "
            f"box size is {model._box_size} (Mpc/h) yielding a analyzed volume of {volume} (Mpc/h)^3."
        )

        return volume

    def read_sage_params(self, sage_file_path: str) -> Dict[str, Any]:
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

    def _get_galaxy_struct(self):
        """
        Sets the ``numpy`` structured array for holding the galaxy data.
        """

        galdesc_full = [
            ("SnapNum", np.int32),
            ("Type", np.int32),
            ("GalaxyIndex", np.int64),
            ("CentralGalaxyIndex", np.int64),
            ("SAGEHaloIndex", np.int32),
            ("SAGETreeIndex", np.int32),
            ("SimulationHaloIndex", np.int64),
            ("mergeType", np.int32),
            ("mergeIntoID", np.int32),
            ("mergeIntoSnapNum", np.int32),
            ("dT", np.float32),
            ("Pos", (np.float32, 3)),
            ("Vel", (np.float32, 3)),
            ("Spin", (np.float32, 3)),
            ("Len", np.int32),
            ("Mvir", np.float32),
            ("CentralMvir", np.float32),
            ("Rvir", np.float32),
            ("Vvir", np.float32),
            ("Vmax", np.float32),
            ("VelDisp", np.float32),
            ("ColdGas", np.float32),
            ("StellarMass", np.float32),
            ("BulgeMass", np.float32),
            ("HotGas", np.float32),
            ("EjectedMass", np.float32),
            ("BlackHoleMass", np.float32),
            ("IntraClusterStars", np.float32),
            ("MetalsColdGas", np.float32),
            ("MetalsStellarMass", np.float32),
            ("MetalsBulgeMass", np.float32),
            ("MetalsHotGas", np.float32),
            ("MetalsEjectedMass", np.float32),
            ("MetalsIntraClusterStars", np.float32),
            ("SfrDisk", np.float32),
            ("SfrBulge", np.float32),
            ("SfrDiskZ", np.float32),
            ("SfrBulgeZ", np.float32),
            ("DiskRadius", np.float32),
            ("Cooling", np.float32),
            ("Heating", np.float32),
            ("QuasarModeBHaccretionMass", np.float32),
            ("TimeOfLastMajorMerger", np.float32),
            ("TimeOfLastMinorMerger", np.float32),
            ("OutflowRate", np.float32),
            ("infallMvir", np.float32),
            ("infallVvir", np.float32),
            ("infallVmax", np.float32)
        ]
        names = [galdesc_full[i][0] for i in range(len(galdesc_full))]
        formats = [galdesc_full[i][1] for i in range(len(galdesc_full))]
        galdesc = np.dtype({"names": names, "formats": formats}, align=True)

        self.galaxy_struct = galdesc

    def determine_num_gals(self, model: Model, *args):
        """
        Determines the number of galaxies in all files for this
        :py:class:`~sage_analysis.model.Model`.

        Parameters
        ----------

        model: :py:class:`~sage_analysis.model.Model` class
            The :py:class:`~sage_analysis.model.Model` we're reading data for.

        *args : Any
            Extra arguments to allow other data class to pass extra arguments to their version of
            ``determine_num_gals``.
        """

        num_gals = 0

        for file_num in range(model._first_file_to_analyze, model._last_file_to_analyze+1):

            fname = self._check_for_file(model, file_num)
            if fname is None:
                logger.debug(f"File\t{fname} \tdoes not exist!")
                continue

            with open(fname, "rb") as f:
                _ = np.fromfile(f, np.dtype(np.int32), 1)[0]
                num_gals_file = np.fromfile(f, np.dtype(np.int32), 1)[0]

                num_gals += num_gals_file

        model._num_gals_all_files = num_gals

    def read_gals(
        self,
        model: Model,
        file_num: int,
        snapshot: int,
        pbar: Optional[tqdm] = None,
        plot_galaxies: bool = False,
        debug: bool = False
    ):
        """
        Reads the galaxies of a model file at snapshot specified by
        :py:attr:`~sage_analysis.model.Model.snapshot`.

        Parameters
        ----------

        model: :py:class:`~sage_analysis.model.Model` class
            The :py:class:`~sage_analysis.model.Model` we're reading data for.

        file_num: int
            Suffix number of the file we're reading.

        pbar: ``tqdm`` class instance, optional
            Bar showing the progress of galaxy reading.  If ``None``, progress bar will
            not show.

        plot_galaxies: bool, optional
            If set, plots and saves the 3D distribution of galaxies for this file.

        debug: bool, optional
            If set, prints out extra useful debug information.

        Returns
        -------
        gals : ``numpy`` structured array with format given by :py:method:`~_get_galaxy_struct`
            The galaxies for this file.

        Notes
        -----
        ``tqdm`` does not play nicely with printing to stdout. Hence we disable
        the ``tqdm`` progress bar if ``debug=True``.
        """

        fname = self._check_for_file(model, file_num)
        if fname is None:
            logger.debug(f"File\t{fname} \tdoes not exist!")
            return None

        with open(fname, "rb") as f:
            # First read the header information.
            Ntrees = np.fromfile(f, np.dtype(np.int32), 1)[0]
            num_gals = np.fromfile(f, np.dtype(np.int32), 1)[0]
            _ = np.fromfile(f, np.dtype((np.int32, Ntrees)), 1)

            # If there aren't any galaxies, exit here.
            if num_gals == 0:
                return None

            # Then the actual galaxies.
            gals = np.fromfile(f, self.galaxy_struct, num_gals)

            # If we're using the `tqdm` package, update the progress bar.
            if pbar is not None:
                pbar.set_postfix(file=fname, refresh=False)
                pbar.update(num_gals)

        if debug:
            print("")
            print(f"File {fname} contained {Ntrees} trees with {num_gals} galaxies")

            w = np.where(gals["StellarMass"] > 1.0)[0]
            print(f"{len(w)} of these galaxies have mass greater than 10^10Msun/h")

        if plot_galaxies:

            from sage_analysis.plots import plot_spatial_3d

            # Show the distribution of galaxies in 3D.
            pos = gals["Pos"][:]
            output_file = f"./galaxies_{file_num}.{model.plot_output_format}"
            plot_spatial_3d(pos, output_file, self.box_size)

        # For the HDF5 file, some data sets have dimensions Nx1 rather than Nx3
        # (e.g., Position). To ensure the galaxy data format is identical to the binary
        # output, we will split the binary fields into Nx1. This is simpler than creating
        # a new dataset within the HDF5 regime.
        from numpy.lib import recfunctions as rfn
        multidim_fields = ["Pos", "Vel", "Spin"]
        dim_names = ["x", "y", "z"]

        for field in multidim_fields:
            for dim_num, dim_name in enumerate(dim_names):
                dim_field = f"{field}{dim_name}"
                gals = rfn.rec_append_fields(gals, dim_field,
                                             gals[field][:, dim_num])

        return gals

    def update_snapshot_and_data_path(self, model: Model, snapshot: int, use_absolute_path: bool = False):
        """
        Updates the :py:attr:`~sage_analysis.model.Model._sage_data_path` to point to a new redshift file. Uses the
        redshift array :py:attr:`~sage_analysis.model.Model.redshifts`.

        Parameters
        ----------
        snapshot : int
            Snapshot we're updating :py:attr:`~sage_analysis.model.Model._sage_data_path` to
            point to.

        use_absolute_path : bool
            If specified, will use the absolute path to the **SAGE** output data. Otherwise, will use the path that is
            relative to the **SAGE** parameter file.  This is hand because the **SAGE** parameter file can contain
            either relative or absolute paths.
        """

        model._snapshot = snapshot

        new_redshift = model.redshifts[snapshot]

        # The parameter file could refer to the absolute path or the relative path, so be careful.
        if use_absolute_path:
            base_path = model._base_sage_output_path_absolute
        else:
            base_path = model._base_sage_output_path_relative

        model._sage_data_path = f"{base_path}_z{new_redshift:.3f}"

    def _check_for_file(self, model: Model, file_num: int) -> Optional[str]:
        """
        Checks to see if a file for the given file number exists.  Importantly, we check assuming that the path given
        in the **SAGE** parameter file is **relative** and **absolute**.

        Parameters
        ----------
        file_num : int
            The file number that we're checking for files.

        Returns
        -------
        fname or ``None``
            If a file exists, the name of that file.  Otherwise, if the file does not exist (using either relative or
            absolute paths), then ``None``.
        """

        # Try relative, and then absolute paths.
        for use_absolute_path in [False, True]:

            self.update_snapshot_and_data_path(model, model._snapshot, use_absolute_path=use_absolute_path)

            fname = f"{model._sage_data_path}_{file_num}"
            if os.path.isfile(fname):
                return fname
        return None

    def close_file(self, model: Model):
        """
        An empty method to ensure consistency with the HDF5 data class. This is empty because snapshots are saved over
        different files by default in the binary format.
        """
        pass
