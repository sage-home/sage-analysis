#!/usr/bin/env python
"""
This module defines the ``SageBinaryData`` class. This class interfaces with the
:py:class:`~sage_analysis.model.Model` class to read in binary data written by **SAGE**.
The value of :py:attr:`~sage_analysis.model.Model.sage_output_format` is generally
``sage_binary`` if it is to be read with this class.

We refer to :doc:`../user/custom_data_classes` for more information about adding your
own Data Class to ingest data.

Author: Jacob Seiler.
"""

from sage_analysis.model import Model

import numpy as np
import os

class SageBinaryData():
    """
    Class intended to inteface with the :py:class:`~sage_analysis.model.Model` class to
    ingest the data written by **SAGE**. It includes methods for reading the output
    galaxies, setting cosmology etc. It is specifically written for when
    :py:attr:`~sage_analysis.model.Model.sage_output_format` is ``sage_binary``.
    """

    def __init__(self, model, num_output_files, sage_file_to_read=None,
                 snapshot_to_use=None):
        """
        Instantiates the Data Class for reading in **SAGE** binary data. In particular,
        generates the ``numpy`` structured array to read the output galaxies.

        model: :py:class:`~sage_analysis.model.Model` instance
            The model that this data class is associated with; this class will read the
            data for this model.

        num_output_files: int
            The number of files that **SAGE** wrote out when simulating ``model``. This is
            generally equal to the number of processors used to run **SAGE**.

        sage_file_to_read: string, optional
            Specifies the **SAGE** file to be read and used to update the
            ``sage_model_dict`` attribute with the parameters specified inside.  If set
            to ``None``, does not update this attribute.  Instead, the user must provide
            all the parameters to analyze the data to the
            :py:meth:`~sage_analysis.model.Model.update_attributes`.

        snapshot_to_use: int, optional
            The snapshot number being analysed for this ``model``. If reading a **SAGE**
            file (``sage_file_to_read`` is not ``None``), this must be specified, or
            :py:meth:`~sage_analysis.sage_binary.SageBinaryData.update_snapshot` must
            be called before any reading of data is done.
        """

        self.get_galaxy_struct()

        # To properly scale properties that use the simulation volume (e.g., SMF), we need
        # to know how much of the volume this model is analysing.  SAGE is formulated such
        # that every processor writes out a single file.  However, each model here can
        # analyse fewer files than were simulated by SAGE.

        # For example, SAGE could have run on 4 processors, and hence 4 files would be
        # produced.  To allow the quick inspection of results, we may only be running our
        # analysis on one file. Hence, we should scale some of the properties by a factor
        # of 4.

        # Importantly: There is no way for the binary output of SAGE to know this factor!
        # Hence, if the user is running in binary mode, they MUST specify the total number
        # of files that SAGE output (i.e., the number of processors they ran SAGE with).
        model._num_output_files = num_output_files

        # Need the snapshot to specify the name of the SAGE output file (if we're reading
        # the SAGE ini file).
        if sage_file_to_read:
            self._snapshot = snapshot_to_use

        # Use the SAGE parameter file to generate a bunch of attributes.
        if sage_file_to_read:
            sage_dict = self.read_sage_params(sage_file_to_read, snapshot_to_use)
            self.sage_model_dict = sage_dict


    def read_sage_params(self, sage_file_path, snapshot_to_use=None):
        """
        Reads the **SAGE** parameter file values.

        Parameters
        ----------

        sage_file_path: string
            Path to the **SAGE** parameter file.

        snapshot_to_use: int, optional
            The snapshot that this model is reading.  If this is not specified,
            :py:meth:`~sage_analysis.sage_binary.SageBinaryData.update_snapshot` must be
            called before any reading of data is done.

        Returns
        -------

        model_dict: dict [string, variable], optional
            Dictionary containing the parameter values for this class instance. Attributes
            of the class are set with name defined by the key with corresponding values.

        FileNotFoundError
            Raised if the specified **SAGE** parameter file is not found.
        """

        # Fields that we will be reading from the ini file.
        SAGE_fields = [ "FileNameGalaxies",
                        "OutputDir",
                        "FirstFile",
                        "LastFile",
                        "OutputFormat",
                        "NumSimulationTreeFiles",
                        "FileWithSnapList",
                        "Hubble_h",
                        "BoxSize",
                        "PartMass"]
        SAGE_dict = {}

        # Ignore lines starting with one of these.
        comment_characters = [";", "%", "-"]

        try:
            with open(sage_file_path, "r") as SAGE_file:
                data = SAGE_file.readlines()

                # Each line in the parameter file is of the form...
                # parameter_name       parameter_value.
                for line in range(len(data)):

                    # Remove surrounding whitespace from the line.
                    stripped = data[line].strip()

                    # May have been an empty line.
                    try:
                        first_char = stripped[0]
                    except IndexError:
                        continue

                    # Check for comment.
                    if first_char in comment_characters:
                        continue

                    # Split into [name, value] list.
                    split = stripped.split()

                    # Then check if the field is one we care about.
                    if split[0] in SAGE_fields:

                        SAGE_dict[split[0]] = split[1]

        except FileNotFoundError:
            raise FileNotFoundError("Could not file SAGE ini file {0}".format(fname))

        # Now we have all the fields, rebuild the dictionary to be exactly what we need for
        # initialising the model.
        model_dict = {}

        alist = np.loadtxt(SAGE_dict["FileWithSnapList"])
        redshifts = 1.0 / alist - 1.0
        model_dict["redshifts"] = redshifts

        if snapshot_to_use:
            output_tag = f"_z{redshifts[snapshot_to_use]:.3f}"
        else:
            print("A data class was instantiated without specifying an initial "
                  "snapshot to read. This is allowed, but `Data_Class.update_snapshot` "
                  "must be called before any reading is done.")
            output_tag == "NOT_SET"

        model_path = f"{SAGE_dict['OutputDir']}/{SAGE_dict['FileNameGalaxies']}" \
                     f"{output_tag}"
        model_dict["_model_path"] = model_path

        model_dict["_output_path"] = f"{SAGE_dict['OutputDir']}/plots/"

        model_dict["_hubble_h"] = float(SAGE_dict["Hubble_h"])
        model_dict["_box_size"] = float(SAGE_dict["BoxSize"])
        model_dict["_num_sim_tree_files"] = int(SAGE_dict["NumSimulationTreeFiles"])

        return model_dict


    def get_galaxy_struct(self):
        """
        Sets the ``numpy`` structured array for holding the galaxy data.
        """

        galdesc_full = [
            ("SnapNum"                      , np.int32),
            ("Type"                         , np.int32),
            ("GalaxyIndex"                  , np.int64),
            ("CentralGalaxyIndex"           , np.int64),
            ("SAGEHaloIndex"                , np.int32),
            ("SAGETreeIndex"                , np.int32),
            ("SimulationHaloIndex"          , np.int64),
            ("mergeType"                    , np.int32),
            ("mergeIntoID"                  , np.int32),
            ("mergeIntoSnapNum"             , np.int32),
            ("dT"                           , np.float32),
            ("Pos"                          , (np.float32, 3)),
            ("Vel"                          , (np.float32, 3)),
            ("Spin"                         , (np.float32, 3)),
            ("Len"                          , np.int32),
            ("Mvir"                         , np.float32),
            ("CentralMvir"                  , np.float32),
            ("Rvir"                         , np.float32),
            ("Vvir"                         , np.float32),
            ("Vmax"                         , np.float32),
            ("VelDisp"                      , np.float32),
            ("ColdGas"                      , np.float32),
            ("StellarMass"                  , np.float32),
            ("BulgeMass"                    , np.float32),
            ("HotGas"                       , np.float32),
            ("EjectedMass"                  , np.float32),
            ("BlackHoleMass"                , np.float32),
            ("IntraClusterStars"            , np.float32),
            ("MetalsColdGas"                , np.float32),
            ("MetalsStellarMass"            , np.float32),
            ("MetalsBulgeMass"              , np.float32),
            ("MetalsHotGas"                 , np.float32),
            ("MetalsEjectedMass"            , np.float32),
            ("MetalsIntraClusterStars"      , np.float32),
            ("SfrDisk"                      , np.float32),
            ("SfrBulge"                     , np.float32),
            ("SfrDiskZ"                     , np.float32),
            ("SfrBulgeZ"                    , np.float32),
            ("DiskRadius"                   , np.float32),
            ("Cooling"                      , np.float32),
            ("Heating"                      , np.float32),
            ("QuasarModeBHaccretionMass"    , np.float32),
            ("TimeOfLastMajorMerger"        , np.float32),
            ("TimeOfLastMinorMerger"        , np.float32),
            ("OutflowRate"                  , np.float32),
            ("infallMvir"                   , np.float32),
            ("infallVvir"                   , np.float32),
            ("infallVmax"                   , np.float32)
            ]
        names = [galdesc_full[i][0] for i in range(len(galdesc_full))]
        formats = [galdesc_full[i][1] for i in range(len(galdesc_full))]
        galdesc = np.dtype({"names":names, "formats":formats}, align=True)

        self.galaxy_struct = galdesc


    def determine_num_gals(self, model):
        """
        Determines the number of galaxies in all files for this
        :py:class:`~sage_analysis.model.Model`.

        Parameters
        ----------

        model: :py:class:`~sage_analysis.model.Model` class
            The :py:class:`~sage_analysis.model.Model` we're reading data for.
        """

        num_gals = 0

        for file_num in range(model.first_file, model.last_file+1):

            fname = "{0}_{1}".format(model.model_path, file_num)

            if not os.path.isfile(fname):
                print("File\t{0} \tdoes not exist!".format(fname))
                raise FileNotFoundError

            with open(fname, "rb") as f:
                Ntrees = np.fromfile(f, np.dtype(np.int32),1)[0]
                num_gals_file = np.fromfile(f, np.dtype(np.int32),1)[0]

                num_gals += num_gals_file

        model.num_gals_all_files = num_gals


    def read_gals(self, model, file_num, pbar=None, plot_galaxies=False, debug=False):
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

        gals : ``numpy`` structured array with format given by ``get_galaxy_struct()``
            The galaxies for this file.

        Notes
        -----

        ``tqdm`` does not play nicely with printing to stdout. Hence we disable
        the ``tqdm`` progress bar if ``debug=True``.
        """

        fname = "{0}_{1}".format(model.model_path, file_num)

        # We allow the skipping of files.  If we skip, don't increment a counter.
        if not os.path.isfile(fname):
            print("File\t{0} \tdoes not exist!".format(fname))
            return None

        with open(fname, "rb") as f:
            # First read the header information.
            Ntrees = np.fromfile(f, np.dtype(np.int32),1)[0]
            num_gals = np.fromfile(f, np.dtype(np.int32),1)[0]
            gals_per_tree = np.fromfile(f, np.dtype((np.int32, Ntrees)), 1)

            # If there aren't any galaxies, exit here.
            if num_gals == 0:
                return None

            # Then the actual galaxies.
            gals = np.fromfile(f, self.galaxy_struct, num_gals)

            # If we're using the `tqdm` package, update the progress bar.
            if pbar:
                pbar.set_postfix(file=fname, refresh=False)
                pbar.update(num_gals)

        if debug:
            print("")
            print("File {0} contained {1} trees with {2} galaxies".format(fname, Ntrees, num_gals))

            w = np.where(gals["StellarMass"] > 1.0)[0]
            print("{0} of these galaxies have mass greater than 10^10Msun/h".format(len(w)))

        if plot_galaxies:

            from sage_analysis.plots import plot_spatial_3d

            # Show the distribution of galaxies in 3D.
            pos = gals["Pos"][:]
            output_file = "./galaxies_{0}.{1}".format(file_num, model.plot_output_format)
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
                dim_field = "{0}{1}".format(field, dim_name)
                gals = rfn.rec_append_fields(gals, dim_field,
                                             gals[field][:, dim_num])

        return gals


    def update_snapshot(self, model, snapshot):
        """
        Updates the :py:attr:`~sage_analysis.model.Model.model_path` to point to a new
        redshift file. Uses the redshift array
        :py:attr:`~sage_analysis.model.Model.redshifts`.

        Parameters
        ----------

        snapshot: int
            Snapshot we're updating :py:attr:`~sage_analysis.model.Model.model_path` to
            point to.
        """

        model._snapshot = snapshot

        # First get the new redshift.
        new_redshift = model.redshifts[snapshot]

        # We can't be guaranteed the `model_path` has been set just yet.
        try:
            _ = self.set_redshift
        except AttributeError:
            # model_path hasn't been set so use the full form.
            model.model_path = "{0}_z{1:.3f}".format(self.model_path, new_redshift)
        else:
            # model_path is of the form "<Initial/Path/To/File_zXXX.XXX>"
            # The number of decimal places is always 3.
            # The number of characters before "z" is arbitrary, we could be at z8.539 or z127.031.
            # Hence walk backwards through the model path until we reach a "z".
            letters_from_end = 0
            letter = self.model_path[-(letters_from_end+1)]
            while letter != "z":
                letters_from_end += 1
                letter = self.model_path[-(letters_from_end+1)]

            # Then truncate there and append the new redshift.
            model_path_prefix = self.model_path[:-(letters_from_end+1)]
            model.model_path = "{0}z{1:.3f}".format(model_path_prefix, new_redshift)
