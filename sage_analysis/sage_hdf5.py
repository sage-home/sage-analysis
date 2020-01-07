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

import numpy as np
import h5py

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

    def __init__(self, model, sage_file_to_read=None, model_path=None):
        """
        Instantiates the Data Class for reading in **SAGE** HDF5 data. In particular,
        opens up the file and ensures the data version matches the expected value.

        Parameters
        ----------

        model: :py:class:`~sage_analysis.model.Model` instance
            The model that this data class is associated with; this class will read the
            data for this model.

        sage_file_to_read: string, optional
            Specifies the **SAGE** file to be read and used to update the
            ``sage_model_dict`` attribute with the parameters specified inside.  If set
            to ``None``, does not update this attribute.  Instead, the user must provide
            all the parameters to analyze the data to the
            :py:meth:`~sage_analysis.model.Model.update_attributes`.

        model_path: string, optional
            Path to the master **SAGE** output file.  This must be specified
            only if the **SAGE** ``.ini`` file is not read (``sage_file_to_read`` is
            ``None``).
        """

        # Use the SAGE parameter file to generate a bunch of attributes.
        if sage_file_to_read:
            sage_dict = self.read_sage_params(sage_file_to_read)
            self.sage_model_dict = sage_dict

        # The user may not have read the SAGE parameter file.  If so, they needed to have
        # specifed the path to the model file.
        if sage_file_to_read:
            my_model_path = sage_dict["_model_path"]
        elif model_path:
            my_model_path = model_path
        else:
            print("A Model data class was instantiated without either specifying the "
                  "SAGE parameter file or a specific output file to read.  One of these "
                  "MUST be specified.")
            print(model)

        model.hdf5_file = h5py.File(my_model_path, "r")

        # Due to how attributes are created in C, they will need to be decoded to get cast
        # to a string.
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

        model.num_output_files = model.hdf5_file["Header"]["Misc"].attrs["num_cores"]


    def read_sage_params(self, sage_file_path):
        """
        Reads the **SAGE** parameter file values.

        Parameters
        ----------

        sage_file_path: string
            Path to the **SAGE** parameter file.

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

        model_path = f"{SAGE_dict['OutputDir']}/{SAGE_dict['FileNameGalaxies']}" \
                     f".hdf5"

        model_dict["_model_path"] = model_path

        model_dict["_output_path"] = f"{SAGE_dict['OutputDir']}/plots/"

        model_dict["_hubble_h"] = float(SAGE_dict["Hubble_h"])
        model_dict["_box_size"] = float(SAGE_dict["BoxSize"])
        model_dict["_num_sim_tree_files"] = int(SAGE_dict["NumSimulationTreeFiles"])

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

        for core_idx in range(model.first_file, model.last_file + 1):

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


    def update_snapshot(self, model, snapshot):
        """
        Updates the :py:attr:`~sage_analysis.Model.snapshot` attribute to ``snapshot``.
        """
        model._snapshot = snapshot


    def close_file(self, model):
        """
        Closes the open HDF5 file.
        """

        model.hdf5_file.close()
