Analysing Across Multiple Snapshots
===================================

We show how to analyse the output of **SAGE** at a single snapshot :doc:`here
<./analysing_sage>`.  On this page, we show how to analyse **SAGE** output
across multiple snapshots.  This is very useful if you wish to analyse the
evolution of (e.g.,) the stellar mass function, the stellar mass density, etc.

This full example is shown in the `history module`_ using the default
parameter file.

Setting Things Up
-----------------

In a similar manner to analysing a :doc:`single snapshot <./analysing_sage>`, we
first specify which properties we wish to analyse and plot.

.. code-block:: python

    # Base specifications.
    plot_output_format = "png"
    plot_output_path = "./plots"  # Will be created if path doesn't exist.

    plot_toggles = {"SMF" : 1,  # Stellar mass function across redshift.
                    "SFRD" : 1,  # Star formation rate density across redshift.
                    "SMD" : 1}  # Stellar mass density across redshift.


Then, we specify the **SAGE** output we wish to analyse and the redshifts at
which we want the properties to be calculated at.

**NOTE:** For all the specified redshifts, **sage-analysis** searches for the
snapshot closest to the redshift.  Also set the entry to ``"All"`` to analyse
across all snapshots.

.. code-block:: python

    millennium = {"SMF_z" : [0.0, 1.0, 2.0, 3.0],  # Redshifts you wish to plot the stellar mass function at.
                  "density_z": "All", # Redshifts to plot the stellar mass/star formation density at.
                  "IMF": "Chabrier",  # Chabrier or Salpeter.
                  "label": "Mini-Millennium",  # Legend label.
                  "sage_file": "../input/millennium.par",
                  "sage_output_format": "sage_hdf5",
                  "first_file": 0,  # File range (or core range for HDF5) to plot.
                  "last_file": 0,  # Closed interval, [first_file, last_file].
                 }

**NOTE:** If the ``sage_output_format`` is ``sage_binary`` (i.e., **SAGE**
wrote as binary output), then you must also specify the number of output files,
:attr:`~sage_analysis.model.Model.num_output_files`.

Initializing the Model
--------------------

Setting up the Class
~~~~~~~~~~~~~~~~~~~~

The analysis of **SAGE** models is done through a specialized
:class:`~sage_analysis.model.Model` class. **Importantly,** the Model class only
handles the calculating properties.  To actually read the **SAGE** output, each
Model requires a data class.  These are specific to
the **SAGE** output format.  For example, we include a data class for
:class:`~sage_analysis.sage_hdf5.SageHdf5Data` and
:class:`~sage_analysis.sage_binary.SageBinaryData`.
Through this data class, the package can be easily extended to ingest
any arbitrary **SAGE** output format.  We show such an example
:doc:`here <./custom_data_classes>`.

.. code-block:: python

    from sage_analysis.model import Model
    from sage_analysis.sage_hdf5 import SageHdf5Data

    model = Model()
    model.plot_output_format = plot_output_format

    model.data_class = SageHdf5Data(model, millennium["sage_file"])

    # The data class has read the SAGE ini file.  Update the model with the parameters
    # read and those specified by the user.
    model.update_attributes(model_dict)

Specifying the Empty Property Containers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We also initialize the Model properties as outlined :doc:`previously
<./analysing_sage.rst>`.

.. code-block:: python

    stellar_properties = ["SMF", "red_SMF", "blue_SMF"]
    my_model.init_binned_properties(8.0, 12.0, 0.1, "stellar_mass_bins",
                                    stellar_properties)

    # Properties that are extended as lists.
    scatter_properties = []
    my_model.init_scatter_properties(scatter_properties)

    # Properties that are stored as a single number.
    single_properties = ["SMFD", "SFRD"]
    my_model.init_single_properties(single_properties)

    # We will store the values of each snapshot in a dictionary.
    model.properties["SMF_dict"] = {}
    model.properties["SFRD_dict"] = {}
    model.properties["SMD_dict"] = {}

Calculation and Plotting Dictionaries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Again, as outlined :doc:`previously <./analysing_sage.rst>`, we also generate
the dictionaries necessary to analyse and plot properties.

.. code-block:: python

    from sage_analysis.utils import generate_func_dict

    # Search for functions named "calc_<plot_toggle>" in the "example_calcs"
    # module.
    calculation_functions = generate_func_dict(
                                plot_toggles,
                                module_name="sage_analysis.example_calcs",
                                function_prefix="calc"
                                )

    # Search for functions named "plot_<plot_toggles>" in the "example_plots"
    # module.
    plot_functions = generate_func_dict(
                        plot_toggles,
                        module_name="sage_analysis.example_plots",
                        function_prefix="plot_"
                        )

Setting Up The Snapshot Loop
----------------------------

The key difference for this example is that we want to analyse properties over
a number of redshifts.  We hence must determine which snapshots in the model
correspond to the requested redshifts.

.. code-block:: python

    # We may be plotting the density at all snapshots...
    if model_dict["density_z"] == "All":
        model.density_redshifts = model.redshifts
    else:
        model.density_redshifts = model_dict["density_z"]

    # Same for SMF
    if model_dict["SMF_z"] == "All":
        model.SMF_redshifts = model.redshifts
    else:
        model.SMF_redshifts = model_dict["SMF_z"]

    # Find the snapshots that most closely match the requested redshifts.
    model.SMF_snaps = [(np.abs(model.redshifts - SMF_redshift)).argmin() for
                       SMF_redshift in model.SMF_redshifts]

    model.density_snaps = [(np.abs(model.redshifts - density_redshift)).argmin() for
                           density_redshift in model.density_redshifts]

    # Check which snapshots we uniquely need to loop through.
    snaps_to_loop = np.unique(my_model.SMF_snaps + my_model.density_snaps)

Iterating Through Snapshots
---------------------------

Finally, we are poised to iterate through the snapshots and calculate all the
properties required. Importantly, at the end of each snapshot, we must place
the calculate properties into the appropriate dictionary and reset the
property.

.. code-block:: python

    for snap in snap_iter:

        # Each snapshot is unique. So reset the tracking.
        model.properties["SMF"] = np.zeros(len(model.bins["stellar_mass_bins"])-1,
                                           dtype=np.float64)
        model.properties["SFRD"] = 0.0
        model.properties["SMD"] = 0.0

        # Update the snapshot we're reading from. Data Class specific.
        model.data_class.update_snapshot(model, snap)

        # Calculate all the properties. Since we're using a HDF5 file, we want to keep
        # the file open because we read other snapshots from that one file.
        model.calc_properties_all_files(calculation_functions, close_file=False)

        # We need to place the SMF inside the dictionary to carry through.
        if snap in model.SMF_snaps:
            model.properties["SMF_dict"][snap] = model.properties["SMF"]

        # Same with the densities.
        if snap in model.density_snaps:

            model.properties["SFRD_dict"][snap] = model.properties["SFRD"]
            model.properties["SMD_dict"][snap] = model.properties["SMD"]

    # Close the HDF5 file cause we're done with it.
    model.data_class.close_file(model)

Finally, plot the properties!

.. code-block:: python

    # Similar to the calculation functions, all of the plotting functions are in the
    # `example_plots_history.py` module and are labelled `plot_<toggle>`.
    plot_functions = generate_func_dict(plot_toggles,
                                        module_name="sage_analysis.example_plots_history",
                                        function_prefix="plot_")

    # Now do the plotting.
    for func_name in plot_functions.keys():
        func = plot_functions[func_name][0]
        keyword_args = plot_functions[func_name][1]

        func(models, plot_output_path, plot_output_format, **keyword_args)

This produces the stellar mass function, star formation rate density, and
stellar mass density over the various redshifts.

|SMF_redshift| |SFRD| |SMD|

.. |SMF_redshift| image:: ../figs/SMF_redshift.png
.. |SFRD| image:: ../figs/SFRD.png
.. |SMD| image:: ../figs/SMD.png

.. _history module: https://github.com/sage-home/sage-model/plotting/galaxy_properties.py
