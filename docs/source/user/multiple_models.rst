Analysing and Plotting Multiple Models
======================================

A key feature of the **sage-analysis** package is its ability to easily and
succinctly analyse and plotting a number of different **SAGE** Models.  This
makes it easy to compare different models very quickly.

We have covered the basics of plotting a single model :doc:`here
<./analysing_sage>`.  In this example, we touch upon those things that must
change to visualize multiple models.

Setting Things Up
-----------------

Again, specify what you want to plot in an identical manner to previously.

.. code-block:: python

    # Going to just plot the stellar mass function.
    plot_toggles = {"SMF": 1}

    plot_output_format = "png"
    plot_output_path = "./plots"

    if not os.path.exists(plot_output_path):
        os.makedirs(plot_output_path)

Defining the Model Dictionaries
-------------------------------

Each Model requires its own dictionary to define its input files and other
parameters.  For this example, let's suppose we wish to compare the results of
running **SAGE** on the default Mini-Millennium simulation and Mini-Millennium
without supernova feedback.

.. code-block:: python

    millennium = { "snapshot": 63,   # Snapshot we're plotting properties at.
                   "IMF": "Chabrier",  # Chabrier or Salpeter.
                   "label": "Mini-Millennium",  # Legend label.
                   "sage_file": "../input/millennium.par",
                   "sage_output_format": "sage_hdf5",
                   "first_file": 0,  # File range (or core range for HDF5) to plot.
                   "last_file": 0,  # Closed interval, [first_file, last_file].
                 }

    millennium_noSN = { "snapshot": 63,   # Snapshot we're plotting properties at.
                        "IMF": "Chabrier",  # Chabrier or Salpeter.
                        "label": "Mini-Millennium-NoSN",  # Legend label.
                        "sage_file": "../input/millennium_noSN.par",
                        "sage_output_format": "sage_binary",
                        "first_file": 0,  # File range (or core range for HDF5) to plot.
                        "last_file": 0,  # Closed interval, [first_file, last_file].
                        "num_output_files": 1,
                 }

    sims_to_plot = [millennium, millennium_noSN]
    models = []

The important line here is the very last line.  We place all the models we wish
to analyse into a single list, ``sims_to_plot``.  Then, in the following code
block, we will iterate through each of these simulations and place the Model
classes in the ``model`` list.

**NOTE:** For example purposes, we use the binary **SAGE** output option for
the no supernovae simulation.  Hence, we must specify the number of output
files that **SAGE** wrote to, :attr:`~sage_analysis.model.Model.num_output_files`.

Iterating Through Models
------------------------

Analysing multiple models is down in an identical manner as shown
:doc:`previously <./analysing_sage>`.  The only difference is that we wrap the
code in a ``for`` loop and iterate through ``models_to_plot``. We defer to the
:doc:`previous page <./analysing_sage>` for a full explanation of the following
code block.

.. code-block:: python

    # Set up the dictionaries for the calculation and plotting functions.
    from sage_analysis.utils import generate_func_dict

    # Search for functions named "calc_<plot_toggle>" in the "example_calcs"
    # module and "plot_<plot_toggle>" in the "example_plots" module.
    calculation_functions = generate_func_dict(
                                plot_toggles,
                                module_name="sage_analysis.example_calcs",
                                function_prefix="calc"
                                )

    plot_functions = generate_func_dict(
                        plot_toggles,
                        module_name="sage_analysis.example_plots",
                        function_prefix="plot_"
                        )

    # Iterate through the simulations and set up a Model for each.
    from sage_analysis.model import Model
    from sage_analysis.sage_hdf5 import SageHdf5Data
    from sage_analysis.sage_binary import SageBinaryData

    for model_dict in sims_to_plot:
        model = Model()
        model.plot_output_format = plot_output_format

        # This switch case should be extended if you're defining your own
        # custom data class.
        if model_dict["sage_output_format"] == "sage_hdf5":
            model.data_class = SageHdf5Data(model, millennium["sage_file"])
        elif model_dict["sage_output_format"] == "sage_binary":
            model.data_class = SageBinaryData(model, model_dict["num_output_files"],
                                              model_dict["sage_file"],
                                              model_dict["snapshot"])
        else:
            raise ValueError

        # The data class has read the SAGE ini file.  Update the model with the parameters
        # read and those specified by the user.
        model.update_attributes(model_dict)

        # Initialize the properties for this model. Only plotting the SMF.
        stellar_properties = ["SMF", "red_SMF", "blue_SMF"]
        min_mass = 8.0  # log10(Msun).
        max_mass = 12.0  # log10(Msun).
        bin_width = 0.1  # log10(Msun).
        bin_name = "stellar_mass_bins"
        model.init_binned_properties(min_mass, max_mass, bin_width, bin_name,
                                     stellar_properties)

        # Calculate all the properties.
        model.calc_properties_all_files(calculations_functions)

        # Append this model to the list.
        models.append(model)

Plotting Multiple Models
------------------------

All the hard work has already been done!  The plotting functions defined in
``plot_functions`` accept a list of models as their first argument.  Hence,
lets just pass in ``models``!

.. code-block:: python

    for func_name in plot_functions.keys():
        func = plot_functions[func_name][0]
        func(models, plot_output_path, plot_output_format)

This produces the stellar mass function for multiple models.  From this, we can
better understand the role that supernova feedback plays in regulating galaxy
growth.

|SMF_noSN|

.. |SMF_noSN| image:: ../figs/SMF_noSN.png

Using Keyword Arguments
-----------------------

Both the calculation and plotting functions support the use of optional keyword
arguments.  This allows finer control over plotting (e.g.) galaxy
sub-populations. The procedure for using this option for multiple models identical to the
single Model scenario which we show in :ref:`optional-kwargs`.


Analysing Multiple Models Over Redshift
---------------------------------------

To analyse multiple models over a number of redshifts, one simply needs to wrap
the code outlined :doc:`here <./history>` in the ``for model in sim_to_plot`` loop.

We defer to the `history module`_ for a full example of analysing multiple
models over redshifts.

|SMF_redshift_noSN| |SFRD_noSN| |SMD_noSN|

.. |SMF_redshift_noSN| image:: ../figs/SMF_redshift_noSN.png
.. |SFRD_noSN| image:: ../figs/SFRD_noSN.png
.. |SMD_noSN| image:: ../figs/SMD_noSN.png
.. _history module: https://github.com/sage-home/sage-model/plotting/history.py
