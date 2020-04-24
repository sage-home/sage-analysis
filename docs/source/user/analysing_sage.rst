# Showing how to plot sub-populations.

Analysing **SAGE** Output
=========================

The output from **SAGE** is analyzed using the respective parameter files used to run **SAGE** itself. Here, we will
assume the parameter file is located in ``/home/Desktop/sage-model/input/millennium.ini``.

Basic Analysis
--------------

Out of the box, **sage-analysis** will analyse the latest snapshot (i.e., the lowest redshift) and save the plots in
the ``./plots`` directory.

.. code-block:: python

    from sage_analysis.galaxy_analysis import GalaxyAnalysis

    par_fnames = ["/home/Desktop/sage-model/input/millennium.ini"]

    galaxy_analysis = GalaxyAnalysis(par_fnames)
    galaxy_analysis.analyze_galaxies()
    galaxy_analysis.generate_plots()

Changing the Snapshot
---------------------

The

Setting Things Up
-----------------

First things first, we need to specify exactly what we want to analyse/plot,
how the plots will be saved, and where the plots will be saved.

.. code-block:: python

    # Going to just plot the stellar mass function.
    plot_toggles = {"SMF": 1}

    plot_output_format = "png"
    plot_output_path = "./plots"

    if not os.path.exists(plot_output_path):
        os.makedirs(plot_output_path)

Model Dictionary
----------------

Each model that you wish to analyse is specifed through a dictionary.  This
defines properties such as the snapshot you wish to analyse, the location of
the **SAGE** parameter file, etc.

.. code-block:: python

    millennium = { "snapshot": 63,   # Snapshot we're plotting properties at.
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

.. _func-dict:

Setting up the Calculation and Plotting Dictionaries
----------------------------------------------------

To ensure that **sage-analysis** does not perform extraneous computations, the
properties for each Model are calculated depending upon the ``plot_toggles``
specified.  For example, the black hole mass of each galaxy will only be read
if the black hole-bulge relationship plot toggle is set. We refer to :doc:`this
page <./custom_calculations` for a full list of the default plot toggles.

To achieve this, we search for all functions in a module that are named
``calc_<plot_toggle>``.  We build these functions into a dictionary that are
passed into :py:meth:`~sage_analysis.model.Model.calc_properties_all_files`.

.. code-block:: python

    from sage_analysis.utils import generate_func_dict

    # Search for functions named "calc_<plot_toggle>" in the "example_calcs"
    # module.
    calculation_functions = generate_func_dict(
                                plot_toggles,
                                module_name="sage_analysis.example_calcs",
                                function_prefix="calc"
                                )

**NOTE:** All functions must have the signature
``calc_<plot_toggle>(model, galaxies, **optional keyword arguments)``.  We
expand on this more in :ref:`optional-kwargs`.

In a similar manner, we search for all the functions in a module that are named
``plot_<plot_toggle>``.  From this dictionary, we can then iterate over and
make all the plots!

.. code-block:: python

    # Search for functions named "plot_<plot_toggles>" in the "example_plots"
    # module.
    plot_functions = generate_func_dict(
                        plot_toggles,
                        module_name="sage_analysis.example_plots",
                        function_prefix="plot_"
                        )

**NOTE:** All functions must have the signature
``calc_<plot_toggle>(list of models, plot_output_path, **optional keyword arguments)``.
We expand on this more in :ref:`optional-kwargs`.

Initializing a Model
--------------------

With the calculation functions prepped, we are now poised to perform the actual
analysis. The analysis of **SAGE** models is done through a specialized
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

Storing Galaxy Properties
-------------------------

When performing calculations, **sage-analysis** stores all the calculating
properties in the :py:attr:`~sage_analysis.model.Model.properties` attribute of the Model instance.
This attribute is a dictionary and can be used to access any of the properties
pertaining to the Model; for example, ``model.properties["SMF"]`` stores the
array representing the stellar mass function.

These properties must first be initialized. **sage-analysis** offers three ways
to compute and store galaxy properties.

Binned Properties
~~~~~~~~~~~~~~~~~

These are properties that are binned on some value.  For example: the stellar
mass function is binned depending upon the galaxy stellar mass; the fraction of
quiescent galaxies is binned upon the galaxy stellar mass; the mass of gas in
each **SAGE** reservoir (cold gas/hot gas/stars/etc) is binned upon the
friends-of-friends halo mass.  The bins themselves are
accessed through the :py:attr:`~sage_analysis.model.Model.bins` attribute of the model instance.
This attribute is a dictionary and can be used to access any of the bins for
the Model; for example, ``model.bins["stellar_mass_bins"]`` would return the
stellar mass bins used for the stellar mass function.

.. code-block:: python

    # Properties binned on stellar mass.
    stellar_properties = ["SMF", "red_SMF", "blue_SMF"]
    min_mass = 8.0  # log10(Msun).
    max_mass = 12.0  # log10(Msun).
    bin_width = 0.1  # log10(Msun).
    bin_name = "stellar_mass_bins"
    model.init_binned_properties(min_mass, max_mass, bin_width, bin_name,
                                 stellar_properties)

    # Properties binned on FoF halo mass.
    component_properties = [f"halo_{component}_fraction_sum" for component in
                            ["baryon", "stars", "cold", "hot", "ejected", "ICS", "bh"]]
    min_mass = 10.0  # log10(Msun)
    max_mass = 14.0  # log10(Msun)
    bin_width = 0.1  # log10(Msun)
    bin_name = "halo_mass_bins"
    model.init_binned_properties(min_mass, max_mass, bin_width, bin_name,
                                 component_properties)


Scatter Properties
~~~~~~~~~~~~~~~~~~

In many instances, we don't want to fit an exact line to the properties, but
rather just get a sense of the typical data point values.  For these, we want
to compute lists of ``(x, y)`` coordinates that we will plot later.  For
example, the black hole bulge relationship will show a number of black hole
masses and the corresponding bulge mass.  The (maximum) number of data points
shown on each plot can be set through the :py:attr:`~sage_analysis.model.Model.sample_size`
attribute.

.. code-block:: python

    # For each of these, we need a list for both x and y points. E.g., the
    # black hole bulge needs both "bh_mass" and "bulge_mass".
    scatter_properties = ["bh_mass", "bulge_mass", "BTF_mass", "BTF_vel"]
    model.init_scatter_properties(scatter_properties)

Single Properties
~~~~~~~~~~~~~~~~~

Finally, often we want to use a single number to summarize a property for all
galaxies across a single snapshot.  This is most useful when analyzing galaxy
properties over a range of snapshots through the `history module`_. These
properties are initialized with a value of ``0.0``.

.. code-block:: python

    single_properties = ["SMFD", "SFRD"]
    model.init_single_properties(single_properties)

Doing the Analysis and Plotting
-------------------------------

We have set up the dictionary for the plotting functions in :ref:`func-dict`.
Once all the properties have been calculated, we can finally do the plotting!

.. code-block:: python

    # Calculate all the properties.
    model.calc_properties_all_files(calculations_functions)

    # Now do the plotting.
    for func_name in plot_functions.keys():
        func = plot_functions[func_name][0]
        func([model], plot_output_path, plot_output_format)

**NOTE:** The plotting scripts accept a list of Model classes as the first
argument.  For this scenario, we only have one model and so we cast it to a
list first.

The above code snippets produce the glorious stellar mass function!

|SMF|
.. |SMF| image:: ../figs/SMF.png

.. _optional-kwargs:

Using Keyword Arguments
-----------------------

:py:func:`~sage_analysis.utils.generate_func_dict` accepts an optional
argument to allow the calculation or plotting functions to handle keyword
arugments. This argument is a dictionary with keys equal to the plot toggles.
The value of each entry is another dictionary containing all of the keyword
arguments and their corresponding value.

.. code-block:: python

    from sage_analysis.utils import generate_func_dict

    # By default, the stellar mass function is not computed for the red and blue
    # galaxy populations. Let's turn it on.
    keyword_args = {"SMF": {"calc_sub_populations": True}}

    calculation_functions = generate_func_dict(
                                plot_toggles,
                                module_name="sage_analysis.example_calcs",
                                function_prefix="calc",
                                keyword_args=keyword_args
                                )
    model.calc_properties_all_files(calculations_functions)

    # Then we can adjust "plot_SMF" to also plot these extra populations.
    keyword_args = {"SMF": {"plot_sub_populations": True}}

    plot_functions = generate_func_dict(
                        plot_toggles,
                        module_name="sage_analysis.example_plots",
                        function_prefix="plot_",
                        keyword_args=keyword_args
                        )

    # Now do the plotting with the extra kwargs.
    for func_name in plot_functions.keys():
        func = plot_functions[func_name][0]
        keyword_args = plot_functions[func_name][1]
        func(models, plot_output_path, plot_output_format, **keyword_args)

|SMF_pop|
.. |SMF_pop| image:: ../figs/SMF_pop.png

.. _repository: https://github.com/sage-home/sage-model
.. _galaxy_properties module: https://github.com/sage-home/sage-model/plotting/galaxy_properties.py
.. _history module: https://github.com/sage-home/sage-model/plotting/history.py
.. _**SAGE**: https://github.com/sage-home/sage-model
