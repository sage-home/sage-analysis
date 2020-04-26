Defining Custom Calculations
============================

# Custom calculation using a pre-defined property.

How Do Calculations Work?
-------------------------

The :py:class:`~sage_analysis.galaxy_analysis.GalaxyAnalysis` constructor accepts two key paremters:
::py:attr:`~sage_analysis.model.Model.calculation_functions` and :py:attr:`~sage_analysis.model.Model.plot_functions`.
From these two dictionaries, the exact functions that need to be run for each galaxy file and the functions that
produce the final plots are defined. Under the hood, **sage-analysis** operates by looping over
``calculation_functions`` and calling the constituent functions with the galaxies loaded for each file.  To plot, each
function in ``plot_functions`` is called using the model data that has been previously analyzed.

Hence, to define your own custom calculations, you must update the ``calculation_functions`` and ``plot_functions`` and
pass it to the :py:class:`~sage_analysis.galaxy_analysis.GalaxyAnalysis` constructor.

Building New Dictionaries
-------------------------

The :py:func:`~sage_analysis.utils.generate_func_dict` function can be used to build a new function dictionary that is
passed to the :py:class:`~sage_analysis.galaxy_analysis.GalaxyAnalysis` constructor.  This function searches the
specified module for the functions named ``function_prefix`` ``plot_toggle_key``.

For example, let's say we wish to implement our own function to plot the stellar mass function, as the default figure
does not capture all the physics we're interested in.  We first write our own function, ``my_plot_SMF``, which accepts
the arguments ``(models, snapshots, plot_output_path, plot_output_format)``, and save it in the module
``custom_functions.py``.

.. note::

   The signature of **ALL** plotting functions must be ``func(models: List[`` :py:class:`~sage_analysis.model.Model`
    ``], snapshots: List[List[int]], plot_output_path: str, plot_output_format: str)``.

   The signature of **ALL** calculations functions must be ``func(model: `` :py:class:`~sage_analysis.model.Model`
    ``, gals: Any, snapshot: int)``.

.. code-block:: python

    from sage_analysis.utils import generate_func_dict

    plot_toggles = {"SMF": True}

    plot_functions = generate_func_dict(
