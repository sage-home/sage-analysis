Analyzing **SAGE** Output
=========================

The output from **SAGE** is analyzed using the respective parameter files used to run **SAGE** itself. Here, we will
assume the parameter file is located in ``/home/Desktop/sage-model/input/millennium.ini``.

On this page, we outline some of the basic features of **sage-analysis** that can be used to analyze and plot **SAGE**
output.

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

The output path can be changed by specifying ``plot_output_path`` as
``galaxy_analysis.generate_plots(plot_output_path="<new>/<output>/<path>/")``.

If you ran **SAGE** using ``sage-binary`` output, you will need to specify the
:py:attr:`~sage_analysis.model.Model.first_file_to_analyze`, :py:attr:`~sage_analysis.model.Model.last_file_to_analyze`, and
:py:attr:`~sage_analysis.model.Model.num_sage_output_files` for each model.  These will need to be specified for all
the following examples. For brevity, we will omit them in the following and assume that **SAGE** has been run using
``sage-hdf5`` output (recommended).

.. code-block:: python

    from sage_analysis.galaxy_analysis import GalaxyAnalysis

    par_fnames = ["/home/Desktop/sage-model/input/millennium.ini"]
    first_files_to_analyze = [0]  # The first files that you wish to analyze + plot.
    last_files_to_analyze = [0]  # The last files that you wish to analyze + plot.
    num_sage_output_files = [1]  # The number of files that SAGE produced; usually the number of processors it ran on.

    galaxy_analysis = GalaxyAnalysis(
        par_fnames,
        first_files_to_analyze=first_files_to_analyze,
        last_files_to_analyze=last_files_to_analyze,
        num_sage_output_files=num_sage_output_files
    )
    galaxy_analysis.analyze_galaxies()
    galaxy_analysis.generate_plots()

Analyze Only a Subset of Files
------------------------------

For extremey large simulations, it may be prudent to analyze only a subset of files. For example, if **SAGE** run in
parallel across 32 processors, we may only wish to analyze a quarter of these. This can be achieved by specifying the
:py:attr:`~sage_analysis.model.Model.first_file_to_analyze` and
:py:attr:`~sage_analysis.model.Model.last_file_to_analyze` for each model.

.. code-block:: python

    from sage_analysis.galaxy_analysis import GalaxyAnalysis

    par_fnames = ["/home/Desktop/sage-model/input/millennium.ini"]
    first_files_to_analyze = [0]
    last_files_to_analyze = [7]

    galaxy_analysis = GalaxyAnalysis(
        par_fnames,
        first_files_to_analyze=first_files_to_analyze,
        last_files_to_analyze=last_files_to_analyze,
    )
    galaxy_analysis.analyze_galaxies()
    galaxy_analysis.generate_plots()

Turning On and Off Properties
-----------------------------

Properties are analyzed and plotted according to the values in
:py:attr:`~sage_analysis.galaxy_analysis.GalaxyAnalysis.plot_toggles`. The default values of this dictionary are set to
analyze all basic properties, with the exception of properties tracked over time.

.. code-block:: python

    from sage_analaysis.default_analysis_arguments import default_plot_toggles
    print(default_plot_toggles)

    >>> {
            'SMF': True,
            'BMF': True,
            'GMF': True,
            'BTF': True,
            'sSFR': True,
            'gas_fraction': True,
            'metallicity': True,
            'bh_bulge': True,
            'quiescent': True,
            'bulge_fraction': True,
            'baryon_fraction': True,
            'reservoirs': True,
            'spatial': True,
            'SMF_history': False,
            'SFRD_history': False,
            'SMD_history': False
        }

By adjusting these properties, or specifying a custom set, you can control which properties you want to analyze.

.. code-block:: python

    from sage_analaysis.default_analysis_arguments import default_plot_toggles
    from sage_analysis.galaxy_analysis import GalaxyAnalysis

    par_fnames = ["/home/Desktop/sage-model/input/millennium.ini"]

    # Plot only the stellar mass function and black hole-bulge relationship.
    plot_toggles = {"SMF": True, "bh_bulge": True}
    galaxy_analysis = GalaxyAnalysis(par_fnames, plot_toggles=plot_toggles)
    galaxy_analysis.analyze_galaxies()
    galaxy_analysis.generate_plots()

    # Plot all properties EXCEPT the mass-metallicity relationship.
    plot_toggles = default_plot_toggles.copy()  # Copy to ensure ``default_plot_toggles`` aren't overwritten.
    plot_toggles["metallicity"] = False

    galaxy_analysis = GalaxyAnalysis(par_fnames, plot_toggles=plot_toggles)
    galaxy_analysis.analyze_galaxies()
    galaxy_analysis.generate_plots()

Analyzing Basic Properties Over Redshift
----------------------------------------

It can also be very useful to investigate how properties evolve over many snapshots.  By default, **sage-analysis**
supports analyzing the stellar mass function, stellar mass density, and star formation rate density.

.. note::
    Ensure that **SAGE** has outputs for multiple snapshots.  Try setting ``NumOutputs`` to ``-1`` and re-running
    **SAGE**.

These extra properties can be set by turning their respective entries in ``plot_toggles``.

.. code-block:: python

    from sage_analysis.galaxy_analysis import GalaxyAnalysis

    par_fnames = ["/home/Desktop/sage-model/input/millennium.ini"]
    plot_toggles = {"SMF_history": True, "SMD_history": True, "SFRD_history": True}

    galaxy_analysis = GalaxyAnalysis(par_fnames, plot_toggles=plot_toggles)
    galaxy_analysis.analyze_galaxies()
    galaxy_analysis.generate_plots()

By default, these extra properties are analyzed and plotted for all available redshifts.  You can also specify which
redshifts you want to analyze, with **sage-analysis** selecting the snapshots that are closest to the desired redshifts
specified.  This is especially useful for the stellar mass function where we often want to investigate the evolution at
specific redshifts.

.. code-block:: python

    from sage_analysis.galaxy_analysis import GalaxyAnalysis

    par_fnames = ["/home/Desktop/sage-model/input/millennium.ini"]
    plot_toggles = {"SMF_history": True}
    history_redshifts = {"SMF_history": [0.0, 0.5, 1.0, 2.0, 3.0]}

    galaxy_analysis = GalaxyAnalysis(par_fnames, plot_toggles=plot_toggles, history_redshifts=history_redshifts)
    galaxy_analysis.analyze_galaxies()
    galaxy_analysis.generate_plots()

To analyse and plot these properties in addition to the other properties (e.g., the baryon fraction, quiescent
fraction, etc), use and update the ``default_plot_toggles`` value.

.. code-block:: python

    from sage_analaysis.default_analysis_arguments import default_plot_toggles

    plot_toggles = default_plot_toggles.copy()  # Copy to ensure ``default_plot_toggles`` aren't overwritten.

    plot_toggles["SMF_history"] = True
    plot_toggles["SMD_history"] = True
    plot_toggles["SFRD_history"] = True

    galaxy_analysis = GalaxyAnalysis(par_fnames, plot_toggles=plot_toggles)
    galaxy_analysis.analyze_galaxies()
    galaxy_analysis.generate_plots()

Changing the Snapshot
---------------------

By default, **sage-analysis** will analyze the lowest redshift snapshot for each model.  This behaviour can be adjusted
to analyze any arbitrary snapshot.

.. code-block:: python

    from sage_analysis.galaxy_analysis import GalaxyAnalysis

    par_fnames = ["/home/Desktop/sage-model/input/millennium.ini"]
    snapshots = [[50]]

    galaxy_analysis = GalaxyAnalysis(par_fnames)
    galaxy_analysis.analyze_galaxies(snapshots=snapshots)
    galaxy_analysis.generate_plots(snapshots=snapshots)

Changing the Redshift
---------------------

Alternatively, rather than specifying the snapshot for each model, one can specify the redshift.  **sage-analysis**
will analyze the snapshot closest to these redshifts.

.. code-block:: python

    from sage_analysis.galaxy_analysis import GalaxyAnalysis

    par_fnames = ["/home/Desktop/sage-model/input/millennium.ini"]
    redshifts = [[1.0]]

    galaxy_analysis = GalaxyAnalysis(par_fnames)
    galaxy_analysis.analyze_galaxies(redshifts=redshifts)
    galaxy_analysis.generate_plots(redshifts=redshifts)

.. note::
   The ``snapshots`` and ``redshifts`` parameters **cannot both be** specified, only one may be used.

Multiple Models
---------------

**sage-analysis** supports analyzing and plotting of multiple **SAGE** model outputs.  For example, let's say we want
to compare the stellar mass function for **SAGE** run with and without supernovae feedback.  This model has been run
using a parameter file ``/home/Desktop/sage-model/input/millennium_no_SN.ini``

.. code-block:: python

    from sage_analysis.galaxy_analysis import GalaxyAnalysis

    par_fnames = ["/home/Desktop/sage-model/input/millennium.ini", "/home/Desktop/sage-model/input/millennium_no_SN.ini"]
    labels = ["Supernovae feedback on", "Supernovae feedback off"]

    galaxy_analysis = GalaxyAnalysis(par_fnames, labels=labels)
    galaxy_analysis.analyze_galaxies()
    galaxy_analysis.generate_plots()

Multiple Simulations
--------------------

In the above example, we ran **SAGE** on the same underlying N-body simulation.  However, we may wish to analyze how
**SAGE** performs on different simulations, at the same redshift; e.g., we may wish to compare the stellar mass
function at z = 1 for *Millennium* and *Bolshoi*.

.. code-block:: python

    from sage_analysis.galaxy_analysis import GalaxyAnalysis

    par_fnames = ["/home/Desktop/sage-model/input/millennium.ini", "/home/Desktop/sage-model/input/bolshoi.ini"]
    labels = ["Millennium", "Bolshoi"]

    galaxy_analysis = GalaxyAnalysis(par_fnames, labels=labels)

    redshifts = [[1.0], [1.0]]  # Specify the redshift for each model.
    galaxy_analysis.analyze_galaxies(redshifts=redshifts)
    galaxy_analysis.generate_plots(redshifts=redshifts)

Or perhaps we wish to see how the stellar mass density evolves for the different simulations...

.. code-block:: python

    from sage_analysis.galaxy_analysis import GalaxyAnalysis

    par_fnames = ["/home/Desktop/sage-model/input/millennium.ini", "/home/Desktop/sage-model/input/bolshoi.ini"]
    labels = ["Millennium", "Bolshoi"]
    plot_toggles = {"SFRD_history": True}

    galaxy_analysis = GalaxyAnalysis(par_fnames, plot_toggles=plot_toggles)

    galaxy_analysis.analyze_galaxies()
    galaxy_analysis.generate_plots()


Adding Extra Keywords for Analysis and Plotting
-----------------------------------------------

Some properties can be broken down into sub-populations and analyzed separately. For example, the stellar mass function
can be split into red and blue galaxies or the baryon fraction can be split into its constituent reservoirs.  To access
these extra functionalities, the :py:attr:`~sage_analysis.model.Model.calculation_functions` and
:py:attr:`~sage_analysis.model.Model.plot_functions` dictionaries passed to the
:py:class:`~sage_analysis.galaxy_analysis.GalaxyAnalysis` constructor need to be adjusted.

.. code-block:: python

    from sage_analysis.utils import generate_func_dict
    from sage_analysis.galaxy_analysis import GalaxyAnalysis

    par_fnames = ["/home/Desktop/sage-model/input/millennium.ini"]
    plot_toggles = {"SMF": True, "baryon_fraction": True}

    # For each toggle, specify the extra keyword arguments and their values.

    # The calculation and plotting step can each have different keywords.
    extra_keywords_calculations = {"SMF": {"calc_sub_populations": True}}
    extra_keywords_plotting = {
        "SMF": {"plot_sub_populations": True},
        "baryon_fraction": {"plot_sub_populations": True}
    }

    # Now build a dictionary with these extra arguments.
    calculation_functions = generate_func_dict(
        plot_toggles, "sage_analysis.example_calcs", "calc_", extra_keywords_calculations
    )
    plot_functions = generate_func_dict(
        plot_toggles, "sage_analysis.example_plots", "plot_", extra_keywords_plotting
    )

    # Then construct with these new dictionaries.
    galaxy_analysis = GalaxyAnalysis(
        par_fnames,
        plot_toggles=plot_toggles,
        calculation_functions=calculation_functions,
        plot_functions=plot_functions
    )

    galaxy_analysis.analyze_galaxies()
    galaxy_analysis.generate_plots()

.. _**SAGE**: https://github.com/sage-home/sage-model
