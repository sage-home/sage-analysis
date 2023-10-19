Analyzing Custom Properties
===========================

We show here a worked example of defining a custom property, writing a custom function to compute its value as the
galaxies are processed, and then plotting the output.  We refer to :doc:`defining_custom_properties` for further detail
on the available property types that can be defined.

Things To Be Aware Of When Analyzing Custom Properties
------------------------------------------------------

**SAGE** operates by allowing each processor to write to its own file as galaxies are evolved through cosmic time,
with **sage-analysis** processing the galaxy properties for each of these files individually and separately.
Consequently, each property should **MUST** have an entry in ``model.properties["snapshot_<snapshot_number>]`` that is
carried **across** the different files.

For ``binned`` properties (see :doc:`defining_custom_properties`), the entry in
``model.properties["snapshot_<snapshot_number>"]`` is a list that is continuously updated for each file.  For example,
the stellar mass function for snapshot 63 is stored in ``model.properties["snapshot_63"]["SMF"]``. When the galaxies
are processed for file 0, the stellar mass function at snapshot 63 is computed and added to
``model.properties["snapshot_63"]["SMF"]``.  These galaxies are discarded and new ones read in for file 1, with the
stellar mass function at snapshot 63 for these new galaxies computed and added to
``model.properties["snapshot_63"]["SMF"]``, and so on.

For ``scatter`` properties (see :doc:`defining_custom_properties`), the entry in
``model.properties["snapshot_<snapshot_number>"]`` is an expanding list.  For example, 10 galaxies at snapshot 63 from
file 0 are appended to ``model.properties["snapshot_63"]["BTF_mass"]``, 10 galaxies from file 1, 10 galaxies from file 2, etc.

For ``single`` properties (see :doc:`defining_custom_properties`), the entry in
``model.properties["snapshot_<snapshot_number>"]`` is a single number that is adjusted for each file.  For example, the
sum of stellar mass divided by the volume at snapshot 63 in file 0 is added to
``model.properties["snapshot_63"]["SMD_history"]``.  The stellar mass density at snapshot 63 in file 1 is then added,
and so on.

Worked Examples
---------------

We show here how to compute the number of particles in the background FoF halo (as a binned property), the mass of hot
gas as a function of cold gas (as a scatter property), and the time of last major merger (as a single property) tracked
over redshift.

Number of Particles
~~~~~~~~~~~~~~~~~~~

Firstly, we need to tell **sage_analysis** the properties that we are analyzing and plotting.

.. code-block:: python

   plot_toggles = {"halo_parts": True}

Now, lets define the properties that will be used to store all of our results. As outlined in
:doc:`defining_custom_properties`, each property type is defined in a slightly different manner.

.. code-block:: python

    galaxy_properties_to_analyze = {
        "number_particles_bins": {
            "type": "binned",
            "bin_low": 0,
            "bin_high": 5,
            "bin_width": 0.1,
            "property_names": ["particle_mass_function"],
        },

Next, we need the function that will compute the values relevant for this property.

.. code-block:: python

    # Saved in ``my_calculations_functions.py``.
    from typing import Any

    import numpy as np

    from sage_analysis.model import Model

    def calc_halo_parts(model: Model, gals: Any, snapshot: int) -> None:

        non_zero_parts = np.where(gals["Len"][:] > 0)[0]
        halo_len = np.log10(gals["Len"][:][non_zero_parts])  # Ensure that the data is the same units as bins.
        gals_per_bin, _ = np.histogram(halo_len, bins=model.bins["number_particles_bins"])

        # Update properties to keep persistent across files.
        model.properties[f"snapshot_{snapshot}"]["particle_mass_function"] += gals_per_bin

Then, the function that will plot the results.

.. code-block:: python

    # Save as ``my_plot_functions.py``.
    from typing import List

    from sage_analysis.model import Model

    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np

    colors = ["r", "g", "b", "c"]
    linestyles = ["--", "-.", "."]
    markers = ["x", "o"]


    def plot_halo_parts(
        models: List[Model], snapshots: List[List[int]], plot_output_path: str, plot_output_format: str = "png",
    ) -> matplotlib.figure.Figure:

        fig = plt.figure()
        ax = fig.add_subplot(111)

        # Go through each of the models and plot.
        for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

            # Set the x-axis values to be the centre of the bins.
            bin_widths = model.bins["number_particles_bins"][1::] - model.bins["number_particles_bins"][0:-1]
            bin_middles = model.bins["number_particles_bins"][:-1] + bin_widths

            # Colour will be used for the snapshot, linestyle for the model.
            ls = linestyles[model_num]
            label = model.label

            for snapshot_num, snapshot in enumerate(model_snapshots):
                color = colors[snapshot_num]
                ax.plot(
                    bin_middles,
                    model.properties[f"snapshot_{snapshot}"]["particle_mass_function"],
                    color=color,
                    ls=ls,
                    label=f"{label} - z = {model._redshifts[snapshot]:.2f}",
                )

        ax.set_xlabel(r"$\log_{10} Number Particles in Halo$")
        ax.set_ylabel(r"$N$")

        ax.set_yscale("log", nonpositive="clip")
        ax.legend()

        fig.tight_layout()

        output_file = f"{plot_output_path}particles_in_halos.{plot_output_format}"
        fig.savefig(output_file)
        print(f"Saved file to {output_file}")
        plt.close()

        return fig

With everything defined and our functions written, we are now ready to execute **sage-analysis** itself.

.. code-block:: python

    from sage_analysis.galaxy_analysis import GalaxyAnalysis
    from sage_analysis.utils import generate_func_dict

    par_fnames = ["/home/Desktop/sage-model/input/millennium.ini"]

    # Generate the dictionaries with our custom functions.
    calculation_functions = generate_func_dict(plot_toggles, __name__, "calc_")
    plot_functions = generate_func_dict(plot_toggles, __name__, "plot_")

    # We're good to go now!
    galaxy_analysis = GalaxyAnalysis(
        par_fnames,
        plot_toggles=plot_toggles,
        galaxy_properties_to_analyze=galaxy_properties_to_analyze,
        history_redshifts=history_redshifts,
        calculation_functions=calculation_functions,
        plot_functions=plot_functions
    )

    galaxy_analysis.analyze_galaxies()
    galaxy_analysis.generate_plots()

Mass of Hot Gas as Function of Cold Gas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    plot_toggles = {"hot_cold": True}
    galaxy_properties_to_analyze = {
        "hot_cold_scatter": {
            "type": "scatter",
            "property_names": ["hot_gas", "cold_gas"],
        },
    }


    def calc_hot_cold(model: Model, gals: Any, snapshot: int) -> None:

        non_zero_stellar = np.where(gals["StellarMass"][:] > 0.0)[0]

        # Remember that mass is kept in units of 1.0e10 Msun/h. Convert to log10(Msun).
        hot_gas_mass = np.log10(gals["HotGas"][:][non_zero_stellar] * 1.0e10 / model.hubble_h)
        cold_gas_mass = np.log10(gals["ColdGas"][:][non_zero_stellar] * 1.0e10 / model.hubble_h)

        # Append to properties to keep persistent across files.
        model.properties[f"snapshot_{snapshot}"]["hot_gas"] = np.append(
            model.properties[f"snapshot_{snapshot}"]["hot_gas"], hot_gas_mass
        )

        model.properties[f"snapshot_{snapshot}"]["cold_gas"] = np.append(
            model.properties[f"snapshot_{snapshot}"]["cold_gas"], cold_gas_mass
        )


    def plot_hot_cold(
        models: List[Model], snapshots: List[List[int]], plot_output_path: str, plot_output_format: str = "png",
    ) -> matplotlib.figure.Figure:

        fig = plt.figure()
        ax = fig.add_subplot(111)

        # Go through each of the models and plot.
        for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

            # Colour will be used for the snapshot, marker style for the model.
            marker = markers[model_num]
            label = model.label

            for snapshot_num, snapshot in enumerate(model_snapshots):
                color = colors[snapshot_num]

                ax.scatter(
                    model.properties[f"snapshot_{snapshot}"]["cold_gas"],
                    model.properties[f"snapshot_{snapshot}"]["hot_gas"],
                    marker=marker,
                    s=1,
                    color=color,
                    alpha=0.5,
                    label=f"{label} - z = {model._redshifts[snapshot]:.2f}",
                )

        ax.set_xlabel(r"$\log_{10} Cold Gas Mass [M_\odot]$")
        ax.set_ylabel(r"$\log_{10} Hot Gas Mass [M_\odot]$")

        ax.legend()

        fig.tight_layout()

        output_file = f"{plot_output_path}hot_cold.{plot_output_format}"
        fig.savefig(output_file)
        print(f"Saved file to {output_file}")
        plt.close()

        return fig

Defining the Properties
~~~~~~~~~~~~~~~~~~~~~~~

Now, lets define the properties that will be used to store all of our results. As outlined in
:doc:`defining_custom_properties`, each property type is defined in a slightly different manner.

.. code-block:: python

    galaxy_properties_to_analyze = {
        "number_particles_bins": {
            "type": "binned",
            "bin_low": 0,
            "bin_high": 5,
            "bin_width": 0.1,
            "property_names": ["particle_mass_function"],
        },
        "hot_cold_scatter": {
            "type": "scatter",
            "property_names": ["hot_gas", "cold_gas"],
        },
        "time_major_merger": {
            "type": "single",
            "property_names": ["sum_time_since_major_merger", "var_time_since_major_merger", "num_galaxies"]
        }
    }

For the first property, we have used log-spaced bins. For the last property, we will track the sum and variance of the
time since major merger alongside the total number of galaxies.  Then, when we plot, we can compute the mean and show
the mean plus variance trend.

Tracking a Property Over Redshift
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We want to track the time since major merger over redshift explicitly.  To do so, we need to specify the redshifts we
wish to track it over, :py:attr:`~sage_analysis.galaxy_analysis.GalaxyAnalysis.history_redshifts`.

.. code-block:: python

    history_redshifts = {"major_merger_history": "All"}

.. note::

    The key names in this dictionary must exactly match the key name in ``plot_toggles``.

Defining the Functions
~~~~~~~~~~~~~~~~~~~~~~

The :py:class:`~sage_analysis.galaxy_analysis.GalaxyAnalysis` constructor accepts two key parameters:
::py:attr:`~sage_analysis.model.Model.calculation_functions` and :py:attr:`~sage_analysis.model.Model.plot_functions`.
From these two dictionaries, the exact functions that need to be run for each galaxy file and the functions that
produce the final plots are defined. Under the hood, **sage-analysis** operates by looping over
``calculation_functions`` and calling the constituent functions with the galaxies loaded for each file.  To plot, each
function in ``plot_functions`` is called using the model data that has been previously analyzed.

Hence, to define your own custom properties, we must first update the ``calculation_functions`` and ``plot_functions``
and pass it to the :py:class:`~sage_analysis.galaxy_analysis.GalaxyAnalysis` constructor.

Let's write the functions that will define the calculation functions that will be saved to module
``my_calculation_functions.py``.  These will use our properties defined above to keep the values across different
files.

.. code-block:: python

    # Saved in ``my_calculations_functions.py``.
    from typing import Any

    import numpy as np

    from sage_analysis.model import Model

    def calc_halo_parts(model: Model, gals: Any, snapshot: int) -> None:

        non_zero_parts = np.where(gals["Len"][:] > 0)[0]
        halo_len = np.log10(gals["Len"][:][non_zero_parts])  # Ensure that the data is the same units as bins.
        gals_per_bin, _ = np.histogram(halo_len, bins=model.bins["number_particles_bins"])

        # Update properties to keep persistent across files.
        model.properties[f"snapshot_{snapshot}"]["particle_mass_function"] += gals_per_bin


    def calc_hot_cold(model: Model, gals: Any, snapshot: int) -> None:

        non_zero_stellar = np.where(gals["StellarMass"][:] > 0.0)[0]

        # Remember that mass is kept in units of 1.0e10 Msun/h. Convert to log10(Msun).
        hot_gas_mass = np.log10(gals["HotGas"][:][non_zero_stellar] * 1.0e10 / model.hubble_h)
        cold_gas_mass = np.log10(gals["ColdGas"][:][non_zero_stellar] * 1.0e10 / model.hubble_h)

        # Append to properties to keep persistent across files.
        model.properties[f"snapshot_{snapshot}"]["hot_gas"] = np.append(
            model.properties[f"snapshot_{snapshot}"]["hot_gas"], hot_gas_mass
        )

        model.properties[f"snapshot_{snapshot}"]["cold_gas"] = np.append(
            model.properties[f"snapshot_{snapshot}"]["cold_gas"], cold_gas_mass
        )


    def calc_major_merger_history(model: Model, gals: Any, snapshot: int) -> None:

        non_zero_stellar = np.where(gals["StellarMass"][:] > 0.0)[0]

        time_since_major_merger = gals["TimeOfLastMajorMerger"][:][non_zero_stellar]

        # A galaxy that has not experienced a major merger will have a value of -1. Lets filter these out.
        time_since_major_merger = time_since_major_merger[time_since_major_merger > 0.0]

        # We will handle dividing out the number of galaxies and the number of samples (i.e., number of files) when it
        # comes time to plot.
        model.properties[f"snapshot_{snapshot}"]["sum_time_since_major_merger"] += np.sum(time_since_major_merger)
        model.properties[f"snapshot_{snapshot}"]["var_time_since_major_merger"] += np.var(time_since_major_merger)
        model.properties[f"snapshot_{snapshot}"]["num_galaxies"] += len(time_since_major_merger)

With our calculation functions defined, we now need to define the plot functions.  These functions will be used by
**sage-analysis** to generate the plots themselves.  We will save these functions to the module
``my_plot_functions.py``.

.. code-block:: python

    # Save as ``my_plot_functions.py``.
    from typing import List

    from sage_analysis.model import Model

    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np

    colors = ["r", "g", "b", "c"]
    linestyles = ["--", "-.", "."]
    markers = ["x", "o"]


    def plot_halo_parts(
        models: List[Model], snapshots: List[List[int]], plot_output_path: str, plot_output_format: str = "png",
    ) -> matplotlib.figure.Figure:

        fig = plt.figure()
        ax = fig.add_subplot(111)

        # Go through each of the models and plot.
        for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

            # Set the x-axis values to be the centre of the bins.
            bin_widths = model.bins["number_particles_bins"][1::] - model.bins["number_particles_bins"][0:-1]
            bin_middles = model.bins["number_particles_bins"][:-1] + bin_widths

            # Colour will be used for the snapshot, linestyle for the model.
            ls = linestyles[model_num]
            label = model.label

            for snapshot_num, snapshot in enumerate(model_snapshots):
                color = colors[snapshot_num]
                ax.plot(
                    bin_middles,
                    model.properties[f"snapshot_{snapshot}"]["particle_mass_function"],
                    color=color,
                    ls=ls,
                    label=f"{label} - z = {model._redshifts[snapshot]:.2f}",
                )

        ax.set_xlabel(r"$\log_{10} Number Particles in Halo$")
        ax.set_ylabel(r"$N$")

        ax.set_yscale("log", nonpositive="clip")
        ax.legend()

        fig.tight_layout()

        output_file = f"{plot_output_path}particles_in_halos.{plot_output_format}"
        fig.savefig(output_file)
        print(f"Saved file to {output_file}")
        plt.close()

        return fig


    def plot_hot_cold(
        models: List[Model], snapshots: List[List[int]], plot_output_path: str, plot_output_format: str = "png",
    ) -> matplotlib.figure.Figure:

        fig = plt.figure()
        ax = fig.add_subplot(111)

        # Go through each of the models and plot.
        for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

            # Colour will be used for the snapshot, marker style for the model.
            marker = markers[model_num]
            label = model.label

            for snapshot_num, snapshot in enumerate(model_snapshots):
                color = colors[snapshot_num]

                ax.scatter(
                    model.properties[f"snapshot_{snapshot}"]["cold_gas"],
                    model.properties[f"snapshot_{snapshot}"]["hot_gas"],
                    marker=marker,
                    s=1,
                    color=color,
                    alpha=0.5,
                    label=f"{label} - z = {model._redshifts[snapshot]:.2f}",
                )

        ax.set_xlabel(r"$\log_{10} Cold Gas Mass [M_\odot]$")
        ax.set_ylabel(r"$\log_{10} Hot Gas Mass [M_\odot]$")

        ax.legend()

        fig.tight_layout()

        output_file = f"{plot_output_path}hot_cold.{plot_output_format}"
        fig.savefig(output_file)
        print(f"Saved file to {output_file}")
        plt.close()

        return fig


    def plot_major_merger_history(
        models: List[Model], snapshots: List[List[int]], plot_output_path: str, plot_output_format: str = "png",
    ) -> matplotlib.figure.Figure:

        fig = plt.figure()
        ax = fig.add_subplot(111)

        for (model_num, model) in enumerate(models):

            label = model.label
            color = colors[model_num]
            linestyle = linestyles[model_num]
            marker = markers[model_num]

            sum_time_since_major_merger = np.array(
                [model.properties[f"snapshot_{snap}"]["sum_time_since_major_merger"] for snap in range(len(model.redshifts))]
            )
            var_time_since_major_merger = np.array(
                [model.properties[f"snapshot_{snap}"]["var_time_since_major_merger"] for snap in range(len(model.redshifts))]
            )
            num_galaxies = np.array(
                [model.properties[f"snapshot_{snap}"]["num_galaxies"] for snap in range(len(model.redshifts))]
            )
            redshifts = model.redshifts

            # mean =  sum / number of samples.
            mean_time_since_major_merger = sum_time_since_major_merger / num_galaxies

            # Need to divide out the number of samples for the variance. This is the number of files that we analyzed.
            var_time_since_major_merger /= (model.last_file_to_analyze - model.first_file_to_analyze + 1)

            # All snapshots are initialized with zero values, we only want to plot those non-zero values.
            non_zero_inds = np.where(mean_time_since_major_merger > 0.0)[0]

            # Only use a line if we have enough snapshots to plot.
            if len(non_zero_inds) > 20:
                ax.plot(
                    redshifts[non_zero_inds],
                    mean_time_since_major_merger[non_zero_inds],
                    label=label,
                    color=color,
                    ls=linestyle
                )
            else:
                ax.scatter(
                    redshifts[non_zero_inds],
                    mean_time_since_major_merger[non_zero_inds],
                    label=label,
                    color=color,
                    marker=marker,
                )

        ax.set_xlabel(r"$\mathrm{redshift}$")
        ax.set_ylabel(r"$Time Since Last Major Merger [Myr]$")

        ax.set_xlim([0.0, 8.0])
        #ax.set_ylim([-3.0, -0.4])

        ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
        #ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))

        ax.legend()

        fig.tight_layout()

        output_file = f"{plot_output_path}time_since_last_major_merger.{plot_output_format}"
        fig.savefig(output_file)
        print("Saved file to {0}".format(output_file))
        plt.close()

        return fig

Putting it Together
~~~~~~~~~~~~~~~~~~~

With everything defined and our functions written, we are now ready to execute **sage-analysis** itself.

.. code-block:: python

    import my_calculation_functions, my_plot_functions

    from sage_analysis.galaxy_analysis import GalaxyAnalysis
    from sage_analysis.utils import generate_func_dict

    par_fnames = ["/home/Desktop/sage-model/input/millennium.ini"]

    # Generate the dictionaries with our custom functions.
    calculation_functions = generate_func_dict(plot_toggles, "my_calculation_functions", "calc_")
    plot_functions = generate_func_dict(plot_toggles, "my_plot_functions", "plot_")

    # We're good to go now!
    galaxy_analysis = GalaxyAnalysis(
        par_fnames,
        plot_toggles=plot_toggles,
        galaxy_properties_to_analyze=galaxy_properties_to_analyze,
        history_redshifts=history_redshifts,
        calculation_functions=calculation_functions,
        plot_functions=plot_functions
    )

    galaxy_analysis.analyze_galaxies()
    galaxy_analysis.generate_plots()

And these are our plots that are generated...
