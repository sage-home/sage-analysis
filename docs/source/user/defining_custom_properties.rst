Defining Custom Properties
==========================

Default Properties
------------------

Out of the box, **sage-analysis** supports the analysis of a number of different properties.

+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------+
| Property                     | Plot Toggle Name | Description                                                                                                                          | Property Type                     |
+==============================+==================+======================================================================================================================================+===================================+
| Stellar mass function        | SMF              | Number of galaxies with a given stellar mass.                                                                                        | Binned (on stellar mass).         |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------+
| Baryonic mass function       | BMF              | Number of galaxies with a given stellar plus cold gas mass.                                                                          | Binned (on stellar mass).         |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------+
| Gas mass function            | GMF              | Number of galaxies with a given cold gas mass.                                                                                       | Binned (on stellar mass).         |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------+
| Baryonic Tully-Fisher        | BTF              | Maximum velocity of a galaxy as a function of baryonic (stellar plus cold gas) mass.                                                 | Scatter.                          |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------+
| Specific star formation rate | sSFR             | Specific star formation rate as a function of stellar mass.                                                                          | Scatter.                          |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------+
| Gas fraction                 | gas_frac         | Fraction of baryons (stellar plus cold gas) in the form of cold gas as a function of stellar mass.                                   | Scatter.                          |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------+
| Mass metallicity             | metallicity      | Metallicity as a function of stellar mass.                                                                                           | Scatter.                          |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------+
| Black hole bulge             | bh_bulge         | Mass of galaxy black hole as a function of galaxy bulge mass.                                                                        | Scatter.                          |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------+
| Quiescent galaxy population  | quiescent        | Fraction of galaxies that are quiescent as a function of stellar mass.                                                               | Binned (on stellar mass).         |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------+
| Bulge fraction               | bulge_fraction   | Fraction of stellar mass in the form of bulge/disk as a function of stellar mass.                                                    | Scatter.                          |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------+
| Baryon fraction              | baryon_fraction  | Baryon fraction in each reservoir (cold, hot, stellar, ejected, intracluster, and black hole) as a function of FoF halo virial mass. | Binned (on FoF halo virial mass). |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------+
| Reservoir mass               | reservoirs       | Amount of mass in each reservoir (cold, hot, stellar, ejected, intracluster, and black hole) as a function of FoF halo virial mass.  | Scatter.                          |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------+
| Spatial distribution         | spatial          | Spatial distribution of galaxies across the simulation box.                                                                          | Scatter.                          |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------+

There are also a handful of toggles available to analyse properties over a number of redshifts.

+-----------------------------+------------------+-------------------------------------------------------------------------------------------+---------------------------+
| Property                    | Plot Toggle Name | Description                                                                               | Binning Type              |
+=============================+==================+===========================================================================================+===========================+
| Stellar mass function       | SMF_history      | Number of galaxies with a given stellar mass over multiple redshifts for each model.      | Binned (on stellar mass). |
+-----------------------------+------------------+-------------------------------------------------------------------------------------------+---------------------------+
| Star formation rate density | SFRD_history     | Total star formation rate density across entire simulation box as a function of redshift. | Single.                   |
+-----------------------------+------------------+-------------------------------------------------------------------------------------------+---------------------------+
| Stellar mass density        | SMD_history      | Total stellar mass density across entire simulation box as a function of redshift.        | Single.                   |
+-----------------------------+------------------+-------------------------------------------------------------------------------------------+---------------------------+

Property Types
--------------

**sage-analysis** supports three different property types.

binned
~~~~~~

These properties are binned against another variable.  For example, the stellar mass function counts the number of
galaxy in stellar mass bins, the baryon fraction measures the fraction of baryons in each reservoir in
friends-of-friend halo virial mass.

A property of this types requires the following fields:

.. code-block:: python

    <string denoting the name of the bins>: {
        "type": "binned",
        "bin_low": <float denoting the lower bound of the bins>,
        "bin_high": <float denoting the upper bound of the bins>,
        "bin_width": <float denoting the width of each bin>,
        "property_names": [<list of strings denoting the name of properties to be initialized>],
    }

For example, the stellar mass bins needed for default operation are initialized using:

.. code-block:: python

    "stellar_mass_bins": {
        "type": "binned",
        "bin_low": 8.0,
        "bin_high": 12.0,
        "bin_width": 0.1,
        "property_names": [
            "SMF", "red_SMF", "blue_SMF", "BMF", "GMF",
            "centrals_MF", "satellites_MF", "quiescent_galaxy_counts",
            "quiescent_centrals_counts", "quiescent_satellites_counts",
            "fraction_bulge_sum", "fraction_bulge_var",
            "fraction_disk_sum", "fraction_disk_var", "SMF_history",
        ],
    }

The bins are accessed using ``model.properties["<bin_name>"]`` (e.g., ``model.properties["stellar_mass_bins"]``) and
the properties themselves as ``model.properties["<propety_name>"]`` (e.g., ``model.properties["SMF"]`` or
``model.properties["quiescent_galaxy_counts"]``.  Each property is initialized as a list of 0s.

scatter
~~~~~~~

To get a better picutre of some properties, it is useful to display them as a scatter plot.  For example, the
``metallicity`` property shows the stellar mass vs metallicity for a number of randomly selected galaxies.

A property of this types requires the following fields:

.. code-block:: python

    <string denoting a unique name>: {
        "type": "scatter",
        "property_names": [<list of strings denoting the name of properties to be initialized>],
    }

For example, the default scatter properties are initialized using:

.. code-block:: python

    "scatter_properties": {
        "type": "scatter",
        "property_names": [
            "BTF_mass", "BTF_vel", "sSFR_mass", "sSFR_sSFR",
            "gas_frac_mass", "gas_frac", "metallicity_mass",
            "metallicity", "bh_mass", "bulge_mass", "reservoir_mvir",
            "reservoir_stars", "reservoir_cold", "reservoir_hot",
            "reservoir_ejected", "reservoir_ICS", "x_pos",
            "y_pos", "z_pos"
        ],
    }

The properties are accessed as ``model.properties["<propety_name>"]`` (e.g., ``model.properties["BTF_mass"]`` or
``model.properties["BTF_vel"]``.  Each property is initialized as an empty list.

single
~~~~~~

Finally, we may wish to summarize a property using a single number over an entire snapshot.  For example, the stellar
mass density is the sum of stellar mass divided by the volume for a single snapshot.  This is useful for tracking
properties over a number of snapshots as they can then be depicted as a line on a stellar mass density vs redshift
plot.

A property of this types requires the following fields:

.. code-block:: python

    <string denoting a unique name>: {
        "type": "single",
        "property_names": [<list of strings denoting the name of properties to be initialized>],
    }

For example, the default single properties are initialized using:

.. code-block:: python

    "scatter_properties": {
        "type": "single",
        "property_names": ["SMD_history", "SFRD_history"],
    }

The properties are accessed as ``model.properties["<propety_name>"]`` (e.g., ``model.properties["SMD_history"]`` or
``model.properties["SFRD_history"]``.  Each property is initialized with a value of ``0.0``.
