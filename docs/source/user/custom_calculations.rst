Defining Custom Properties
==========================

We have previously given examples how to analyse **SAGE** properties at a :doc:`single
snapshot <./analysing_sage>` or :doc:`over multiple redshifts <./history>`.

Default Properties
------------------

The **sage-analysis** package includes a number of properties that can be analysed
and plotted by default.

+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Property                     | Plot Toggle Name | Description                                                                                                                          | Property Type |
+==============================+==================+======================================================================================================================================+===============+
| Stellar mass function        | SMF              | Number of galaxies with a given stellar mass.                                                                                        | Binned.       |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Baryonic mass function       | BMF              | Number of galaxies with a given stellar plus cold gas mass.                                                                          | Binned        |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Gas mass function            | GMF              | Number of galaxies with a given cold gas mass.                                                                                       | Binned.       |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Baryonic Tully-Fisher        | BTF              | Maximum velocity of a galaxy as a function of baryonic (stellar plus cold gas) mass.                                                 | Scatter.      |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Specific star formation rate | sSFR             | Specific star formation rate as a function of stellar mass.                                                                          | Scatter.      |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Gas fraction                 | gas_frac         | Fraction of baryons (stellar plus cold gas) in the form of cold gas as a function of stellar mass.                                   | Scatter.      |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Mass metallicity             | metallicity      | Metallicity as a function of stellar mass.                                                                                           | Scatter.      |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Black hole bulge             | bh_bulge         | Mass of galaxy black hole as a function of galaxy bulge mass.                                                                        | Scatter.      |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Quiescent galaxy population  | quiescent        | Fraction of galaxies that are quiescent as a function of stellar mass.                                                               | Binned.       |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Bulge fraction               | bulge_fraction   | Fraction of stellar mass in the form of bulge/disk as a function of stellar mass.                                                    | Scatter.      |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Baryon fraction              | baryon_fraction  | Baryon fraction in each reservoir (cold, hot, stellar, ejected, intracluster, and black hole) as a function of FoF halo virial mass. | Binned.       |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Reservoir mass               | reservoirs       | Amount of mass in each reservoir (cold, hot, stellar, ejected, intracluster, and black hole) as a function of FoF halo virial mass.  | Scatter.      |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+---------------+
| Spatial distribution         | spatial          | Spatial distribution of galaxies across the simulation box.                                                                          | Scatter.      |
+------------------------------+------------------+--------------------------------------------------------------------------------------------------------------------------------------+---------------+

There are also a handful of toggles available to analyse properties over a
number of redshifts.

+-----------------------------+------------------+-------------------------------------------------------------------------------------------+--------------+
| Property                    | Plot Toggle Name | Description                                                                               | Binning Type |
+=============================+==================+===========================================================================================+==============+
| Stellar mass function       | SMF_z            | Number of galaxies with a given stellar mass over multiple redshifts for each model.      | Binned.      |
+-----------------------------+------------------+-------------------------------------------------------------------------------------------+--------------+
| Star formation rate density | SFRD_z           | Total star formation rate density across entire simulation box as a function of redshift. | Single.      |
+-----------------------------+------------------+-------------------------------------------------------------------------------------------+--------------+
| Stellar mass density        | SMD_z            | Total stellar mass density across entire simulation box as a function of redshift.        | Single.      |
+-----------------------------+------------------+-------------------------------------------------------------------------------------------+--------------+

Adding Your Own Properties
--------------------------

**SAGE** operates by allowing each processor to write to its own file as
galaxies are evolved through cosmic time.  **sage-analysis** processes galaxies
properties on each of these files individually.  After calculating each
property, they are stored in the
:py:attr:`~sage_analysis.model.Model.properties` attribute and carried across
files.  The pseudo-code looks like this:

.. code-block:: python

    for file in num_files:

        compute stellar mass function for file
        add stellar mass function to Model.properties["SMF"] array.

        calculate black hole bulge relationship for file
        extend the Model.properties["bh_mass"] and Model.properties["bulge_mass"] lists

        ...complete for other properties...


