Analyzing Custom Properties
===========================

We show here a worked example of defining a custom property, writing a custom function to compute its value as the
galaxies are processed, and then plotting the output.  We refer to :doc:`defining_custom_properties` and
:doc:`defining_custom_functions` for further detail on the available property types that can be defined and how
calculations are handled within **sage-analysis**.

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

Worked Example
--------------

Here we compute the number of particles in the background FoF halo as a function of (sub-)halo virial mass (as a binned
property), the mass of hot gas as a function of cold gas (as a scatter property), and the time of last major merger (as
a single property) as a tracked property over redshift.

.. code-block:: python

    from sage_analysis
