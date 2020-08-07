"""
Here we show a myriad of functions that can be used to calculate properties from the
**SAGE** output.  By setting the correct plot toggles and calling
:py:func:`~sage_analysis.utils.generate_func_dict`, a dictionary containing these
functions can be generated and passed to
:py:meth:`~sage_analysis.model.Model.calc_properties_all_files` to calculate the
properties.

The properties are stored (and updated) in the
:py:attr:`~sage_analysis.model.Model.properties` attribute.

We refer to :doc:`../user/analysing_sage` for more information on how the calculations are
handled.

Author: Jacob Seiler
"""

import warnings

import numpy as np
from scipy import stats

from sage_analysis.model import Model


def calc_SMF(
    model: Model,
    gals,
    snapshot: int,
    calc_sub_populations: bool = False,
    smf_property_name: str = "SMF"
):
    """
    Calculates the stellar mass function of the given galaxies.  That is, the number of galaxies at a given stellar
    mass.

    The ``Model.properties["snapshot_<snapshot>"]"SMF"]`` array will be updated. We also split the galaxy population
    into "red" and "blue" based on the value of :py:attr:`~sage_analysis.model.Model.sSFRcut` and update the
    ``Model.properties["snapshot_<snapshot>"]["red_SMF"]`` and ``Model.properties["snapshot_<snapshot>"]["blue_SMF"]``
    arrays.

    Parameters
    ----------
    snapshot : int
        The snapshot the SMF is being calculated at.

    plot_sub_populations : boolean, optional
        If ``True``, calculates the stellar mass function for red and blue sub-populations.

    smf_property_name : string, optional
        The name of the property used to store the stellar mass function. Useful if different calculations are
        computing the stellar mass function but saving it as a different property.
    """

    non_zero_stellar = np.where(gals["StellarMass"][:] > 0.0)[0]

    if len(non_zero_stellar) == 0:
        logger.info(f"Could not find any galaxies with non-zero stellar mass for the stellar mass function.")
        return

    stellar_mass = np.log10(gals["StellarMass"][:][non_zero_stellar] * 1.0e10 / model.hubble_h)

    gals_per_bin, _ = np.histogram(stellar_mass, bins=model.bins["stellar_mass_bins"])
    model.properties[f"snapshot_{snapshot}"][f"{smf_property_name}"] += gals_per_bin

    # We often want to plot the red and blue subpopulations. So bin them if requested.
    if calc_sub_populations:

        sSFR = (gals["SfrDisk"][:][non_zero_stellar] + gals["SfrBulge"][:][non_zero_stellar]) / \
            (gals["StellarMass"][:][non_zero_stellar] * 1.0e10 / model.hubble_h)
        red_gals = np.where(sSFR < 10.0**model._sSFRcut)[0]
        red_mass = stellar_mass[red_gals]
        counts, _ = np.histogram(red_mass, bins=model.bins["stellar_mass_bins"])
        model.properties[f"snapshot_{snapshot}"]["red_SMF"] += counts

        blue_gals = np.where(sSFR > 10.0**model._sSFRcut)[0]
        blue_mass = stellar_mass[blue_gals]
        counts, _ = np.histogram(blue_mass, bins=model.bins["stellar_mass_bins"])
        model.properties[f"snapshot_{snapshot}"]["blue_SMF"] += counts


def calc_BMF(model, gals, snapshot: int):
    """
    Calculates the baryon mass function of the given galaxies.  That is, the number of galaxies at a given baryon
    (stellar + cold gas) mass.

    The ``Model.properties["snapshot_<snapshot>"]["BMF"]`` array will be updated.
    """

    non_zero_baryon = np.where(gals["StellarMass"][:] + gals["ColdGas"][:] > 0.0)[0]
    baryon_mass = np.log10(
        (gals["StellarMass"][:][non_zero_baryon] + gals["ColdGas"][:][non_zero_baryon])
        * 1.0e10 / model.hubble_h
    )

    (counts, _) = np.histogram(baryon_mass, bins=model.bins["stellar_mass_bins"])
    model.properties[f"snapshot_{snapshot}"]["BMF"] += counts


def calc_GMF(model, gals, snapshot: int):
    """
    Calculates the gas mass function of the given galaxies.  That is, the number of galaxies at a given cold gas mass.

    The ``Model.properties["snapshot_<snapshot>"]["GMF"]`` array will be updated.
    """

    non_zero_cold = np.where(gals["ColdGas"][:] > 0.0)[0]
    cold_mass = np.log10(gals["ColdGas"][:][non_zero_cold] * 1.0e10 / model.hubble_h)

    (counts, _) = np.histogram(cold_mass, bins=model.bins["stellar_mass_bins"])
    model.properties[f"snapshot_{snapshot}"]["GMF"] += counts


def calc_BTF(model, gals, snapshot: int):
    """
    Calculates the baryonic Tully-Fisher relation for spiral galaxies in the given set of galaxies.

    The number of galaxies added to ``Model.properties["snapshot_<snapshot>"]["BTF_mass"]`` and
    ``Model.properties["snapshot_<snapshot>"]["BTF_vel"]`` arrays is given by
    :py:attr:`~sage_analysis.model.Model.sample_size` weighted by ``number_spirals_passed /``
    :py:attr:`~sage_analysis.model.Model._num_gals_all_files`. If this value is greater than ``number_spirals_passed``,
    then all spiral galaxies will be used.
    """

    # Make sure we're getting spiral galaxies. That is, don't include galaxies that are too bulgy.
    w = np.where(gals["StellarMass"][:] > 0.0)[0]  # This will ensure we don't get divide by 0 errors.
    spirals = np.where((gals["Type"][:][w] == 0) & (gals["StellarMass"][:][w] + gals["ColdGas"][:][w] > 0.0) &
                       (gals["StellarMass"][:][w] > 0.0) & (gals["ColdGas"][:][w] > 0.0) &
                       (gals["BulgeMass"][:][w] / gals["StellarMass"][:][w] > 0.1) &
                       (gals["BulgeMass"][:][w] / gals["StellarMass"][:][w] < 0.5))[0]

    if len(spirals) == 0:
        logger.info(f"Could not find any spiral galaxies for analysis of the baryonic Tully-Fisher relationship.")
        return

    # Careful here, ``spirals`` is selecting on ``w``.  We want to select on ``gals``.
    spirals = w[spirals]

    # Select a random subset of galaxies (if necessary).
    spirals = model.select_random_galaxy_indices(spirals, len(model.properties[f"snapshot_{snapshot}"]["BTF_mass"]))

    baryon_mass = np.log10((gals["StellarMass"][:][spirals] + gals["ColdGas"][:][spirals]) * 1.0e10 / model.hubble_h)
    velocity = np.log10(gals["Vmax"][:][spirals])

    model.properties[f"snapshot_{snapshot}"]["BTF_mass"] = np.append(
            model.properties[f"snapshot_{snapshot}"]["BTF_mass"], baryon_mass
    )
    model.properties[f"snapshot_{snapshot}"]["BTF_vel"] = np.append(
        model.properties[f"snapshot_{snapshot}"]["BTF_vel"], velocity
    )


def calc_sSFR(model, gals, snapshot: int):
    """
    Calculates the specific star formation rate (star formation divided by the stellar mass of the galaxy) as a
    function of stellar mass.

    The number of galaxies added to ``Model.properties["snapshot_<snapshot>"]["sSFR_mass"]`` and
    ``Model.properties["snapshot_<snapshot>"]["sSFR_sSFR"]`` arrays is given by
    :py:attr:`~sage_analysis.model.Model.sample_size` weighted by ``number_gals_passed /``
    :py:attr:`~sage_analysis.model.Model._num_gals_all_files`. If this value is greater than ``number_gals_passed``,
    then all galaxies with non-zero stellar mass will be used.
    """

    non_zero_stellar = np.where(
        (gals["StellarMass"][:] > 0.0) & (gals["SfrDisk"][:] + gals["SfrBulge"][:] > 0.0)
    )[0]

    # Select a random subset of galaxies (if necessary).
    random_inds = model.select_random_galaxy_indices(non_zero_stellar, len(model.properties[f"snapshot_{snapshot}"]["sSFR_mass"]))

    stellar_mass = np.log10(gals["StellarMass"][:][random_inds] * 1.0e10 / model.hubble_h)
    sSFR = (gals["SfrDisk"][:][random_inds] + gals["SfrBulge"][:][random_inds]) / \
        (gals["StellarMass"][:][random_inds] * 1.0e10 / model.hubble_h)

    model.properties[f"snapshot_{snapshot}"]["sSFR_mass"] = np.append(
        model.properties[f"snapshot_{snapshot}"]["sSFR_mass"], stellar_mass
    )
    model.properties[f"snapshot_{snapshot}"]["sSFR_sSFR"] = np.append(
        model.properties[f"snapshot_{snapshot}"]["sSFR_sSFR"], np.log10(sSFR)
    )


def calc_gas_fraction(model, gals, snapshot: int):
    """
    Calculates the fraction of baryons that are in the cold gas reservoir as a function of stellar mass.

    The number of galaxies added to ``Model.properties["snapshot_<snapshot>"]["gas_frac_mass"]`` and
    ``Model.properties["snapshot_<snapshot>"]["gas_frac"]`` arrays is given by
    :py:attr:`~sage_analysis.model.Model.sample_size` weighted by ``number_spirals_passed /``
    :py:attr:`~sage_analysis.model.Model._num_gals_all_files`. If this value is greater than ``number_spirals_passed``,
    then all spiral galaxies will be used.
    """

    # Make sure we're getting spiral galaxies. That is, don't include galaxies that are too bulgy.
    w = np.where(gals["StellarMass"][:] > 0.0)[0]  # This will ensure we don't get divide by 0 errors.
    spirals = np.where((gals["Type"][:][w] == 0) & (gals["StellarMass"][:][w] + gals["ColdGas"][:][w] > 0.0) &
                       (gals["BulgeMass"][:][w] / gals["StellarMass"][:][w] > 0.1) &
                       (gals["BulgeMass"][:][w] / gals["StellarMass"][:][w] < 0.5))[0]

    # Careful here, ``spirals`` is selecting on ``w``.  We want to select on ``gals``.
    spirals = w[spirals]

    # Select a random subset of galaxies (if necessary).
    spirals = model.select_random_galaxy_indices(
        spirals, len(model.properties[f"snapshot_{snapshot}"]["gas_frac_mass"])
    )

    stellar_mass = np.log10(gals["StellarMass"][:][spirals] * 1.0e10 / model.hubble_h)
    gas_fraction = gals["ColdGas"][:][spirals] / (gals["StellarMass"][:][spirals] + gals["ColdGas"][:][spirals])

    model.properties[f"snapshot_{snapshot}"]["gas_frac_mass"] = np.append(
        model.properties[f"snapshot_{snapshot}"]["gas_frac_mass"], stellar_mass
    )
    model.properties[f"snapshot_{snapshot}"]["gas_frac"] = np.append(
        model.properties[f"snapshot_{snapshot}"]["gas_frac"], gas_fraction
    )


def calc_metallicity(model, gals, snapshot: int):
    """
    Calculates the metallicity as a function of stellar mass.

    The number of galaxies added to ``Model.properties["snapshot_<snapshot>"]["metallicity_mass"]`` and
    ``Model.properties["snapshot_<snapshot>"]["metallicity"]`` arrays is given by
    :py:attr:`~sage_analysis.model.Model.sample_size` weighted by ``number_centrals_passed /``
    :py:attr:`~sage_analysis.model.Model._num_gals_all_files`. If this value is greater than
    ``number_centrals_passed``, then all central galaxies will be used.
    """

    # Only care about central galaxies (Type 0) that have appreciable mass.
    w = np.where(gals["StellarMass"][:] > 0.0)[0]  # This will ensure we don't get divide by 0 errors.
    centrals = np.where(
        (gals["Type"][:][w] == 0) &
        (gals["ColdGas"][:][w] > 0.0) &
        (gals["MetalsColdGas"][:][w] > 0.0) &
        (gals["ColdGas"][:][w] / (gals["StellarMass"][:][w] + gals["ColdGas"][:][w]) > 0.1)
    )[0]

    # Careful here, ``centrals`` is selecting on ``w``.  We want to select on ``gals``.
    centrals = w[centrals]

    # Select a random subset of galaxies (if necessary).
    centrals = model.select_random_galaxy_indices(
        centrals, len(model.properties[f"snapshot_{snapshot}"]["metallicity_mass"])
    )

    stellar_mass = np.log10(gals["StellarMass"][:][centrals] * 1.0e10 / model.hubble_h)
    Z = np.log10((gals["MetalsColdGas"][:][centrals] / gals["ColdGas"][:][centrals]) / 0.02) + 9.0

    model.properties[f"snapshot_{snapshot}"]["metallicity_mass"] = np.append(
        model.properties[f"snapshot_{snapshot}"]["metallicity_mass"], stellar_mass
    )
    model.properties[f"snapshot_{snapshot}"]["metallicity"] = np.append(
        model.properties[f"snapshot_{snapshot}"]["metallicity"], Z
    )


def calc_bh_bulge(model, gals, snapshot: int):
    """
    Calculates the black hole mass as a function of bulge mass.

    The number of galaxies added to ``Model.properties["snapshot_<snapshot>"]["BlackHoleMass"]`` and
    ``Model.propertiesp["snapshot_<snapshot>"]["BulgeMass"]`` arrays is given by
    :py:attr:`~sage_analysis.model.Model.sample_size` weighted by ``number_galaxies_passed /``
    :py:attr:`~sage_analysis.model.Model._num_gals_all_files`. If this value is greater than
    ``number_galaxies_passed``, then all galaxies will be used.

    Notes
    -----
    We only consider galaxies with bulge mass greater than 10^8 Msun/h and a black hole mass greater than 10^5 Msun/h.
    """

    # Only care about galaxies that have appreciable masses.
    my_gals = np.where((gals["BulgeMass"][:] > 0.01) & (gals["BlackHoleMass"][:] > 0.00001))[0]

    # Select a random subset of galaxies (if necessary).
    my_gals = model.select_random_galaxy_indices(
        my_gals, len(model.properties[f"snapshot_{snapshot}"]["bh_mass"])
    )

    bh = np.log10(gals["BlackHoleMass"][:][my_gals] * 1.0e10 / model.hubble_h)
    bulge = np.log10(gals["BulgeMass"][:][my_gals] * 1.0e10 / model.hubble_h)

    model.properties[f"snapshot_{snapshot}"]["bh_mass"] = np.append(
        model.properties[f"snapshot_{snapshot}"]["bh_mass"], bh
    )
    model.properties[f"snapshot_{snapshot}"]["bulge_mass"] = np.append(
        model.properties[f"snapshot_{snapshot}"]["bulge_mass"], bulge
    )


def calc_quiescent(model, gals, snapshot: int):
    """
    Calculates the quiescent galaxy fraction as a function of stellar mass.  The galaxy population is also split into
    central and satellites and the quiescent fraction of these are calculated.

    The ``Model.properties["snapshot_<snapshot>"]["centrals_MF"]``,
    ``Model.properties["snapshot_<snapshot>"]["satellites_MF"]``,
    ``Model.properties["snapshot_<snapshot>"]["quiescent_galaxy_counts"]``,
    ``Model.properties["snapshot_<snapshot>"]["quiescent_centrals_counts"]``, and
    ``Model.properties["snapshot_<snapshot>"]["quiescent_satellites_counts"]`` arrays will be updated.

    Notes
    -----
    We only **count** the number of quiescent galaxies in each stellar mass bin.  When converting this to the quiescent
    fraction, one must divide by the number of galaxies in each stellar mass bin, the stellar mass function
    ``Model.properties["snapshot_<snapshot>"]["SMF"]``. See :func:`~sage_analysis.example_plots.plot_quiescent` for an
    example implementation.
    """

    non_zero_stellar = np.where(gals["StellarMass"][:] > 0.0)[0]

    mass = np.log10(gals["StellarMass"][:][non_zero_stellar] * 1.0e10 / model.hubble_h)
    gal_type = gals["Type"][:][non_zero_stellar]

    # For the sSFR, we will create a mask that is True for quiescent galaxies and
    # False for star-forming galaxies.
    sSFR = (gals["SfrDisk"][:][non_zero_stellar] + gals["SfrBulge"][:][non_zero_stellar]) / \
           (gals["StellarMass"][:][non_zero_stellar] * 1.0e10 / model.hubble_h)
    quiescent = sSFR < 10.0 ** model._sSFRcut

    # Mass function for number of centrals/satellites.
    centrals_counts, _ = np.histogram(mass[gal_type == 0], bins=model.bins["stellar_mass_bins"])
    model.properties[f"snapshot_{snapshot}"]["centrals_MF"] += centrals_counts

    satellites_counts, _ = np.histogram(mass[gal_type == 1], bins=model.bins["stellar_mass_bins"])
    model.properties[f"snapshot_{snapshot}"]["satellites_MF"] += satellites_counts

    # Then bin those galaxies/centrals/satellites that are quiescent.
    quiescent_counts, _ = np.histogram(mass[quiescent], bins=model.bins["stellar_mass_bins"])
    model.properties[f"snapshot_{snapshot}"]["quiescent_galaxy_counts"] += quiescent_counts

    quiescent_centrals_counts, _ = np.histogram(mass[(gal_type == 0) & quiescent],
                                                bins=model.bins["stellar_mass_bins"])
    model.properties[f"snapshot_{snapshot}"]["quiescent_centrals_counts"] += quiescent_centrals_counts

    quiescent_satellites_counts, _ = np.histogram(mass[(gal_type == 1) & quiescent],
                                                  bins=model.bins["stellar_mass_bins"])

    model.properties[f"snapshot_{snapshot}"]["quiescent_satellites_counts"] += quiescent_satellites_counts


def calc_bulge_fraction(model, gals, snapshot: int):
    """
    Calculates the ``bulge_mass / stellar_mass`` and ``disk_mass / stellar_mass`` ratios as a function of stellar mass.

    The ``Model.properties["snapshot_<snapshot>"]["fraction_bulge_sum"]``,
    ``Model.properties["snapshot_<snapshot>"]["fraction_disk_sum"]``,
    ``Model.properties["snapshot_<snapshot>"]["fraction_bulge_var"]``,
    ``Model.properties["snapshot_<snapshot>"]["fraction_disk_var"]`` arrays will be updated.

    Notes
    -----

    We only **sum** the bulge/disk mass in each stellar mass bin.  When converting this to the mass fraction, one must
    divide by the number of galaxies in each stellar mass bin, the stellar mass function
    ``Model.properties["snapshot_<snapshot>"]["SMF"]``. See :func:`~sage_analysis.example_plots.plot_bulge_fraction`
    for full implementation.
    """

    non_zero_stellar = np.where(gals["StellarMass"][:] > 0.0)[0]

    stellar_mass = np.log10(gals["StellarMass"][:][non_zero_stellar] * 1.0e10 / model.hubble_h)
    fraction_bulge = gals["BulgeMass"][:][non_zero_stellar] / gals["StellarMass"][:][non_zero_stellar]
    fraction_disk = 1.0 - (gals["BulgeMass"][:][non_zero_stellar] / gals["StellarMass"][:][non_zero_stellar])

    # We want the mean bulge/disk fraction as a function of stellar mass. To allow
    # us to sum across each file, we will record the sum in each bin and then average later.
    fraction_bulge_sum, _, _ = stats.binned_statistic(stellar_mass, fraction_bulge,
                                                      statistic=np.sum,
                                                      bins=model.bins["stellar_mass_bins"])
    model.properties[f"snapshot_{snapshot}"]["fraction_bulge_sum"] += fraction_bulge_sum

    fraction_disk_sum, _, _ = stats.binned_statistic(stellar_mass, fraction_disk,
                                                     statistic=np.sum,
                                                     bins=model.bins["stellar_mass_bins"])
    model.properties[f"snapshot_{snapshot}"]["fraction_disk_sum"] += fraction_disk_sum

    # For the variance, weight these by the total number of samples we will be
    # averaging over (i.e., number of files).
    fraction_bulge_var, _, _ = stats.binned_statistic(stellar_mass, fraction_bulge,
                                                      statistic=np.var,
                                                      bins=model.bins["stellar_mass_bins"])
    model.properties[f"snapshot_{snapshot}"]["fraction_bulge_var"] += fraction_bulge_var / \
        (model._last_file_to_analyze - model._first_file_to_analyze + 1)

    fraction_disk_var, _, _ = stats.binned_statistic(
        stellar_mass, fraction_disk, statistic=np.var, bins=model.bins["stellar_mass_bins"]
    )

    model.properties[f"snapshot_{snapshot}"]["fraction_disk_var"] += fraction_disk_var / \
        (model._last_file_to_analyze - model._first_file_to_analyze + 1)


def calc_baryon_fraction(model, gals, snapshot: int):
    """
    Calculates the ``mass_baryons / halo_virial_mass`` as a function of halo virial mass for each baryon reseroivr
    (stellar, cold, hot, ejected, intra-cluster stars and black hole). Also calculates the ratio for the total baryonic
    mass.

    The ``Model.properties["snapshot_<snapshot>"]["halo_<reservoir_name>_fraction_sum"]`` arrays are updated for
    each reservoir. In addition, ``Model.properties["snapshot_<snapshot>"]["halo_baryon_fraction_sum"]`` is updated.

    Notes
    -----
    The halo virial mass we use is the **background FoF halo**, not the immediate host halo of each galaxy.

    We only **sum** the baryon mass in each stellar mass bin.  When converting this to the mass fraction, one must
    divide by the number of halos in each halo mass bin, ``Model.properties["snapshot_<snapshot>"]["fof_HMF"]``. See
    :func:`~sage_analysis.example_plots.plot_baryon_fraction` for full implementation.

    If the ``Model.properties["snapshot_<snapshot>"]["fof_HMF"]`` property, with associated bins
    ``Model.bins["halo_mass"bin"]`` have not been initialized, a ``ValueError`` is thrown.
    """

    # Careful here, our "Halo Mass Function" is only counting the *BACKGROUND FOF HALOS*.
    centrals = np.where((gals["Type"][:] == 0) & (gals["Mvir"][:] > 0.0))[0]
    centrals_fof_mass = np.log10(gals["Mvir"][:][centrals] * 1.0e10 / model.hubble_h)
    try:
        halos_binned, _ = np.histogram(centrals_fof_mass, bins=model.bins["halo_mass_bins"])
    except KeyError:
        print("The `halo_mass_bins` bin array was not initialised.")
        raise ValueError

    try:
        model.properties[f"snapshot_{snapshot}"]["fof_HMF"] += halos_binned
    except KeyError:
        print("The `fof_HMF` property was not iniitalized.")
        raise ValueError

    non_zero_mvir = np.where((gals["CentralMvir"][:] > 0.0))[0]  # Will only be dividing by this value.

    # These are the *BACKGROUND FOF HALO* for each galaxy.
    fof_halo_mass = gals["CentralMvir"][:][non_zero_mvir]
    fof_halo_mass_log = np.log10(gals["CentralMvir"][:][non_zero_mvir] * 1.0e10 / model.hubble_h)

    # We want to calculate the fractions as a function of the FoF mass. To allow
    # us to sum across each file, we will record the sum in each bin and then
    # average later.
    components = ["StellarMass", "ColdGas", "HotGas", "EjectedMass", "IntraClusterStars", "BlackHoleMass"]
    attrs_different_name = ["stars", "cold", "hot", "ejected", "ICS", "bh"]

    for (component_key, attr_name) in zip(components, attrs_different_name):

        # The bins are defined in log. However the other properties are all in 1.0e10 Msun/h.
        fraction_sum, _, _ = stats.binned_statistic(
            fof_halo_mass_log,
            gals[component_key][:][non_zero_mvir] / fof_halo_mass,
            statistic=np.sum,
            bins=model.bins["halo_mass_bins"]
        )

        dict_key = "halo_{0}_fraction_sum".format(attr_name)
        model.properties[f"snapshot_{snapshot}"][dict_key] += fraction_sum

    # Finally want the sum across all components.
    baryons = sum(gals[component_key][:][non_zero_mvir] for component_key in components)
    baryon_fraction_sum, _, _ = stats.binned_statistic(
        fof_halo_mass_log, baryons / fof_halo_mass, statistic=np.sum, bins=model.bins["halo_mass_bins"]
    )
    model.properties[f"snapshot_{snapshot}"]["halo_baryon_fraction_sum"] += baryon_fraction_sum


def calc_reservoirs(model, gals, snapshot: int):
    """
    Calculates the mass in each reservoir as a function of halo virial mass.

    The number of galaxies added to ``Model.properties["snapshot_<snapshot>"]["reservoir_mvir"]`` and
    ``Model.properties["snapshot_<snapshot>"]["reservoir_<reservoir_name>"]`` arrays is given by
    :py:attr:`~sage_analysis.model.Model.sample_size` weighted by ``number_centrals_passed /``
    :py:attr:`~sage_analysis.model.Model._num_gals_all_files`. If this value is greater than
    ``number_centrals_passed``, then all central galaxies will be used.
    """

    # To reduce scatter, only use galaxies in halos with mass > 1.0e10 Msun/h.
    centrals = np.where(
        (gals["Type"][:] == 0) & (gals["Mvir"][:] > 1.0) & (gals["StellarMass"][:] > 0.0)
    )[0]

    # Select a random subset of galaxies (if necessary).
    centrals = model.select_random_galaxy_indices(
        centrals, len(model.properties[f"snapshot_{snapshot}"]["reservoir_mvir"])
    )

    reservoirs = ["Mvir", "StellarMass", "ColdGas", "HotGas", "EjectedMass", "IntraClusterStars"]
    attribute_names = ["mvir", "stars", "cold", "hot", "ejected", "ICS"]

    for (reservoir, attribute_name) in zip(reservoirs, attribute_names):

        # Some galaxies will have a zero reservoir mass which cannot be (theroetically) logged. However, to keep the
        # length of all arrays equal, we will take the log and use the entry of ``-np.inf``.  WHen plotting, these will
        # be naturally cutoff when adjusting the axis.
        # We know this will throw an error for these galaxies, so lets temporarily disable warnings.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mass = np.log10(gals[reservoir][:][centrals] * 1.0e10 / model.hubble_h)

        # Extend the previous list of masses with these new values.
        dict_key = "reservoir_{0}".format(attribute_name)
        model.properties[f"snapshot_{snapshot}"][dict_key] = np.append(
            model.properties[f"snapshot_{snapshot}"][dict_key], mass
        )


def calc_spatial(model, gals, snapshot: int):
    """
    Calculates the spatial position of the galaxies.

    The number of galaxies added to ``Model.properties["snapshot_<snapshot>"]["<x/y/z>_pos"]`` arrays is given by
    :py:attr:`~sage_analysis.model.Model.sample_size` weighted by ``number_galaxies_passed /``
    :py:attr:`~sage_analysis.model.Model._num_gals_all_files`. If this value is greater than
    ``number_galaxies_passed``, then all galaxies will be used.
    """

    non_zero_stellar = np.where((gals["Mvir"][:] > 0.0) & (gals["StellarMass"][:] > 0.1))[0]

    # Select a random subset of galaxies (if necessary).
    non_zero_stellar = model.select_random_galaxy_indices(
        non_zero_stellar, len(model.properties[f"snapshot_{snapshot}"]["x_pos"])
    )

    attribute_names = ["x_pos", "y_pos", "z_pos"]
    data_names = ["Posx", "Posy", "Posz"]

    for (attribute_name, data_name) in zip(attribute_names, data_names):

        # Units are Mpc/h.
        pos = gals[data_name][:][non_zero_stellar]
        model.properties[f"snapshot_{snapshot}"][attribute_name] = np.append(
            model.properties[f"snapshot_{snapshot}"][attribute_name], pos
        )


def calc_SMF_history(model, gals, snapshot: int):
    """
    Calculates the stellar mass function of the given galaxies.  That is, the number of galaxies at a given stellar
    mass.

    The ``Model.properties["SMF"_history]`` array will be updated.
    """
    calc_SMF(model, gals, snapshot, calc_sub_populations=False, smf_property_name="SMF_history")


def calc_SFRD_history(model, gals, snapshot: int):
    """
    Calculates the sum of the star formation across all galaxies. This will be normalized by the simulation volume to
    determine the density. See :func:`~sage_analysis.example_plots.plot_SFRD` for full implementation.

    The ``Model.properties["snapshot_<snapshot>"]["SFRD"]`` value is updated.
    """

    SFR = gals["SfrDisk"][:] + gals["SfrBulge"][:]
    model.properties[f"snapshot_{snapshot}"]["SFRD_history"] += np.sum(SFR)


def calc_SMD_history(model, gals, snapshot: int):
    """
    Calculates the sum of the stellar mass across all galaxies. This will be normalized by the simulation volume to
    determine the density. See :func:`~sage_analysis.example_plots.plot_SMD` for full implementation.

    The ``Model.properties["snapshot_<snapshot>"]["SMD"]`` value is updated.
    """

    non_zero_stellar = np.where(gals["StellarMass"][:] > 0.0)[0]
    stellar_mass = gals["StellarMass"][:][non_zero_stellar] * 1.0e10 / model.hubble_h

    model.properties[f"snapshot_{snapshot}"]["SMD_history"] += np.sum(stellar_mass)
