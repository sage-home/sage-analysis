#!/usr/bin/env python
"""
Here we show a myriad of functions that can be used to plot properties calculated from the
**SAGE** output.

We refer to :doc:`../user/plot` for more information on how plotting is handled.

Authors: (Jacob Seiler, Manodeep Sinha)
"""

from typing import List

import matplotlib
# Make the plotting scripts function without a
# valid DISPLAY variable -- MS 17/03/2020
matplotlib.use('Agg')

import numpy as np  # noqa: E402
from matplotlib import pyplot as plt  # noqa: E402

import sage_analysis.observations as obs  # noqa: E402
from sage_analysis.model import Model  # noqa: E402
from sage_analysis.plot_helper import PlotHelper  # noqa: E402


def plot_SMF(
    models: List[Model],
    snapshots: List[List[int]],
    plot_helper: PlotHelper,
    plot_sub_populations: bool = False
) -> matplotlib.figure.Figure:
    """
    Plots the stellar mass function for the specified models.

    Parameters
    ----------
    models : List of :py:class:`~sage_analysis.model.Model` class instance
        Models that will be plotted. These instances contain the properties necessary to create the plot, accessed via
        ``Model.properties["snapshot_<snapshot>"]["property_name"]``.

    snapshots : nested list of ints
        The snapshots to be plotted for each :py:class:`~sage_analysis.model.Model` in ``models``.

        The length of the outer list **MUST** be equal to the length of ``models``. For each model, the stellar mass
        function of all snapshots are plotted on the figure.

    plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`
        A helper class that contains attributes and methods to assist with plotting. In particular, the path where
        the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
        how to initialize this class and its use.

    plot_sub_populations : Boolean, default False
        If ``True``, plots the stellar mass function for red and blue sub-populations.

    Generates
    ---------
    The plot will be saved as "<output_path>1.StellarMassFunction.<output_format>"
    """

    fig = plt.figure(figsize=plot_helper.figsize)
    ax = fig.add_subplot(111)

    # Go through each of the models and plot.
    for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

        # Set the x-axis values to be the centre of the bins.
        bin_widths = model.bins["stellar_mass_bins"][1::] - model.bins["stellar_mass_bins"][0:-1]
        bin_middles = model.bins["stellar_mass_bins"][:-1] + bin_widths

        # The SMF is normalized by the simulation volume which is in Mpc/h.
        normalization_factor = model._volume / pow(model.hubble_h, 3) * bin_widths

        # Colour will be used for the model, linestyle for the snapshot.
        color = plot_helper.colors[model_num]
        label = model.label

        for snapshot_num, snapshot in enumerate(model_snapshots):
            ls = plot_helper.linestyles[snapshot_num]

            norm_SMF = model.properties[f"snapshot_{snapshot}"]["SMF"] / normalization_factor
            ax.plot(
                bin_middles,
                norm_SMF,
                color=color,
                ls=ls,
                label=f"{label} - z = {model._redshifts[snapshot]:.2f} - All",
            )

            if plot_sub_populations:
                norm_red = model.properties[f"snapshot_{snapshot}"]["red_SMF"] / normalization_factor
                norm_blue = model.properties[f"snapshot_{snapshot}"]["blue_SMF"] / normalization_factor

                ax.plot(
                    bin_middles,
                    norm_red,
                    color=color,
                    ls=plot_helper.linestyles[model_num+1],
                    lw=2,
                    label=f"{label} - z = {model._redshifts[snapshot]:.2f} - Red"
                )
                ax.plot(
                    bin_middles,
                    norm_blue,
                    color=color,
                    ls=plot_helper.linestyles[model_num+2],
                    lw=2,
                    label=f"{label} - z = {model._redshifts[snapshot]:.2f} - Blue"
                )

    # For scaling the observational data, we use the values of the zeroth
    # model.
    zeroth_hubble_h = models[0].hubble_h
    zeroth_IMF = models[0].IMF
    ax = obs.plot_smf_data(ax, zeroth_hubble_h, zeroth_IMF)

    ax.set_xlabel(r"$\log_{10} M_{\mathrm{stars}}\ (M_{\odot})$")
    ax.set_ylabel(r"$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$")

    ax.set_yscale("log", nonpositive="clip")

    ax.set_xlim([8.0, 12.0])
    ax.set_ylim([1.0e-6, 1.0e-1])

    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

    plot_helper.adjust_legend(ax, location="lower left", scatter_plot=0)

    fig.tight_layout()

    output_file = f"{plot_helper.output_path}1.StellarMassFunction.{plot_helper.output_format}"
    fig.savefig(output_file)
    print(f"Saved file to {output_file}")
    plt.close()

    return fig


def plot_BMF(
    models: List[Model],
    snapshots: List[List[int]],
    plot_helper: PlotHelper,
) -> matplotlib.figure.Figure:
    """
    Plots the baryonic mass function for the specified models. This is the mass
    function for the stellar mass + cold gas.

    Parameters
    ----------
    models : List of :py:class:`~sage_analysis.model.Model` class instance
        Models that will be plotted. These instances contain the properties necessary to create the plot, accessed via
        ``Model.properties["snapshot_<snapshot>"]["property_name"]``.

    snapshots : nested list of ints
        The snapshots to be plotted for each :py:class:`~sage_analysis.model.Model` in ``models``.

        The length of the outer list **MUST** be equal to the length of ``models``. For each model, the baryonic mass
        function of all snapshots are plotted on the figure.

    plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`
        A helper class that contains attributes and methods to assist with plotting. In particular, the path where
        the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
        how to initialize this class and its use.

    Generates
    ---------
    The plot will be saved as "<output_path>2.BaryonicMassFunction.<output_format>"
    """

    fig = plt.figure(figsize=plot_helper.figsize)
    ax = fig.add_subplot(111)

    for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

        # Set the x-axis values to be the centre of the bins.
        bin_widths = model.bins["stellar_mass_bins"][1::] - model.bins["stellar_mass_bins"][0:-1]
        bin_middles = model.bins["stellar_mass_bins"][:-1] + bin_widths

        # The MF is normalized by the simulation volume which is in Mpc/h.
        normalization_factor = model._volume / pow(model.hubble_h, 3) * bin_widths

        # Colour will be used for the snapshot, linestyle for the model.
        ls = plot_helper.linestyles[model_num]
        label = model.label

        for snapshot_num, snapshot in enumerate(model_snapshots):
            color = plot_helper.colors[snapshot_num]
            ax.plot(
                bin_middles,
                model.properties[f"snapshot_{snapshot}"]["BMF"] / normalization_factor,
                color=color,
                ls=ls,
                label=f"{label} - z = {model._redshifts[snapshot]:.2f} - All",
            )

    # For scaling the observational data, we use the values of the zeroth
    # model.
    zeroth_hubble_h = models[0].hubble_h
    zeroth_IMF = models[0].IMF
    ax = obs.plot_bmf_data(ax, zeroth_hubble_h, zeroth_IMF)

    ax.set_xlabel(r"$\log_{10}\ M_{\mathrm{bar}}\ (M_{\odot})$")
    ax.set_ylabel(r"$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$")

    ax.set_yscale("log", nonpositive="clip")

    ax.set_xlim([7.8, 12.2])
    ax.set_ylim([1.0e-6, 1.0e-1])

    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

    plot_helper.adjust_legend(ax, location="lower left", scatter_plot=0)

    fig.tight_layout()

    output_file = f"{plot_helper.output_path}2.BaryonicMassFunction.{plot_helper.output_format}"
    fig.savefig(output_file)
    print(f"Saved file to {output_file}")
    plt.close()

    return fig


def plot_GMF(
    models: List[Model],
    snapshots: List[List[int]],
    plot_helper: PlotHelper,
) -> matplotlib.figure.Figure:
    """
    Plots the gas mass function for the specified models. This is the mass function for the cold gas.

    Parameters
    ----------
    models : List of :py:class:`~sage_analysis.model.Model` class instance
        Models that will be plotted. These instances contain the properties necessary to create the plot, accessed via
        ``Model.properties["snapshot_<snapshot>"]["property_name"]``.

    snapshots : nested list of ints
        The snapshots to be plotted for each :py:class:`~sage_analysis.model.Model` in ``models``.

        The length of the outer list **MUST** be equal to the length of ``models``. For each model, the gas mass
        function of all snapshots are plotted on the figure.

    plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`
        A helper class that contains attributes and methods to assist with plotting. In particular, the path where
        the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
        how to initialize this class and its use.

    Generates
    ---------
    The plot will be saved as "<output_path>3.GasMassFunction.<output_format>"
    """

    fig = plt.figure(figsize=plot_helper.figsize)
    ax = fig.add_subplot(111)

    for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

        # Set the x-axis values to be the centre of the bins.
        bin_widths = model.bins["stellar_mass_bins"][1::] - model.bins["stellar_mass_bins"][0:-1]
        bin_middles = model.bins["stellar_mass_bins"][:-1] + bin_widths

        # The GMF is normalized by the simulation volume which is in Mpc/h.
        normalization_factor = model._volume / pow(model.hubble_h, 3) * bin_widths

        # Colour will be used for the snapshot, linestyle for the model.
        ls = plot_helper.linestyles[model_num]
        label = model.label

        for snapshot_num, snapshot in enumerate(model_snapshots):
            color = plot_helper.colors[snapshot_num]

            ax.plot(
                bin_middles,
                model.properties[f"snapshot_{snapshot}"]["GMF"] / normalization_factor,
                color=color,
                ls=ls,
                label=f"{label} - z = {model._redshifts[snapshot]:.2f} - Cold Gas",
            )

    # For scaling the observational data, we use the values of the zeroth
    # model.
    zeroth_hubble_h = models[0].hubble_h
    obs.plot_gmf_data(ax, zeroth_hubble_h)

    ax.set_xlabel(r"$\log_{10} M_{\mathrm{X}}\ (M_{\odot})$")
    ax.set_ylabel(r"$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$")

    ax.set_yscale("log", nonpositive="clip")

    # Find the models that have the smallest/largest stellar mass bin.
    ax.set_xlim([7.8, 12.2])
    ax.set_ylim([1.0e-6, 1.0e-1])

    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

    plot_helper.adjust_legend(ax, location="lower left", scatter_plot=0)

    fig.tight_layout()

    output_file = f"{plot_helper.output_path}3.GasMassFunction.{plot_helper.output_format}"
    fig.savefig(output_file)  # Save the figure
    print(f"Saved file to {output_file}")
    plt.close()

    return fig


def plot_BTF(
    models: List[Model],
    snapshots: List[List[int]],
    plot_helper: PlotHelper,
) -> matplotlib.figure.Figure:
    """
    Plots the baryonic Tully-Fisher relationship for the specified models.

    Parameters
    ----------
    models : List of :py:class:`~sage_analysis.model.Model` class instance
        Models that will be plotted. These instances contain the properties necessary to create the plot, accessed via
        ``Model.properties["snapshot_<snapshot>"]["property_name"]``.

    snapshots : nested list of ints
        The snapshots to be plotted for each :py:class:`~sage_analysis.model.Model` in ``models``.

        The length of the outer list **MUST** be equal to the length of ``models``. For each model, the baryonic
        Tully-Fisher relationship of all snapshots are plotted on the figure.

    plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`
        A helper class that contains attributes and methods to assist with plotting. In particular, the path where
        the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
        how to initialize this class and its use.

    Generates
    ---------
    The plot will be saved as "<output_path>4.BaryonicTullyFisher.<output_format>"
    """

    fig = plt.figure(figsize=plot_helper.figsize)
    ax = fig.add_subplot(111)

    # Go through each of the models and plot.
    for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

        # Colour will be used for the model, marker style for the snapshot.
        color = plot_helper.colors[model_num]
        label = model.label

        for snapshot_num, snapshot in enumerate(model_snapshots):
            marker = plot_helper.markers[snapshot_num]

            ax.scatter(
                model.properties[f"snapshot_{snapshot}"]["BTF_vel"],
                model.properties[f"snapshot_{snapshot}"]["BTF_mass"],
                marker=marker,
                s=5,
                color=color,
                alpha=0.8,
                label=f"{label} - z = {model._redshifts[snapshot]:.2f} - Sb/c galaxies",
            )

    ax.set_xlim([1.4, 2.6])
    ax.set_ylim([8.0, 12.0])

    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))

    ax.set_xlabel(r"$\log_{10}V_{max}\ (km/s)$")
    ax.set_ylabel(r"$\log_{10}\ M_{\mathrm{bar}}\ (M_{\odot})$")

    ax = obs.plot_btf_data(ax)

    plot_helper.adjust_legend(ax, location="upper left", scatter_plot=1)

    fig.tight_layout()

    output_file = f"{plot_helper.output_path}4.BaryonicTullyFisher.{plot_helper.output_format}"
    fig.savefig(output_file)
    print(f"Saved file to {output_file}")
    plt.close()

    return fig


def plot_sSFR(
    models: List[Model],
    snapshots: List[List[int]],
    plot_helper: PlotHelper,
) -> matplotlib.figure.Figure:
    """
    Plots the specific star formation rate as a function of stellar mass for the specified
    models.

    Parameters
    ----------
    models : List of :py:class:`~sage_analysis.model.Model` class instance
        Models that will be plotted. These instances contain the properties necessary to create the plot, accessed via
        ``Model.properties["snapshot_<snapshot>"]["property_name"]``.

    snapshots : nested list of ints
        The snapshots to be plotted for each :py:class:`~sage_analysis.model.Model` in ``models``.

        The length of the outer list **MUST** be equal to the length of ``models``. For each model, the specific star
        formation rate of all snapshots are plotted on the figure.

    plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`
        A helper class that contains attributes and methods to assist with plotting. In particular, the path where
        the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
        how to initialize this class and its use.

    Generates
    ---------
    The plot will be saved as "<output_path>5.SpecificStarFormationRate.<output_format>"
    """

    fig = plt.figure(figsize=plot_helper.figsize)
    ax = fig.add_subplot(111)

    # Go through each of the models and plot.
    for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

        # Colour will be used for the model, marker style for the snapshot.
        color = plot_helper.colors[model_num]
        label = model.label

        for snapshot_num, snapshot in enumerate(model_snapshots):
            marker = plot_helper.markers[snapshot_num]

            ax.scatter(
                model.properties[f"snapshot_{snapshot}"]["sSFR_mass"],
                model.properties[f"snapshot_{snapshot}"]["sSFR_sSFR"],
                marker=marker,
                s=5,
                color=color,
                alpha=0.8,
                label=f"{label} - z = {model._redshifts[snapshot]:.2f}",
            )

    # Overplot a dividing line between passive and SF galaxies.
    w = np.arange(7.0, 13.0, 1.0)
    min_sSFRcut = np.min([model._sSFRcut for model in models])
    ax.plot(w, w/w*min_sSFRcut, "b:", lw=2.0)

    ax.set_xlabel(r"$\log_{10} M_{\mathrm{stars}}\ (M_{\odot})$")
    ax.set_ylabel(r"$\log_{10}\ s\mathrm{SFR}\ (\mathrm{yr^{-1}})$")

    ax.set_xlim([8.0, 12.0])
    ax.set_ylim([-16.0, -8.0])

    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))

    plot_helper.adjust_legend(ax, scatter_plot=1)

    fig.tight_layout()

    output_file = f"{plot_helper.output_path}5.SpecificStarFormationRate.{plot_helper.output_format}"
    fig.savefig(output_file)
    print(f"Saved file to {output_file}")
    plt.close()

    return fig


def plot_gas_fraction(
    models: List[Model],
    snapshots: List[List[int]],
    plot_helper: PlotHelper,
) -> matplotlib.figure.Figure:
    """
    Plots the fraction of baryons that are in the cold gas reservoir as a function of
    stellar mass for the specified models.

    Parameters
    ----------
    models : List of :py:class:`~sage_analysis.model.Model` class instance
        Models that will be plotted. These instances contain the properties necessary to create the plot, accessed via
        ``Model.properties["snapshot_<snapshot>"]["property_name"]``.

    snapshots : nested list of ints
        The snapshots to be plotted for each :py:class:`~sage_analysis.model.Model` in ``models``.

        The length of the outer list **MUST** be equal to the length of ``models``. For each model, the gas fraction
        of all snapshots are plotted on the figure.

    plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`
        A helper class that contains attributes and methods to assist with plotting. In particular, the path where
        the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
        how to initialize this class and its use.

    Generates
    ---------
    The plot will be saved as "<output_path>6.GasFraction.<output_format>"
    """

    fig = plt.figure(figsize=plot_helper.figsize)
    ax = fig.add_subplot(111)

    # Go through each of the models and plot.
    for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

        # Colour will be used for the model, marker style for the snapshot.
        color = plot_helper.colors[model_num]
        label = model.label

        for snapshot_num, snapshot in enumerate(model_snapshots):
            marker = plot_helper.markers[snapshot_num]

            ax.scatter(
                model.properties[f"snapshot_{snapshot}"]["gas_frac_mass"],
                model.properties[f"snapshot_{snapshot}"]["gas_frac"],
                marker=marker,
                s=20,
                color=color,
                alpha=0.5,
                label=f"{label} - z = {model._redshifts[snapshot]:.2f} - Sb/c galaxies",
            )

    ax.set_xlabel(r"$\log_{10} M_{\mathrm{stars}}\ (M_{\odot})$")
    ax.set_ylabel(r"$\mathrm{Cold\ Mass\ /\ (Cold+Stellar\ Mass)}$")

    ax.set_xlim([7.8, 12.2])
    ax.set_ylim([0.0, 1.0])

    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))

    plot_helper.adjust_legend(ax, scatter_plot=1)

    fig.tight_layout()

    output_file = f"{plot_helper.output_path}6.GasFraction.{plot_helper.output_format}"
    fig.savefig(output_file)
    print(f"Saved file to {output_file}")
    plt.close()

    return fig


def plot_metallicity(
    models: List[Model],
    snapshots: List[List[int]],
    plot_helper: PlotHelper,
) -> matplotlib.figure.Figure:
    """
    Plots the metallicity as a function of stellar mass for the speicifed models.

    Parameters
    ----------
    models : List of :py:class:`~sage_analysis.model.Model` class instance
        Models that will be plotted. These instances contain the properties necessary to create the plot, accessed via
        ``Model.properties["snapshot_<snapshot>"]["property_name"]``.

    snapshots : nested list of ints
        The snapshots to be plotted for each :py:class:`~sage_analysis.model.Model` in ``models``.

        The length of the outer list **MUST** be equal to the length of ``models``. For each model, the metallicity
        of all snapshots are plotted on the figure.

    plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`
        A helper class that contains attributes and methods to assist with plotting. In particular, the path where
        the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
        how to initialize this class and its use.

    Generates
    ---------
    The plot will be saved as "<output_path>7.Metallicity.<output_format>"
    """

    fig = plt.figure(figsize=plot_helper.figsize)
    ax = fig.add_subplot(111)

    # Go through each of the models and plot.
    for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

        # Colour will be used for the model, marker style for the snapshot.
        color = plot_helper.colors[model_num]
        label = model.label

        for snapshot_num, snapshot in enumerate(model_snapshots):
            marker = plot_helper.markers[snapshot_num]

            ax.scatter(
                model.properties[f"snapshot_{snapshot}"]["metallicity_mass"],
                model.properties[f"snapshot_{snapshot}"]["metallicity"],
                marker=marker,
                s=5,
                color=color,
                alpha=0.8,
                label=f"{label} - z = {model._redshifts[snapshot]:.2f}",
            )

    # Use the IMF of the zeroth model to scale the observational results.
    zeroth_IMF = models[0].IMF
    ax = obs.plot_metallicity_data(ax, zeroth_IMF)

    ax.set_xlabel(r"$\log_{10} M_{\mathrm{stars}}\ (M_{\odot})$")
    ax.set_ylabel(r"$12\ +\ \log_{10}[\mathrm{O/H}]$")

    ax.set_xlim([7.8, 12.2])
    ax.set_ylim([8.0, 9.5])

    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))

    # Since we're doing a scatter plot, we need to resize the legend points.
    plot_helper.adjust_legend(ax, location="upper right", scatter_plot=1)

    fig.tight_layout()

    output_file = f"{plot_helper.output_path}7.Metallicity.{plot_helper.output_format}"
    fig.savefig(output_file)
    print(f"Saved file to {output_file}")
    plt.close()

    return fig


def plot_bh_bulge(
    models: List[Model],
    snapshots: List[List[int]],
    plot_helper: PlotHelper,
) -> matplotlib.figure.Figure:
    """
    Plots the black-hole bulge relationship for the specified models.

    Parameters
    ----------
    models : List of :py:class:`~sage_analysis.model.Model` class instance
        Models that will be plotted. These instances contain the properties necessary to create the plot, accessed via
        ``Model.properties["snapshot_<snapshot>"]["property_name"]``.

    snapshots : nested list of ints
        The snapshots to be plotted for each :py:class:`~sage_analysis.model.Model` in ``models``.

        The length of the outer list **MUST** be equal to the length of ``models``. For each model, the black hole
        bulge relationship of all snapshots are plotted on the figure.

    plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`
        A helper class that contains attributes and methods to assist with plotting. In particular, the path where
        the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
        how to initialize this class and its use.

    Generates
    ---------
    The plot will be saved as "<output_path>8.BlackHoleBulgeRelationship.<output_format>"
    """

    fig = plt.figure(figsize=plot_helper.figsize)
    ax = fig.add_subplot(111)

    # Go through each of the models and plot.
    for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

        # Colour will be used for the model, marker style for the snapshot.
        color = plot_helper.colors[model_num]
        label = model.label

        for snapshot_num, snapshot in enumerate(model_snapshots):
            marker = plot_helper.markers[snapshot_num]

            ax.scatter(
                model.properties[f"snapshot_{snapshot}"]["bulge_mass"],
                model.properties[f"snapshot_{snapshot}"]["bh_mass"],
                marker=marker,
                s=5,
                color=color,
                alpha=0.8,
                label=f"{label} - z = {model._redshifts[snapshot]:.2f}",
            )

    ax = obs.plot_bh_bulge_data(ax)

    ax.set_xlabel(r"$\log\ M_{\mathrm{bulge}}\ (M_{\odot})$")
    ax.set_ylabel(r"$\log\ M_{\mathrm{BH}}\ (M_{\odot})$")

    ax.set_xlim([7.8, 12.2])
    ax.set_ylim([6.0, 10.0])

    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))

    plot_helper.adjust_legend(ax, location="upper right", scatter_plot=1)

    fig.tight_layout()

    output_file = f"{plot_helper.output_path}8.BlackHoleBulgeRelationship.{plot_helper.output_format}"
    fig.savefig(output_file)
    print(f"Saved file to {output_file}")
    plt.close()

    return fig


def plot_quiescent(
    models: List[Model],
    snapshots: List[List[int]],
    plot_helper: PlotHelper,
    plot_sub_populations: bool = False
) -> matplotlib.figure.Figure:
    """
    Plots the fraction of galaxies that are quiescent as a function of stellar mass for the
    specified models.  The quiescent cut is defined by :py:attr:`~sage_analysis.model.Model.sSFRcut`.

    Parameters
    ----------
    models : List of :py:class:`~sage_analysis.model.Model` class instance
        Models that will be plotted. These instances contain the properties necessary to create the plot, accessed via
        ``Model.properties["snapshot_<snapshot>"]["property_name"]``.

    snapshots : nested list of ints
        The snapshots to be plotted for each :py:class:`~sage_analysis.model.Model` in ``models``.

        The length of the outer list **MUST** be equal to the length of ``models``. For each model, the quiescent
        fraction of all snapshots are plotted on the figure.

    plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`
        A helper class that contains attributes and methods to assist with plotting. In particular, the path where
        the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
        how to initialize this class and its use.

    plot_sub_populations : Boolean, default False
        If ``True``, plots the centrals and satellite sub-populations.

    Generates
    ---------
    The plot will be saved as "<output_path>9.QuiescentFraction.<output_format>"
    """

    fig = plt.figure(figsize=plot_helper.figsize)
    ax = fig.add_subplot(111)

    # Go through each of the models and plot.
    for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

        # Set the x-axis values to be the centre of the bins.
        bin_widths = model.bins["stellar_mass_bins"][1::] - model.bins["stellar_mass_bins"][0:-1]
        bin_middles = model.bins["stellar_mass_bins"][:-1] + bin_widths

        # Colour will be used for the snapshot, linestyle for the model.
        ls = plot_helper.linestyles[model_num]
        label = model.label

        for snapshot_num, snapshot in enumerate(model_snapshots):
            color = plot_helper.colors[snapshot_num]

            non_zero_SMF = np.where(model.properties[f"snapshot_{snapshot}"]["SMF"] > 0)
            quiescent_fraction = model.properties[f"snapshot_{snapshot}"]["quiescent_galaxy_counts"][non_zero_SMF] / \
                model.properties[f"snapshot_{snapshot}"]["SMF"][non_zero_SMF]

            ax.plot(
                bin_middles[non_zero_SMF],
                quiescent_fraction,
                color=color,
                ls=ls,
                label=f"{label} - z = {model._redshifts[snapshot]:.2f} - All",
            )

            # Be careful to not overcrowd the plot.
            if plot_sub_populations:

                non_zero_MF = np.where(model.properties[f"snapshot_{snapshot}"]["centrals_MF"])[0]
                quiescent_central_fraction = \
                    model.properties[f"snapshot_{snapshot}"]["quiescent_centrals_counts"][non_zero_MF] / \
                    model.properties[f"snapshot_{snapshot}"]["centrals_MF"][non_zero_MF]

                ax.plot(
                    bin_middles[non_zero_MF],
                    quiescent_central_fraction,
                    color=color,
                    linestyle="--",
                )

                non_zero_MF = np.where(model.properties[f"snapshot_{snapshot}"]["satellites_MF"])[0]
                quiescent_sat_fraction = \
                    model.properties[f"snapshot_{snapshot}"]["quiescent_satellites_counts"][non_zero_MF] / \
                    model.properties[f"snapshot_{snapshot}"]["satellites_MF"][non_zero_MF]

                ax.plot(
                    bin_middles[non_zero_MF],
                    quiescent_sat_fraction,
                    color=color,
                    linestyle="-.",
                    label=f"{label} - z = {model._redshifts[snapshot]:.2f} - Satellites",
                )

    ax.set_xlabel(r"$\log_{10} M_{\mathrm{stellar}}\ (M_{\odot})$")
    ax.set_ylabel(r"$\mathrm{Quescient\ Fraction}$")

    ax.set_xlim([7.8, 12.2])
    ax.set_ylim([0.0, 1.05])

    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.10))

    plot_helper.adjust_legend(ax, location="upper left", scatter_plot=0)

    fig.tight_layout()

    output_file = f"{plot_helper.output_path}9.QuiescentFraction.{plot_helper.output_format}"
    fig.savefig(output_file)
    print(f"Saved file to {output_file}")
    plt.close()

    return fig


def plot_bulge_fraction(
    models: List[Model],
    snapshots: List[List[int]],
    plot_helper: PlotHelper,
    plot_disk_fraction: bool = False,
    plot_var: bool = False
) -> matplotlib.figure.Figure:
    """
    Plots the fraction of the stellar mass that is located in the bulge/disk as a function
    of stellar mass for the specified models.

    Parameters
    ----------
    models : List of ``Model`` class instance
        Models that will be plotted. These instances contain the properties necessary to
        create the plot, accessed via ``Model.properties["snapshot_<snapshot>"]["property_name"]``.

    snapshots : nested list of ints
        The snapshots to be plotted for each :py:class:`~sage_analysis.model.Model` in ``models``.

        The length of the outer list **MUST** be equal to the length of ``models``. For each model, the bulge fraction
        of all snapshots are plotted on the figure.

    plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`
        A helper class that contains attributes and methods to assist with plotting. In particular, the path where
        the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
        how to initialize this class and its use.

    plot_disk_fraction : bool, optional
        If specified, will also plot the disk fraction.

    plot_var : Boolean, default False
        If ``True``, plots the variance as shaded regions.

    Generates
    ---------
    The plot will be saved as :py:attr:`~sage_analysis.plot_helper.PlotHelper.output_path`10.BulgeMassFraction.
    :py:attr:`~sage_analysis.plot_helper.PlotHelper.output_format`.
    """

    fig = plt.figure(figsize=plot_helper.figsize)
    ax = fig.add_subplot(111)

    # Go through each of the models and plot.
    for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

        # Set the x-axis values to be the centre of the bins.
        bin_widths = model.bins["stellar_mass_bins"][1::] - model.bins["stellar_mass_bins"][0:-1]
        bin_middles = model.bins["stellar_mass_bins"][:-1] + bin_widths

        # Colour will be used for the snapshot, linestyle for the model.
        ls = plot_helper.linestyles[model_num]
        label = model.label

        for snapshot_num, snapshot in enumerate(model_snapshots):
            color = plot_helper.colors[snapshot_num]

            # Remember we need to average the properties in each bin.
            non_zero_SMF = np.where(model.properties[f"snapshot_{snapshot}"]["SMF"] > 0)
            bulge_mean = model.properties[f"snapshot_{snapshot}"]["fraction_bulge_sum"][non_zero_SMF] / \
                model.properties[f"snapshot_{snapshot}"]["SMF"][non_zero_SMF]
            disk_mean = model.properties[f"snapshot_{snapshot}"]["fraction_disk_sum"][non_zero_SMF] / \
                model.properties[f"snapshot_{snapshot}"]["SMF"][non_zero_SMF]

            # The variance has already been weighted when we calculated it.
            bulge_var = model.properties[f"snapshot_{snapshot}"]["fraction_bulge_var"][non_zero_SMF]
            disk_var = model.properties[f"snapshot_{snapshot}"]["fraction_disk_var"][non_zero_SMF]

            ax.plot(
                bin_middles[non_zero_SMF],
                bulge_mean,
                color=color,
                ls=ls,
                label=f"{label} - z = {model._redshifts[snapshot]:.2f} - Bulge",
            )

            if plot_disk_fraction:
                ax.plot(
                    bin_middles,
                    disk_mean,
                    color=color,
                    ls=plot_helper.linestyles[model_num+1],
                    label=f"{label} - z = {model._redshifts[snapshot]:.2f} - Disk",
                )

            if plot_var:
                ax.fill_between(
                    bin_middles,
                    bulge_mean+bulge_var,
                    bulge_mean-bulge_var,
                    facecolor=color,
                    alpha=0.25
                )

                if plot_disk_fraction:
                    ax.fill_between(
                        bin_middles,
                        disk_mean+disk_var,
                        disk_mean-disk_var,
                        facecolor=color,
                        alpha=0.25
                    )

    ax.set_xlabel(r"$\log_{10} M_{\mathrm{stars}}\ (M_{\odot})$")
    ax.set_ylabel(r"$\mathrm{Stellar\ Mass\ Fraction}$")

    ax.set_xlim([7.8, 12.2])
    ax.set_ylim([0.0, 1.05])

    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.10))

    plot_helper.adjust_legend(ax, location="upper left", scatter_plot=0)

    fig.tight_layout()

    output_file = f"{plot_helper.output_path}10.BulgeMassFraction.{plot_helper.output_format}"
    fig.savefig(output_file)
    print(f"Saved file to {output_file}")
    plt.close()

    return fig


def plot_baryon_fraction(
    models: List[Model],
    snapshots: List[List[int]],
    plot_helper: PlotHelper,
    plot_sub_populations: bool = False
) -> matplotlib.figure.Figure:
    """
    Plots the total baryon fraction as afunction of halo mass for the specified models.

    Parameters
    ----------
    models : List of :py:class:`~sage_analysis.model.Model` class instance
        Models that will be plotted. These instances contain the properties necessary to create the plot, accessed via
        ``Model.properties["snapshot_<snapshot>"]["property_name"]``.

    snapshots : nested list of ints
        The snapshots to be plotted for each :py:class:`~sage_analysis.model.Model` in ``models``.

        The length of the outer list **MUST** be equal to the length of ``models``. For each model, the baryon fraction
        of all snapshots are plotted on the figure.

    plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`
        A helper class that contains attributes and methods to assist with plotting. In particular, the path where
        the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
        how to initialize this class and its use.

    plot_sub_populations : Boolean, default False
        If ``True``, plots the baryon fraction for each reservoir. Otherwise, only plots
        the total baryon fraction.

    Generates
    ---------
    The plot will be saved as "<output_path>11.BaryonFraction.<output_format>"
    """

    fig = plt.figure(figsize=plot_helper.figsize)
    ax = fig.add_subplot(111)

    # Go through each of the models and plot.
    for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

        # Set the x-axis values to be the centre of the bins.
        bin_widths = model.bins["stellar_mass_bins"][1::] - model.bins["stellar_mass_bins"][0:-1]
        bin_middles = model.bins["stellar_mass_bins"][:-1] + bin_widths

        # Colour will be used for the snapshot, linestyle for the model.
        ls = plot_helper.linestyles[model_num]
        label = model.label

        for snapshot_num, snapshot in enumerate(model_snapshots):
            color = plot_helper.colors[snapshot_num]

            # Remember we need to average the properties in each bin.
            non_zero_HMF = np.where(model.properties[f"snapshot_{snapshot}"]["fof_HMF"] > 0)
            baryon_mean = model.properties[f"snapshot_{snapshot}"]["halo_baryon_fraction_sum"][non_zero_HMF] / \
                model.properties[f"snapshot_{snapshot}"]["fof_HMF"][non_zero_HMF]

            ax.plot(
                bin_middles[non_zero_HMF],
                baryon_mean,
                label=f"{label} - z = {model._redshifts[snapshot]:.2f} - Total",
                color=color,
                ls=ls,
            )

            if plot_sub_populations:
                attrs = ["stars", "cold", "hot", "ejected", "ICS"]
                res_labels = ["Stars", "Cold", "Hot", "Ejected", "ICS"]
                res_colors = ["k", "b", "r", "g", "y"]

                for (attr, res_label, res_color) in zip(attrs, res_labels, res_colors):
                    dict_key = "halo_{0}_fraction_sum".format(attr)

                    mean = model.properties[f"snapshot_{snapshot}"][dict_key][non_zero_HMF] / \
                        model.properties[f"snapshot_{snapshot}"]["fof_HMF"][non_zero_HMF]

                    ax.plot(
                        bin_middles[non_zero_HMF],
                        mean,
                        label=f"{label} - z = {model._redshifts[snapshot]:.2f} - {res_label}",
                        color=res_color,
                        ls=ls,
                    )

    ax.set_xlabel(r"$\mathrm{Central}\ \log_{10} M_{\mathrm{vir}}\ (M_{\odot})$")
    ax.set_ylabel(r"$\mathrm{Baryon\ Fraction}$")

    ax.set_xlim([9.8, 14.2])
    ax.set_ylim([0.0, 0.23])

    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.05))

    plot_helper.adjust_legend(ax, location="upper left", scatter_plot=0)

    output_file = f"{plot_helper.output_path}11.BaryonFraction.{plot_helper.output_format}"
    fig.savefig(output_file)
    print(f"Saved file to {output_file}")
    plt.close()

    return fig


def plot_reservoirs(
    models: List[Model],
    snapshots: List[List[int]],
    plot_helper: PlotHelper,
) -> List[matplotlib.figure.Figure]:
    """
    Plots the mass in each reservoir as a function of halo mass for the
    specified models.

    Parameters
    ----------
    models : List of :py:class:`~sage_analysis.model.Model` class instance
        Models that will be plotted. These instances contain the properties necessary to create the plot, accessed via
        ``Model.properties["snapshot_<snapshot>"]["property_name"]``.

    snapshots : nested list of ints
        The snapshots to be plotted for each :py:class:`~sage_analysis.model.Model` in ``models``.

        The length of the outer list **MUST** be equal to the length of ``models``. For each model, each snapshot will
        be plotted and saved as a separate figure.

    plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`
        A helper class that contains attributes and methods to assist with plotting. In particular, the path where
        the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
        how to initialize this class and its use.

    Generates
    ---------
    A plot will be saved as ``"<output_path>12.MassReservoirs<model.label>.<output_format>"`` for each mode.
    """

    # This scatter plot will be messy so we're going to make one for each model.
    figs = []
    for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

        label = model.label
        marker = plot_helper.markers[model_num]

        # Furthermore, make one for each snapshot we requested.
        for snapshot_num, snapshot in enumerate(model_snapshots):

            fig = plt.figure(figsize=plot_helper.figsize)
            ax = fig.add_subplot(111)

            attribute_names = ["stars", "cold", "hot", "ejected", "ICS"]
            res_labels = ["Stars", "ColdGas", "HotGas", "EjectedGas", "IntraclusterStars"]
            res_colors = ["k", "b", "r", "g", "y"]

            for (attribute_name, res_label, res_color) in zip(attribute_names, res_labels, res_colors):

                dict_key = f"reservoir_{attribute_name}"
                ax.scatter(
                    model.properties[f"snapshot_{snapshot}"]["reservoir_mvir"],
                    model.properties[f"snapshot_{snapshot}"][dict_key],
                    marker=marker,
                    s=0.3,
                    label=res_label,
                    color=res_color,
                )

            ax.set_xlabel(r"$\log\ M_{\mathrm{vir}}\ (M_{\odot})$")
            ax.set_ylabel(r"$\mathrm{Reservoir\ Mass\ (M_{\odot})}$")

            ax.set_xlim([9.8, 14.2])
            ax.set_ylim([7.5, 12.5])

            ax.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
            ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))

            plot_helper.adjust_legend(ax, location="upper left", scatter_plot=1)

            fig.tight_layout()

            output_file = f"{plot_helper.output_path}12.MassReservoirs_{label}_snapshot{snapshot}.{plot_helper.output_format}"
            fig.savefig(output_file)
            print(f"Saved file to {output_file}")
            plt.close()

            figs.append(fig)

    return figs


def plot_spatial(
    models: List[Model],
    snapshots: List[List[int]],
    plot_helper: PlotHelper,
) -> matplotlib.figure.Figure:
    """
    Plots the spatial distribution of the galaxies for specified models.

    Parameters
    ----------
    models : List of :py:class:`~sage_analysis.model.Model` class instance
        Models that will be plotted. These instances contain the properties necessary to create the plot, accessed via
        ``Model.properties["snapshot_<snapshot>"]["property_name"]``.

    snapshots : nested list of ints
        The snapshots to be plotted for each :py:class:`~sage_analysis.model.Model` in ``models``.

        The length of the outer list **MUST** be equal to the length of ``models``. For each model, the spatial
        position of all snapshots are plotted on the figure.

    plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`
        A helper class that contains attributes and methods to assist with plotting. In particular, the path where
        the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
        how to initialize this class and its use.

    Generates
    ---------

    A plot will be saved as ``"<output_path>13.SpatialDistribution<model.label>.<output_format>"`` for each
    model.
    """

    fig = plt.figure(figsize=plot_helper.figsize)

    # 4-panel plot.
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    # Go through each of the models and plot.
    for model_num, (model, model_snapshots) in enumerate(zip(models, snapshots)):

        # Colour will be used for the model, marker style for the snapshot.
        color = plot_helper.colors[model_num]
        label = model.label

        for snapshot_num, snapshot in enumerate(model_snapshots):
            marker = plot_helper.markers[snapshot_num]

            ax1.scatter(
                model.properties[f"snapshot_{snapshot}"]["x_pos"],
                model.properties[f"snapshot_{snapshot}"]["y_pos"],
                marker=marker,
                s=0.3,
                color=color,
                alpha=0.5
            )

            ax2.scatter(
                model.properties[f"snapshot_{snapshot}"]["x_pos"],
                model.properties[f"snapshot_{snapshot}"]["z_pos"],
                marker=marker,
                s=0.3,
                color=color,
                alpha=0.5
            )

            ax3.scatter(
                model.properties[f"snapshot_{snapshot}"]["y_pos"],
                model.properties[f"snapshot_{snapshot}"]["z_pos"],
                marker=marker,
                s=0.3,
                color=color,
                alpha=0.5
            )

            # The bottom right panel will only contain the legend.
            # For some odd reason, plotting `np.nan` causes some legend entries to not
            # appear. Plot junk and we'll adjust the axis to not show it.
            ax4.scatter(
                -999,
                -999,
                marker=marker,
                color=color,
                label=f"{label} - z = {model._redshifts[snapshot]:.2f}",
            )
            ax4.axis("off")

    ax1.set_xlabel(r"$\mathrm{x}\ [\mathrm{Mpc}/h]$")
    ax1.set_ylabel(r"$\mathrm{y}\ [\mathrm{Mpc}/h]$")

    ax2.set_xlabel(r"$\mathrm{x}\ [\mathrm{Mpc}/h]$")
    ax2.set_ylabel(r"$\mathrm{z}\ [\mathrm{Mpc}/h]$")

    ax3.set_xlabel(r"$\mathrm{y}\ [\mathrm{Mpc}/h]$")
    ax3.set_ylabel(r"$\mathrm{z}\ [\mathrm{Mpc}/h]$")

    # Find the model with the largest box.
    max_box = np.max([model.box_size for model in models])
    buffer = max_box*0.05
    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_xlim([0.0-buffer, max_box+buffer])
        ax.set_ylim([0.0-buffer, max_box+buffer])

        ax.xaxis.set_minor_locator(plt.MultipleLocator(5))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(5))

    plot_helper.adjust_legend(ax4, location="upper left", scatter_plot=1)

    # Make sure everything remains nicely layed out.
    fig.tight_layout()

    output_file = f"{plot_helper.output_path}13.SpatialDistribution.{plot_helper.output_format}"
    fig.savefig(output_file)
    print(f"Saved file to {output_file}")
    plt.close()

    return fig


def plot_spatial_3d(pos, output_file, box_size) -> matplotlib.figure.Figure:
    """
    Plots the 3D spatial distribution of galaxies.

    Parameters
    ==========

    pos : ``numpy`` 3D array with length equal to the number of galaxies
        The position (in Mpc/h) of the galaxies.

    output_file : String
        Name of the file the plot will be saved as.

    Returns
    =======

    None.  A plot will be saved as ``output_file``.
    """

    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    from random import sample

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Generate a subsample if necessary.
    num_gals = len(pos)
    sample_size = 10000
    if num_gals > sample_size:
        w = sample(list(np.arange(num_gals)), sample_size)
    else:
        w = np.arange(num_gals)

    ax.scatter(pos[w, 0], pos[w, 1], pos[w, 2], alpha=0.5)

    ax.set_xlim([0.0, box_size])
    ax.set_ylim([0.0, box_size])
    ax.set_zlim([0.0, box_size])

    ax.set_xlabel(r"$\mathbf{x \: [h^{-1}Mpc]}$")
    ax.set_ylabel(r"$\mathbf{y \: [h^{-1}Mpc]}$")
    ax.set_zlabel(r"$\mathbf{z \: [h^{-1}Mpc]}$")

    fig.tight_layout()

    fig.savefig(output_file)
    print(f"Saved file to {output_file}")
    plt.close()

    return fig


def plot_SMF_history(
    models: List[Model],
    snapshots: List[List[int]],
    plot_helper: PlotHelper,
) -> matplotlib.figure.Figure:
    """
    Plots the evolution of the stellar mass function for the specified models.
    This function loops over the value of ``model.SMF_snaps`` and plots and the SMFs at
    each snapshots.

    Parameters
    ----------
    models : List of :py:class:`~sage_analysis.model.Model` class instance
        Models that will be plotted. These instances contain the properties necessary to create the plot, accessed via
        ``Model.properties["snapshot_<snapshot>"]["property_name"]``.

    snapshots : nested list of ints
        This is a dummy variable that is present to ensure the signature is identical to the other plot functions.

    plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`
        A helper class that contains attributes and methods to assist with plotting. In particular, the path where
        the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
        how to initialize this class and its use.

    Generates
    ---------

    The plot will be saved as "<output_path>A.StellarMassFunction.<output_format>"
    """

    fig = plt.figure(figsize=plot_helper.figsize)
    ax = fig.add_subplot(111)

    # Go through each of the models and plot.
    for (model_num, model) in enumerate(models):

        ls = plot_helper.linestyles[model_num]

        # Set the x-axis values to be the centre of the bins.
        bin_widths = model.bins["stellar_mass_bins"][1::] - model.bins["stellar_mass_bins"][0:-1]
        bin_middles = model.bins["stellar_mass_bins"][:-1] + bin_widths

        # Iterate over the snapshots.
        for snap in model._history_SMF_history_snapshots:

            # Maybe there weren't any galaxies present for this snapshot.
            if np.isclose(np.sum(model.properties[f"snapshot_{snap}"]["SMF_history"]), 0.0):
                continue

            label = f"{model.label} z = {model.redshifts[snap]:.3f}"
            # The SMF is normalized by the simulation volume which is in Mpc/h.
            ax.plot(
                bin_middles,
                model.properties[f"snapshot_{snap}"]["SMF_history"] / model._volume*pow(model.hubble_h, 3)/bin_widths,
                ls=ls,
                label=label
            )

    # For scaling the observational data, we use the values of the zeroth model.
    zeroth_IMF = models[0].IMF
    ax = obs.plot_temporal_smf_data(ax, zeroth_IMF)

    ax.set_xlabel(r"$\log_{10} M_{\mathrm{stars}}\ (M_{\odot})$")
    ax.set_ylabel(r"$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$")

    ax.set_yscale("log", nonpositive="clip")

    ax.set_xlim([8.0, 12.0])
    ax.set_ylim([1.0e-6, 1.0e-1])

    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

    plot_helper.adjust_legend(ax, location="lower left", scatter_plot=0)

    fig.tight_layout()

    output_file = f"{plot_helper.output_path}A.StellarMassFunction.{plot_helper.output_format}"
    fig.savefig(output_file)
    print(f"Saved file to {output_file}")
    plt.close()

    return fig


def plot_SFRD_history(
    models: List[Model],
    snapshots: List[List[int]],
    plot_helper: PlotHelper,
) -> matplotlib.figure.Figure:
    """
    Plots the evolution of star formation rate density for the specified models.

    Parameters
    ----------
    models : List of :py:class:`~sage_analysis.model.Model` class instance
        Models that will be plotted. These instances contain the properties necessary to create the plot, accessed via
        ``Model.properties["snapshot_<snapshot>"]["property_name"]``.

    snapshots : nested list of ints
        This is a dummy variable that is present to ensure the signature is identical to the other plot functions.

    plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`
        A helper class that contains attributes and methods to assist with plotting. In particular, the path where
        the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
        how to initialize this class and its use.

    Generates
    ---------

    The plot will be saved as "<output_path>B.SFRDensity.<output_format>"
    """

    fig = plt.figure(figsize=plot_helper.figsize)
    ax = fig.add_subplot(111)

    for (model_num, model) in enumerate(models):

        label = model.label
        color = plot_helper.colors[model_num]
        linestyle = plot_helper.linestyles[model_num]
        marker = plot_helper.markers[model_num]

        SFRD = np.array(
            [model.properties[f"snapshot_{snap}"]["SFRD_history"] for snap in model._history_SFRD_history_snapshots]
        )
        redshifts = model._history_SFRD_history_redshifts

        # All snapshots are initialized with zero values for "SFRD_history".  We only want to plot those non-zero
        # values.
        non_zero_inds = np.where(SFRD > 0.0)[0]

        # Only use a line if we have enough snapshots to plot.
        if len(non_zero_inds) > 20:
            ax.plot(
                redshifts[non_zero_inds],
                np.log10(SFRD[non_zero_inds] / model._volume*pow(model.hubble_h, 3)),
                label=label,
                color=color,
                ls=linestyle
            )
        else:
            ax.scatter(
                redshifts[non_zero_inds],
                np.log10(SFRD[non_zero_inds] / model._volume*pow(model.hubble_h, 3)),
                label=label,
                color=color,
                marker=marker,
            )

    ax = obs.plot_sfrd_data(ax)

    ax.set_xlabel(r"$\mathrm{redshift}$")
    ax.set_ylabel(r"$\log_{10} \mathrm{SFR\ density}\ (M_{\odot}\ \mathrm{yr}^{-1}\ \mathrm{Mpc}^{-3})$")

    ax.set_xlim([0.0, 8.0])
    ax.set_ylim([-3.0, -0.4])

    ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))

    plot_helper.adjust_legend(ax, location="lower left", scatter_plot=0)

    fig.tight_layout()

    output_file = f"{plot_helper.output_path}B.SFRDensity.{plot_helper.output_format}"
    fig.savefig(output_file)
    print(f"Saved file to {output_file}")
    plt.close()

    return fig


def plot_SMD_history(
    models: List[Model],
    snapshots: List[List[int]],
    plot_helper: PlotHelper,
) -> matplotlib.figure.Figure:
    """
    Plots the evolution of stellar mass density for the specified models.

    Parameters
    ----------
    models : List of :py:class:`~sage_analysis.model.Model` class instance
        Models that will be plotted. These instances contain the properties necessary to create the plot, accessed via
        ``Model.properties["snapshot_<snapshot>"]["property_name"]``.

    snapshots : nested list of ints
        This is a dummy variable that is present to ensure the signature is identical to the other plot functions.

    plot_helper : :py:class:`~sage_analysis.plot_helper.PlotHelper`
        A helper class that contains attributes and methods to assist with plotting. In particular, the path where
        the plots will be saved and the output format.  Refer to :doc:`../user/plot_helper` for more information on
        how to initialize this class and its use.

    Generates
    ---------

    The plot will be saved as "<output_path>C.StellarMassDensity.<output_format>"
    """

    fig = plt.figure(figsize=plot_helper.figsize)
    ax = fig.add_subplot(111)

    for (model_num, model) in enumerate(models):

        label = model.label
        color = plot_helper.colors[model_num]
        linestyle = plot_helper.linestyles[model_num]
        marker = plot_helper.markers[model_num]

        SMD = np.array(
            [model.properties[f"snapshot_{snap}"]["SMD_history"] for snap in model._history_SMD_history_snapshots]
        )
        redshifts = model._history_SMD_history_redshifts

        # All snapshots are initialized with zero values for "SMD_history".  We only want to plot those non-zero
        # values.
        non_zero_inds = np.where(SMD > 0.0)[0]

        # Only use a line if we have enough snapshots to plot.
        if len(non_zero_inds) > 20:
            ax.plot(
                redshifts[non_zero_inds],
                np.log10(SMD[non_zero_inds] / model._volume*pow(model.hubble_h, 3)),
                label=label,
                color=color,
                ls=linestyle
            )
        else:
            ax.scatter(
                redshifts[non_zero_inds],
                np.log10(SMD[non_zero_inds] / model._volume*pow(model.hubble_h, 3)),
                label=label,
                color=color,
                marker=marker,
            )

    # For scaling the observational data, we use the values of the zeroth
    # model.
    zeroth_IMF = models[0].IMF
    ax = obs.plot_smd_data(ax, zeroth_IMF)

    ax.set_xlabel(r"$\mathrm{redshift}$")
    ax.set_ylabel(r'$\log_{10}\ \phi\ (M_{\odot}\ \mathrm{Mpc}^{-3})$')

    ax.set_xlim([0.0, 4.2])
    ax.set_ylim([6.5, 9.0])

    ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))

    plot_helper.adjust_legend(ax, location="lower left", scatter_plot=0)

    fig.tight_layout()

    output_file = f"{plot_helper.output_path}C.StellarMassDensity.{plot_helper.output_format}"
    fig.savefig(output_file)
    print(f"Saved file to {output_file}")
    plt.close()

    return fig
