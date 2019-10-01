#!/usr/bin/env python
"""
Here we show a myriad of functions that can be used to plot properties calculated from the
**SAGE** output across multiple snapshots.

We refer to for more information on how plotting is handled :doc:`here <../docs/source/user/history>`.

Author: Jacob Seiler
"""

import matplotlib
from matplotlib import pyplot as plt
import numpy as np

import sage_analysis.observations as obs

colors = ["r", "b", "g", "m", "c"]
linestyles = ["-", ":", "--", "-.", "-:"]
markers = ["x", "o", "v", "*", "D"]


def setup_matplotlib_options():
    """
    Set the default plotting parameters.
    """

    matplotlib.rcdefaults()
    plt.rc('xtick', labelsize='x-large')
    plt.rc('ytick', labelsize='x-large')
    plt.rc('lines', linewidth='2.0')
    plt.rc('legend', numpoints=1, fontsize='x-large')


def adjust_legend(ax, location="upper right", scatter_plot=0):
    """
    Adjusts the legend of a specified axis.

    Parameters
    ----------

    ax : ``matplotlib`` axes object
        The axis whose legend we're adjusting

    location : String, default "upper right". See ``matplotlib`` docs for full options
        Location for the legend to be placed.

    scatter_plot : {0, 1}
        For plots involved scattered-plotted data, we adjust the size and alpha of the
        legend points.

    Returns
    -------

    None. The legend is placed directly onto the axis.
    """

    legend = ax.legend(loc=location)
    handles = legend.legendHandles

    legend.draw_frame(False)

    # First adjust the text sizes.
    for t in legend.get_texts():
        t.set_fontsize("medium")

    # For scatter plots, we want to increase the marker size.
    if scatter_plot:
        for handle in handles:
            # We may have lines in the legend which we don't want to touch here.
            if isinstance(handle, matplotlib.collections.PathCollection):
                handle.set_alpha(1.0)
                handle.set_sizes([10.0])




def plot_SMF(models, plot_output_path, plot_output_format=".png"):
    """
    Plots the evolution of the stellar mass function for the specified models.
    This function loops over the value of ``model.SMF_snaps`` and plots and the SMFs at
    each snapshots.

    Parameters
    ----------

    models : List of ``Model`` class instance
        Models that will be plotted. These instances contain the properties necessary to
        create the plot, accessed via ``Model.properties["property_name"]``. In
        particular, we acces the ``Model.properties["SMF_dict"][<snap>]`` values.

    plot_output_path : string
        Path to where the plot will be saved.

    plot_output_format : string, default ".png"
        Format the plot will be saved in, includes the full stop.

    Generates
    ---------

    The plot will be saved as "<plot_output_path>/A.StellarMassFunction<plot_output_format>"
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Go through each of the models and plot.
    for (model_num, model) in enumerate(models):

        ls = linestyles[model_num]

        # Set the x-axis values to be the centre of the bins.
        bin_widths = model.bins["stellar_mass_bins"][1::] - model.bins["stellar_mass_bins"][0:-1]
        bin_middles = model.bins["stellar_mass_bins"][:-1] + bin_widths

        # Iterate over the snapshots.
        for snap in model.SMF_snaps:
            label = "{0} z = {1:.3f}".format(model.label, model.redshifts[snap])

            # The SMF is normalized by the simulation volume which is in Mpc/h.
            ax.plot(bin_middles, model.properties["SMF_dict"][snap] / model.volume*pow(model.hubble_h, 3)/bin_widths,
                    ls=ls, label=label)

    # For scaling the observational data, we use the values of the zeroth
    # model.
    zeroth_IMF = models[0].IMF
    ax = obs.plot_temporal_smf_data(ax, zeroth_IMF)

    ax.set_xlabel(r"$\log_{10} M_{\mathrm{stars}}\ (M_{\odot})$")
    ax.set_ylabel(r"$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$")

    ax.set_yscale("log", nonposy="clip")

    ax.set_xlim([8.0, 12.0])
    ax.set_ylim([1.0e-6, 1.0e-1])

    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

    adjust_legend(ax, location="lower left", scatter_plot=0)

    fig.tight_layout()

    output_file = "{0}/A.StellarMassFunction.{1}".format(plot_output_path,
                                                         plot_output_format)
    fig.savefig(output_file)
    print("Saved file to {0}".format(output_file))
    plt.close()


def plot_SFRD(models, plot_output_path, plot_output_format=".png"):
    """
    Plots the evolution of star formation rate density for the specified models.

    Parameters
    ----------

    models : List of ``Model`` class instance
        Models that will be plotted. These instances contain the properties necessary to
        create the plot, accessed via ``Model.properties["property_name"]``. In
        particular, we acces the ``Model.properties["SFRD_dict"][<snap>]`` values.

    plot_output_path : string
        Path to where the plot will be saved.

    plot_output_format : string, default ".png"
        Format the plot will be saved in, includes the full stop.

    Generates
    ---------

    The plot will be saved as "<plot_output_path>/B.SFRDensity<plot_output_format>"
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for (model_num, model) in enumerate(models):

        label = model.label
        color = colors[model_num]
        linestyle = linestyles[model_num]

        # The SFRD is in a dictionary. Pull it out into a array for plotting.
        SFRD = np.array([model.properties["SFRD_dict"][snap] for snap in model.properties["SFRD_dict"].keys()])
        ax.plot(model.redshifts[model.density_snaps], np.log10(SFRD / model.volume*pow(model.hubble_h, 3)),
                label=label, color=color, ls=linestyle)

    ax = obs.plot_sfrd_data(ax)

    ax.set_xlabel(r"$\mathrm{redshift}$")
    ax.set_ylabel(r"$\log_{10} \mathrm{SFR\ density}\ (M_{\odot}\ \mathrm{yr}^{-1}\ \mathrm{Mpc}^{-3})$")

    ax.set_xlim([0.0, 8.0])
    ax.set_ylim([-3.0, -0.4])

    ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))

    adjust_legend(ax, location="lower left", scatter_plot=0)

    fig.tight_layout()

    output_file = "{0}/B.SFRDensity.{1}".format(plot_output_path, plot_output_format)
    fig.savefig(output_file)
    print("Saved file to {0}".format(output_file))
    plt.close()


def plot_SMD(models, plot_output_path, plot_output_format=".png"):
    """
    Plots the evolution of stellar mass density for the specified models.

    Parameters
    ----------

    models : List of ``Model`` class instance
        Models that will be plotted. These instances contain the properties necessary to
        create the plot, accessed via ``Model.properties["property_name"]``. In
        particular, we acces the ``Model.properties["SMD_dict"][<snap>]`` values.

    plot_output_path : string
        Path to where the plot will be saved.

    plot_output_format : string, default ".png"
        Format the plot will be saved in, includes the full stop.

    Generates
    ---------

    The plot will be saved as "<plot_output_path>/C.StellarMassDensity<plot_output_format>"
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for (model_num, model) in enumerate(models):

        label = model.label
        color = colors[model_num]
        linestyle = linestyles[model_num]

        # The SMD is in a dictionary. Pull it out into a array for plotting.
        SMD = np.array([model.properties["SMD_dict"][snap] for snap in model.properties["SMD_dict"].keys()])
        ax.plot(model.redshifts[model.density_snaps], np.log10(SMD / model.volume * pow(model.hubble_h, 3)),
                label=label, color=color, ls=linestyle)

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

    adjust_legend(ax, location="lower left", scatter_plot=0)

    fig.tight_layout()

    output_file = "{0}/C.StellarMassDensity.{1}".format(plot_output_path, plot_output_format)
    fig.savefig(output_file)
    print("Saved file to {0}".format(output_file))
    plt.close()
