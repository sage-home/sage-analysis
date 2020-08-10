from sage_analysis.utils import generate_func_dict

default_plot_toggles = {
    "SMF" : True,  # Stellar mass function.
    "BMF" : True,  # Baryonic mass function.
    "GMF" : True,  # Gas mass function (cold gas).
    "BTF" : True,  # Baryonic Tully-Fisher.
    "sSFR" : True,  # Specific star formation rate.
    "gas_fraction" : True,  # Fraction of galaxy that is cold gas.
    "metallicity" : True,  # Metallicity scatter plot.
    "bh_bulge" : True,  # Black hole-bulge relationship.
    "quiescent" : True,  # Fraction of galaxies that are quiescent.
    "bulge_fraction" : True,  # Fraction of galaxies that are bulge/disc dominated.
    "baryon_fraction" : True,  # Fraction of baryons in galaxy/reservoir.
    "reservoirs" : True,  # Mass in each reservoir.
    "spatial" : True,   # Spatial distribution of galaxies.
    "SMF_history": False,
    "SFRD_history": False,
    "SMD_history": False,
}

default_galaxy_properties_to_analyze = {
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
    },
    "halo_mass_bins": {
        "type": "binned",
        "bin_low": 10.0,
        "bin_high": 14.0,
        "bin_width": 0.1,
        "property_names": ["fof_HMF"] + [f"halo_{component}_fraction_sum"
            for component in ["baryon", "stars", "cold", "hot", "ejected", "ICS", "bh"]
        ],
    },
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
    },
    "single_properties": {
        "type": "single",
        "property_names": ["SMD_history", "SFRD_history"],
    },
}


default_calculation_functions = generate_func_dict(default_plot_toggles, "sage_analysis.example_calcs", "calc_")
default_plot_functions = generate_func_dict(default_plot_toggles, "sage_analysis.example_calcs", "calc_")
