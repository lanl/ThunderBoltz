"""Create specific plots relevant to the ThunderBoltz paper."""
import itertools
import os
import sys
from os.path import join as pjoin

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import thunderboltz as tb
from thunderboltz import read
from thunderboltz import query_tree
from thunderboltz.plotting import pdplot
from thunderboltz.plotting import styles

def plot_He_transport():
    """Plot bulk and flux transport parameters vs. reduced electric field."""
    name = "He_transport"

    # Import data
    dat = load_transport_data()

    # independent styling
    c1, c2, c3, c4 = "red", "blue", "orange", "#0aa"
    ss = {
        "Bulk ThunderBoltz (One Takes All)": {"ls": "", "c": c1, "marker": "v", "markersize": 11},
        "Flux ThunderBoltz (One Takes All)": {"ls": "", "c": c2, "marker": "v", "markersize": 11},
        "Bulk Bolsig (One Takes All)": {"ls": "--", "c": c1},
        "Flux Bolsig (One Takes All)": {"ls": "--", "c": c2},
        "Bulk ThunderBoltz (Equal)": {"ls": "", "c": c3, "marker": "v", "markersize": 11},
        "Flux ThunderBoltz (Equal)": {"ls": "", "c": c4, "marker": "v", "markersize": 11},
        "Bulk Bolsig (Equal)": {"ls": "--", "c": c3},
        "Flux Bolsig (Equal)": {"ls": "--", "c": c4},
    }
    exp_keys = dat[~dat.name.isin(["Bolsig", "ThunderBoltz"])].name.unique()
    ss.update({k: {"markersize": 6, "marker": "o", "ls": "", "alpha": 0.9, "c": c}
               for k, c in zip(exp_keys, ["red", "green", "orange", "blue",
                                          "#afa", "#aa0852", "#dba", "#abd"])})

    y = ["a_n_-21", "mobN"]
    x = "Ered"

    errs = {"mobN": "mobN_std", "a_n-21": "a_n-21_std"}

    bmap = {True: "Bulk", False: "Flux"}
    def create_series_labels(row):
        if not pd.isnull(row.eesd):
            return f"{bmap[row.bulk]} {row['name']} ({row.eesd})"
        return row["name"]

    dat["series"] = dat.apply(create_series_labels, axis=1)
    legend_ordering = [7, 8, 0, 1, 9, 10, 2, 3, 4, 5, 6, 11, 12, 13, 14]
    x0, x_width = 0.52, 0.48
    inset_bounds = [[[x0, .09, x_width, .5]], [None]]

    p = tb.plotting.pdplot.PandaPlot(dat, x=x, y=y, err_bars=errs, series="series", remove_ext=True,
        legend_style="figure bottom", xscale="log", yscale="log", label_map=styles.label_maps,
        value_map=styles.value_maps, series_styles=ss, gui_configure=False, save_spacing_config=True,
        load_spacing_config=False, legend_ordering=legend_ordering, config_dir="simulations", inset_bounds=inset_bounds,
        name="He_transport", snapx=True)

    p.fig.subplots_adjust(left=0.162, bottom=.41, top=0.99, hspace=0.05)
    p.fig.set_size_inches(11/1.49, 24/1.49)
    p.axs[0][0].xaxis.set_tick_params(labelbottom=False)
    p.axs[1][0].set_xlim(.9, 1800)
    p.axs[0][0].set_xlim(.9, 1800)
    p.axs[1][0].set_ylim(1.5, 7)
    p.axs[0][0].set_ylim(.0001, 10)
    p.axs[1][0].set_yscale("linear")

    x1 = invert_coords(p.axs[0][0], x0, 0)[0]
    x2 = invert_coords(p.axs[0][0], x0+x_width, 0)[0]
    p.inset_axs[0][0].set_xlim(x1, x2)
    p.inset_axs[0][0].set_ylim(.1, 6)

    p.save("simulations")
    os.system("open simulations/He_transport.pdf")

def rate_comp():
    name = "rate_comp"
    Ereds = [50, 100, 500]
    # Import data
    tb_dat = load_transport_data()
    tb_dat = tb_dat[tb_dat.name == "ThunderBoltz"].copy()

    rate_keys = [k for k in tb_dat.columns if "k_" in k and k not in ["k_i", "k_tot"] and len(k) < 5]
    rate_std = [k+"_std" for k in rate_keys]
    all_keys = ["name", "growth", "Ered"]+rate_keys+rate_std


    data_store = pd.read_csv(pjoin("simulations", "He_rate_data.csv")) 
    bols_dat = data_store[data_store.name == "Bolsig"].copy()

    dat = pd.concat((bols_dat, tb_dat), ignore_index=True)

    dat = dat[(
        dat.Ered.isin(Ereds)
      & ~dat.k_1.isna()
      & (dat.egen!=False)
      & (dat.eesd!="Equal")
      & ~dat.growth.isin([1])
    )].copy()


    id_keys = [k for k in dat.columns if k not in rate_keys]
    m1 = pd.melt(dat, id_keys, rate_keys, var_name="Reaction #", value_name="k")
    id_keys = [k for k in dat.columns if k not in rate_std]
    m2 = pd.melt(dat, id_keys, rate_std, var_name="Reaction #", value_name="k_std")
    m1.loc[:, "k_std"] = m2.k_std
    dat = m1.copy()

    dat.loc[:, "Reaction #"] = dat["Reaction #"].map(lambda s: int(s.split("_")[1]))

    # Create legend keys
    ML = [("name", {"Bolsig": ["growth", "Ered"],
                    "ThunderBoltz": ["Ered"]})]
    # independent styling
    ss = {"Ered": [Ereds, "c", styles.DCC[:3]],
          # "name": [["Bolsig", "ThunderBoltz"], "ls", ["-", "--"]],
          "growth": [[4.0, 2.0], "ls", ["-", (0, (3, 1, 1, 1))]],
    }

    x = "Reaction #"
    y = ["k"]
    errs = {"k": "k_std"}

    p = pdplot.PandaPlot(dat, x=x, y=y, err_bars=errs, series=ML, multi_legend=True,
        yscale="log", label_map=styles.label_maps, value_map=styles.value_maps,
        legend_style="figure right", series_styles=ss, name=name, snapx=True)

    p.fig.set_size_inches(8.8, 8)
    p.fig.subplots_adjust(left=.17, right=0.9, top=.99)

    # Move figure legend into blank space
    p.figure_legend.set_bbox_to_anchor([0.53, 0.68, 0.3, 0.27])

    p.save("simulations")

    os.system("open simulations/rate_comp.pdf")

def N2_cs_resolve():
    """Plot the cross section features for the N2 gas."""
    name = "N2_cs_resolve"
    # Load Tinity data from LXCat output txt file
    fig, ax = util.plot_LXCat_cs("data/Trinity_N2.txt")

    x1 = 1.59252
    x2 = 3.50687
    y = 7e-20
    # Draw vertical lines
    ylims = ax.get_ylim()
    ax.plot([x1]*3, [1e-20, y, 1e-19], c="black", lw=3)
    ax.plot([x2]*3, [1e-20, y, 1e-19], c="black", lw=3)
    ax.set_ylim(*ylims)
    # Draw connecting arrow
    ax.annotate("", xy=(x1, y),
        textcoords="data", xytext=(x2, y),
        fontsize=17, arrowprops={"arrowstyle": "<|-|>", "ls": "-",
            "lw": 3, "color": "black"})
    # Write symbols
    p = (np.sqrt((x1*x2))-.4, y+1e-20)
    ax.annotate("$2\Delta \epsilon$", xy=p, textcoords="data", xytext=p, fontsize=18)

    p = (.7*x1, 2.5e-20)
    ax.annotate("$\epsilon_0$", xy=p, textcoords="data", xytext=p, fontsize=18)

    ax.set_xlim(2e-2, 9e2)

    # Save to pdf
    fig.savefig("figures/N2_cs_resolve.pdf")

def plot_onsager():
    """Generate kinetic plots to compare to the Onsager relation."""
    # Read thunderboltz calculation directory
    calc = read(pjoin("simulations", "onsager_relation"))
    ts = calc.get_timeseries()
    tab = calc.cs.table

    # This is the species ordering
    sp = ["X", "A", "B", "C"]

    # Grab the densities
    n0 = ts["Ni"]
    ns = pd.concat([calc.particle_tables[i]["Ni"]
                     for i in range(1, 4)], axis=1)
    # And the energies
    es = pd.concat([calc.particle_tables[i]["Mi"]
                     for i in range(0, 4)], axis=1)
    # Rename columns
    ns.columns = ["$n_A$", "$n_B$", "$n_C$"]
    es.columns = [r"$\langle{E\rangle}" + f"_{l}$" for l in "XABC"]

    # Create figure
    fig, axs = plt.subplots(4, 1, figsize=(7.4, 11.8))
    axs[0].set_ylabel("Density (10$^{22}$/m$^3$)")
    axs[1].set_ylabel("$k_{ij}$ (10$^{-16}$m$^3$/s)")
    axs[2].set_ylabel("$n_ik_{ij}$ (10$^6$/s)")
    axs[3].set_ylabel(r"$\langle{E\rangle}$ (eV)")
    axs[3].set_xlabel("Time ($\mu$s)")

    # Plot densities
    for n_key in ns:
        axs[0].plot(ts.t*1e6, ns[n_key]*1e-22, label=n_key)
        axs[0].legend(loc="upper right", fontsize=12)


    # Plot energies
    # Weighted average over all species
    av = (n0*es.iloc[:,0] + (ns*es.iloc[:,1:]).sum(axis=1)) / (n0+ns.sum(axis=1))
    for e_key in es.iloc[:,1:]:
        axs[3].plot(ts.t*1e6, es[e_key], label=e_key)
    axs[3].plot(ts.t*1e6, es.iloc[:,0], c="black", label=r"$\langle{E\rangle}_{X}$")
    axs[3].legend(loc="upper right", fontsize=12)

    # Select the processes with a particle type change
    reacting = tab[tab.p2 != tab.r2].copy()

    style_cycle = itertools.cycle([{}, {"lw": 1, "marker": "x"}])
    colors = 2*["#1f77b4"] + 2*["#ff7f0e"] + 2*["#2ca02c"]

    # Generate the theoretical rates
    AMU_TO_KG = 1.6605e-27 # kg / amu
    QE = 1.60217663e-19 # C or J/eV
    mr = AMU_TO_KG * (14. / 2)
    kbT = 1.
    # the rate formula in terms of cross section area, activation energy
    kf = lambda A, Ea: A*np.sqrt(8*kbT*QE/(np.pi*mr)) * np.exp(-Ea/kbT)
    kr = lambda A: A*np.sqrt(8*kbT*QE/(np.pi*mr))
    analytic_ks = {
        (1, 2): kf(1e-20, 1.),
        (2, 1): kr(1e-20),
        (2, 3): kf(2e-20, 1.),
        (3, 2): kr(2e-20),
        (1, 3): kf(3e-20, 2.),
        (3, 1): kr(3e-20),
    }

    # Plot the rate coefficients and the absolute rates.
    for (i, row), c in zip(reacting.iterrows(), colors):
        # Get theoretical forward rate
        k_th = analytic_ks[(row.r2, row.p2)]
        style = next(style_cycle)
        # Plot rates
        label = f"$k_{{{sp[row.r2]}{sp[row.p2]}}}$"
        axs[1].plot(ts.t*1e6, ts[f"k_{i+1}"]*1e16, label=label, c=c, **style)
        axs[1].plot(ts.t*1e6, len(ts)*[k_th*1e16], c="black", ls=(0, (3, 1, 1, 1)), lw=2)

        # Plot absolute rates
        K_key = f"$k_{{{sp[row.r2]}{sp[row.p2]}}}$"
        n_i = f"$n_{sp[row.r2]}$"
        axs[2].plot(ts.t*1e6, ns[n_i]*ts[f"k_{i+1}"]*1e-6, label=n_i+label, c=c, **style)

    axs[1].plot([], [], c="black", ls=(0, (3, 1, 1, 1)), lw=2, label="benchmark")
    axs[2].legend(loc="upper right", fontsize=10)
    axs[1].legend(loc="upper right", fontsize=8.35)

    fig.subplots_adjust(top=0.97, bottom=0.08, left=0.16, )
    fig.savefig(pjoin("simulations", "onsager_relation.pdf"))

def plot_ikuta_sugai():
    """Generate a velocity moment profile at various
    crossed reduced E/B-fields"""
    # Read output files from simulation base
    calcs = query_tree(pjoin("simulations", "ikuta_sugai"))
    # Extract steady state parameters
    ss = pd.concat([calc.get_ss_params() for calc in calcs], ignore_index=True)
    # Only need B/n, Vxi, Vzi
    ss = ss[["Bred", "Vxi", "Vzi", "MEe"]].copy()
    # Reformat reduced values
    ss.loc[:, "Bred"] = ss.Bred.map(lambda s: float(s.split()[1]))
    ss = ss.sort_values("Bred")

    # Import benchmark data from Ness
    b = pd.read_csv(
        pjoin("simulations", "Ikuta_Sugai_comparison.csv"))

    fig, ax = plt.subplots(figsize=(8, 8))

    ax.set_xlabel("B/n (Hx)")
    l1, = ax.plot(ss.Bred, ss.Vxi*1e-4, ls="", markersize=15, marker="v", label="V$_x$ (10$^4$m/s)")
    l2, = ax.plot(ss.Bred, ss.Vzi*1e-4, ls="", markersize=15, marker="v", label="V$_z$ (10$^4$m/s)")
    l3, = ax.plot(ss.Bred, ss.MEe, ls="", markersize=15, marker="v", label=r"$\langle \epsilon \rangle$ (eV)")

    # Plot benchmark data
    ax.plot(b["B/n"], b.Wx, ls="", markersize=12, marker="+", c="black")
    ax.plot(b["B/n"], b.Wz, ls="", markersize=12, marker="+", c="black")
    ax.plot(b["B/n"], b.e_mean, ls="", markersize=12, marker="+", c="black")
    ax.plot([], [], markersize=12, ls="", marker="+", c="black", label="Ness (benchmark)")
    ax.legend()

    ax.grid(which='major', color='gray', linestyle='-', alpha=0.5)
    ax.grid(which='minor', color='gray', linestyle='-', alpha=0.3)

    fig.subplots_adjust(bottom=0.16)

    fig.savefig(pjoin("simulations", "ikuta_sugai.pdf"))

def view_all_plots():
    # Read output files from simulation base
    calcs = query_tree(pjoin("simulations", "He_transport"))
    calc = calcs[-1]
    calc.plot_timeseries()
    calc.plot_vdfs()
    calc.plot_edfs()
    calc.plot_edf_comps()

    plt.show()


# Helpers

def melt_bulk(dat):
    """Melt the bulk and flux type parameters into separate rows."""
    m = pd.DataFrame()
    bm = {"mobN": "mobN_bulk", "mobN_std": "mobN_bulk_std",
          "a_n": "a_n_bulk", "a_n_std": "a_n_bulk_std"}
    ms = []
    for tp, tpb in bm.items():
        id_vars = [k for k in dat if k not in [tp,tpb]]
        ms.append(pd.melt(dat, id_vars, [tp,tpb], var_name="bulk",
            value_name=tp+"_VALUE"))
        ms[-1] = ms[-1].rename(columns={tp+"_VALUE": tp})
        ms[0].loc[:,tp] = ms[-1][tp]

    dat = ms[0].copy()
    dat.loc[:, "bulk"] = dat.bulk.str.contains("_bulk")
    return dat

def load_tb(csvfile):
    """Format ThunderBoltz Data."""
    df = pd.read_csv(csvfile)
    df = df.rename(columns={"Ered (Td)": "Ered", "EmobN (V^-1m^-1s^-1)": "mobN",
                            "EmobN (V^-1m^-1s^-1)_std": "mobN_std"})
    df.loc[:, "eesd"] = df.eesd.map({"uniform": "Uniform", "default": "One Takes All",
                                     "equal": "Equal"})
    df.loc[:, "eadf"] = df.eadf.map({"default": "Isotropic", "He_Park": "Anisotropic"})
    df["name"] = "ThunderBoltz"
    return df

def invert_coords(ax, x, y):
    disp = ax.transAxes.transform((x, y))
    return ax.transData.inverted().transform(disp)

def load_transport_data():
    # Read simulation output data
    dat = tb.query_tree(
        "simulations/He_transport", # Read all the calculations in this directory
        callback=lambda calc: calc.get_ss_params() # return their steady state values.
    )
    dat = dat[dat.eesd.isin(["default", "equal"])].copy()

    # Format model names
    dat.loc[:, "eesd"] = dat.eesd.map({"uniform": "Uniform", "default": "One Takes All",
                                     "equal": "Equal"})
    dat.loc[:, "eadf"] = dat.eadf.map({"default": "Isotropic", "He_Park": "Anisotropic"})

    # Separate bulk data from flux data.
    dat = melt_bulk(dat)
    dat["name"] = "ThunderBoltz"

    # Load in comparison data
    comp = pd.read_csv("simulations/He_transport_comparison.csv")
    comp = comp[~comp.name.str.contains("Takeda")].reset_index(drop=True).copy()
    comp["data_class"] = "experimental"
    dat = pd.concat((dat, comp), ignore_index=True)
    dat["a_n_-21"] = 1e21*dat.a_n
    dat["a_n_-21_std"] = 1e21*dat.a_n_std
    dat["mobN"]= 1e-24*dat.mobN
    dat["mobN_std"]= 1e-24*dat.mobN_std
    return dat

def format_uncertainty(x, std):
    s = f"{x:.2e}".replace("e", "E")
    if std == 0:
        return s
    i = s.index("E")
    scal = -int(s[i+1:])
    s = s[:i] + f"({std*10**scal:.3f})" + s[i:]
    return s
