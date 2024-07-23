from collections import OrderedDict
import copy
import itertools
import os
from os.path import join as pjoin
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from numpy import vectorize as vec
import seaborn as sns
import scipy

from thunderboltz import parameters
from thunderboltz import parsing
from thunderboltz.plotting import pdplot
from thunderboltz.plotting import styles

ME = 9.1093837e-31 # kg
QE = 1.60217663e-19 # C or J/eV
AMU_TO_KG = 1.6605e-27 # kg / amu

def plot_timeseries(fig, ts, series=["MEe", "mobN", "a_n"], save=None):
    """Plot ThunderBoltz output values as function of time."""
    if len(fig.axes) > 1:
        # Reset the twin axes
        [a.remove() for a in fig.axes[1:]]

    # Map strings to callables
    names = [s[0] if isinstance(s, tuple) else s for s in series]
    user_styles = [s[2] if isinstance(s, tuple) and len(s)>2 else {}
                   for s in series]
    series = [s[1](ts) if isinstance(s, tuple) else ts[s] for s in series]

    # Create twin axes
    [fig.axes[0].twinx() for param in series[1:]]
    ax = fig.axes[0]
    ax.grid(which='major', color='gray', linestyle='-', alpha=0.5)
    ax.grid(which='minor', color='gray', linestyle='-', alpha=0.3)
    # Plot first dataset on left hand side
    ax.plot(ts.t, series[0], label=styles.label_maps(names[0]), c="green", **user_styles[0])

    # Get lims for before plotting verticle lines
    ylm = ax.get_ylim()
    ax.set_ylim(*ylm)
    ax.set_ylabel(styles.label_maps(names[0]))
    ax.set_xlabel("Time (s)")
    colors = ["purple", "red", "orange", "pink", "green", "blue", "#afa", "#33b"]
    offset = 1.
    for i, (s, n, u, c) in enumerate(zip(series[1:], names[1:], user_styles[1:], colors)):
        axtwin = fig.axes[i+1]
        axtwin.spines.right.set_position(("axes", offset))
        offset += 0.18
        axtwin.set_ylabel(styles.label_maps(n))
        axtwin.plot(ts.t, s, c=c, label=styles.label_maps(n), **u)

    handles, labels = figure_legend(fig.axes)
    leg = fig.legend(handles, labels, loc="lower center", ncol=100)
    fig.subplots_adjust(right=max(.2, 2-.9*offset), bottom=.18, top=.80)
    return fig

def plot_eedfs(ax, esamples, bins=150):
    """Plot the electron energy density funciton"""
    pass

def plot_vdfs(vdfs, pair=["vx", "vy"], bins=50):
    """Plot the electron velocity distribution function."""
    # Seaborn color palette
    rocket = sns.color_palette("rocket_r", as_cmap=True)
    # Create joint plot
    jplot = sns.jointplot(x=vdfs[pair[0]]/1e6, y=vdfs[pair[1]]/1e6, kind="hist",
            cbar=True, marginal_ticks=False, cmap=rocket, bins=bins,
            marginal_kws={"color": (209/255, 38/255, 65/255), "bins": bins})
    plt.subplots_adjust(left=0.1, right=0.8, top=0.95, bottom=0.1)
    # Get the positions of the joint and marginal density plots
    jpos = jplot.ax_joint.get_position()
    xmpos = jplot.ax_marg_x.get_position()
    ympos = jplot.ax_marg_y.get_position()
    jplot.ax_joint.set_position([jpos.x0+.02, jpos.y0, xmpos.width, jpos.height])
    jplot.ax_marg_x.set_position([jpos.x0+.02, jpos.y1, xmpos.width, xmpos.height])
    # Move colorbar to the outside
    cbar = jplot.fig.axes[-1]
    cbar.set_position([.83, jpos.y0, .07, jpos.height])
    jax = jplot.fig.axes[0]
    jax.set_xlabel(f"$v_{pair[0][1]}$ (10$^6$ m/s)")
    jax.set_ylabel(f"$v_{pair[1][1]}$ (10$^6$ m/s)")
    jplot.fig.set_size_inches(10, 8.9)
    return jplot.fig

def align_vdfs_on(tbs, params, save=None, bins=100):
    """Align vdfs figures on every unique combination of `params`.

    tbs: list<ThunderBoltz>
        List of ThunderBoltz objects to create VDF plots for.
    params: list | str
        ThunderBoltz parameter key(s) that form groups to align scales on."""
    if isinstance(params, str):
        params = [params]
    def get_val(tb, p):
        tbp = tb.tb_params
        if p in tbp:
            return tbp[p]
        if p in tb.hp:
            return tb.hp[p]
    # Make table of unique sets
    param_table = pd.DataFrame([dict(tb=tb, **{p: get_val(tb, p) for p in params})
                  for tb in tbs])
    # Collect all the figure objects
    figs = []
    # Where to save them, if necessary
    directories = []
    for ids, tb_group in param_table.groupby(params):
        # Get directories
        directories = [tb.get_directory() for tb in tb_group["tb"]]
        # Plot and collect figures for each of the tb objects
        # TODO: implement bins for step grouping
        group_figs = [tb.plot_vdfs(steps="last", save=False, bins=bins)
                      for tb in tb_group["tb"].values]
        xy_figs = [f[0][0] for f in group_figs if f]
        xz_figs = [f[0][1] for f in group_figs if f]
        steps = [f[1][0] for f in group_figs if f]

        # If there are no figures available in this group, continue
        if not steps:
            continue

        # Get the axes and colorbar  limits for each plot
        xy_xlims = [f.axes[0].get_xlim() for f in xy_figs]
        xy_ylims = [f.axes[0].get_ylim() for f in xy_figs]
        xy_clims = [f.axes[0].collections[0].get_clim() for f in xy_figs]
        xz_xlims = [f.axes[0].get_xlim() for f in xz_figs]
        xz_ylims = [f.axes[0].get_ylim() for f in xz_figs]
        xz_clims = [f.axes[0].collections[0].get_clim() for f in xz_figs]

        # Get max lims within group for each pair
        xy_m_x = (min(l[0] for l in xy_xlims), max(l[1] for l in xy_xlims))
        xy_m_y = (min(l[0] for l in xy_ylims), max(l[1] for l in xy_ylims))
        xy_m_c = (min(l[0] for l in xy_clims), max(l[1] for l in xy_clims))
        xz_m_x = (min(l[0] for l in xz_xlims), max(l[1] for l in xz_xlims))
        xz_m_y = (min(l[0] for l in xz_ylims), max(l[1] for l in xz_ylims))
        xz_m_c = (min(l[0] for l in xz_clims), max(l[1] for l in xz_clims))

        # Set each of the lims the these extrema
        for fxy, fxz, step, directory in zip(xy_figs, xz_figs, steps, directories):
            figs.extend([fxy, fxz])
            fxy.axes[0].set_xlim(*xy_m_x)
            fxy.axes[0].set_ylim(*xy_m_y)
            fxy.axes[0].collections[0].set_clim(*xy_m_c)
            fxz.axes[0].set_xlim(*xz_m_x)
            fxz.axes[0].set_ylim(*xz_m_y)
            fxz.axes[0].collections[0].set_clim(*xz_m_c)
            if save:
                sname = Path(directory).parts[-1]
                if not os.path.isdir(save):
                    os.makedirs(save)
                saveas = pjoin(save, sname+f"_xy{step}.pdf")
                print(f"Saving VDF joint plot in {saveas}", flush=True)
                fxy.savefig(saveas)
                plt.close(fxy)
                saveas = pjoin(save, sname+f"_xz{step}.pdf")
                print(f"Saving VDF joint plot in {saveas}", flush=True)
                fxz.savefig(saveas)
                plt.close(fxz)
    return figs

def vdf_pct_err(tbs, params, benchmark, save=None, bins=100):
    """Calculate the % error of the velocity distribution density function.

    tbs: list<ThunderBoltz>
        List of ThunderBoltz objects to create VDF plots for.
    params: list | str
        ThunderBoltz parameters which are being compared in the error calculation.
    benchmark: list | str
        Value of `params` to compare all other calculations against. Must be
        same length as params.
    """

    if isinstance(params, str):
        params = [params]
    def get_val(tb, p):
        tbp = tb.tb_params
        if p in tbp:
            return tbp[p]
        if p in tb.hp:
            return tb.hp[p]
    # Make table of unique sets
    # print(*[dict(tb=tb, **{p: get_val(tb, p) for p in params})
                  # for tb in tbs], sep="\n\n")

    param_table = pd.DataFrame([dict(tb=tb, **{p: get_val(tb, p) for p in params})
                  for tb in tbs])
    # Collect all the figure objects
    figs = []
    # Where to save them, if necessary
    directories = []
    for ids, tb_group in param_table.groupby(params):
        # Get directories
        directories = [tb.get_directory() for tb in tb_group["tb"]]
        # Plot and collect figures for each of the tb objects
        # TODO: implement bins for step grouping
        group_figs = [tb.plot_vdfs(steps="last", save=False, bins=bins)
                      for tb in tb_group["tb"].values]
        xy_figs = [f[0][0] for f in group_figs]
        xz_figs = [f[0][1] for f in group_figs]
        steps = [f[1][0] for f in group_figs]

        # Set each of the lims the these extrema
        for fxy, fxz, step, directory in zip(xy_figs, xz_figs, steps, directories):
            figs.extend([fxy, fxz])
            if False:
                sname = Path(directory).parts[-1]
                if not os.path.isdir(save):
                    os.makedirs(save)
                saveas = pjoin(save, sname+f"_xy{step}.pdf")
                print(f"Saving VDF joint plot in {saveas}", flush=True)
                fxy.savefig(saveas)
                plt.close(fxy)
                saveas = pjoin(save, sname+f"_xz{step}.pdf")
                print(f"Saving VDF joint plot in {saveas}", flush=True)
                fxz.savefig(saveas)
                plt.close(fxz)
    return figs

def plot_bilog(x, y, ylog=False, xlabel=None, xlabel_left=None, ylabel=None, **plotting):
    """Make a log plot with two sides, positive x values are plotting on the
    right and negative xvalues are plotting on the left with a flipped axis.

    Args:
        x (Array-Like): The x coordinate data.
        y (Array-Like): The y coordinate data.
        ylog (bool): If set to ``True``, plot the y data on a log scale.
        xlabel (str): label both x axes.
        xlabel_left (str): Use a custom xlabel for the negative x-coordinate
            axes. Otherwise, use ``xlabel`` with ``"-"`` prepended.
        ylabel(str): label both y axes.
        **plotting: Keyword arguments that will be passed to
            :meth:`matplotlib.axes.Axes.plot`
    """

    # Separate data into positive ande negative
    x_pos = x[x>0]
    y_right = y[x>0]
    x_neg = -x[x<0]
    y_left = y[x<0]

    # Find limits based on data ranges
    xmin = np.min(np.abs(x))
    xmax = np.max(np.abs(x))*1.05
    yrange = np.max(y) - np.min(y)
    ymin = np.min(y) - 1.05*yrange
    ymax = np.max(y) + 1.05*yrange

    # Build figure
    fig, axs = plt.subplots(1, 2)
    fig.subplots_adjust(wspace=0)
    # axs[0].spines["right"].set(visible=False)
    axs[1].spines["left"].set(visible=False)
    # Remove overlapping labels
    axs[0].get_yticklabels()[-1].set(visible=False)
    axs[1].get_yticklabels()[0].set(visible=False)
    axs[1].yaxis.set_label_position("right")
    axs[1].yaxis.tick_right()
    # [lbl.set(visible=False) for lbl in axs[0].get_yticklabels()]
    # [lbl.set(visible=False) for lbl in axs[1].get_yticklabels()]
    # axs[1].set_yticks([])


    axs[0].set_xscale("log")
    axs[1].set_xscale("log")
    axs[0].set_xlim((xmin, xmax))
    axs[1].set_xlim((xmin, xmax))
    axs[0].set_ylim((ymin, ymax))
    axs[1].set_ylim((ymin, ymax))

    axs[0].grid(which='major', color='gray', linestyle='-', alpha=0.5)
    axs[0].grid(which='minor', color='gray', linestyle='-', alpha=0.3)
    axs[1].grid(which='major', color='gray', linestyle='-', alpha=0.5)
    axs[1].grid(which='minor', color='gray', linestyle='-', alpha=0.3)

    # Plot Data
    axs[0].plot(x_neg, y_left, **plotting)
    axs[1].plot(x_pos, y_right, **plotting)
    axs[0].invert_xaxis()

    if xlabel:
        axs[1].set_xlabel(xlabel)
        if xlabel_left:
            axs[0].set_xlabel(xlabel_left)
        else:
            axs[0].set_xlabel("-"+xlabel)
    if ylabel:
        axs[0].set_ylabel(ylabel)
        axs[1].set_ylabel(ylabel)

    return fig


def make_pair_plot(ax, vdfs, bins=150):
    """Plot the marginal and joint correlation distributions of the three
    velocity fields."""
    pass

def get_hashable(h):
    """Return a hashable string that describes the plotting style
    of a handle."""
    props = [h.get_ls(), h.get_lw(), h.get_c(),
        h.get_ms(), h.get_marker()]
    s = "".join(list(map(str, props)))
    return s

def figure_legend(axs, order=None):
    """Take series across all axes and combine them into one legend."""
    # Store unique handle representations.
    hs = set()
    # Map back to the actual line objects.
    lhmap = OrderedDict()
    for ax in axs:
        hand, lab = ax.get_legend_handles_labels()
        lm = {get_hashable(h)+"@"+l: h for l, h in zip(lab, hand)}
        lhmap.update(lm)
    labels = [l.split("@")[1] for l in lhmap.keys()]
    handles = list(lhmap.values())
    if order == None:
        return handles, labels
    return [handles[i] for i in [labels.index(j) for j in order]], order


