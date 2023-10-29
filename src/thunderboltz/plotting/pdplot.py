"""Plotting tool for intuitively comparing multidimensional data."""
from collections import OrderedDict
import copy
from functools import reduce
import itertools
import json
import os
from os.path import join as pjoin

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
import numpy as np
from numpy import logical_and as land
import pandas as pd

# TEMP
from thunderboltz.plotting import styles

# Default color cycle
DCC = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
        "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#2fb", "#00F", "#F0F",
        "#039", "#00a", "#801"]

class PandaPlot(object):

    def __init__(self, df,
            x=None, y=None,
            err_bars={}, # Map from y key to errbar key
            label_map=lambda x: x, # Pretty column names
            value_map=lambda x: x, # Pretty parameter value names
            series=[],
            series_styles={},
            multi_legend=False,
            series_label_join=", ",
            series_label_interjoin="=",
            rows=None,
            cols=None,
            row_label_style="title", # title | annotation | None
            col_label_style="title", # title | annotation | None
            title_label_join=", ",
            suptitle=None,
            legend_style=None, # None | fig bottom | fig right
            broadcast_nans=True, # Series with missing ids are duplicated
            xscale="linear",
            yscale="linear",
            plot=True,
            gui_configure=False,
            save_spacing_config=None, # configuration file for plot spacing
            load_spacing_config=None, # configuration file for plot spacing
            config_dir=".", # place to store configuration files
            name="Figure", # figure name
            snapx=False,
            snapy=False,
            legend_ordering=None,
            remove_ext=False, # Remove extraneous labels.
            inset_bounds=None, # List of positions for the inset boxes
            **default_style,
            ):
        self.df = df
        self.x = x
        self.y = y
        self.err_bars = err_bars
        self.label_map = label_map
        self.value_map = value_map
        self.series = series
        self.series_styles = series_styles
        self.multi_legend = multi_legend
        self.series_label_join = series_label_join
        self.rows = rows
        self.cols = cols
        self.row_label_style = row_label_style
        self.col_label_style = col_label_style
        self.title_label_join = title_label_join
        self.suptitle = suptitle
        self.legend_style = legend_style
        self.broadcast_nans = broadcast_nans
        self.fig = None # Plot handles
        self.axs = None # Plot handles
        self.id_vars = [] # variables not bound to axis values
        self.xscale = xscale
        self.yscale = yscale
        self.default_style = default_style
        self.gui_configure = gui_configure
        self.save_spacing_config = save_spacing_config
        self.load_spacing_config = load_spacing_config
        self.config_dir = config_dir
        self.name = name
        self.snapx = snapx
        self.snapy = snapy
        self.remove_ext = remove_ext
        self.legend_ordering = legend_ordering
        self.inset_bounds=inset_bounds
        self.inset_axs = [[]] # store possible inset axes

        self.row_order = None
        self.col_order = None

        self.cfile = pjoin(config_dir, name)+".json"

        # err_bars should only include data in y
        if isinstance(self.y, list):
            self.err_bars = {k: v for k, v in self.err_bars.items() if k in self.y}
        else:
            self.err_bars = {k: v for k, v in self.err_bars.items() if k == self.y}

        self.validate_input()
        if plot:
            self.plot()

    def validate_input(self):
        # Generalize series into list
        if not isinstance(self.series, list):
            self.series = [self.series]
        # Store unique series keys in case of multi legend
        if self.multi_legend:
            series_keys = self.series
        else:
            series_keys = copy.deepcopy(self.series)
        if self.multi_legend:
            series_keys = []
            for k, v in self.series:
                series_keys.extend([k] + sum(v.values(), []))
        # Add style variables to series_keys
        series_keys.extend([k for k in self.df.columns if k in styles.presets])
        # Drop duplicates
        series_keys = list(set(series_keys))

        if isinstance(self.y, list):
            self.row_order = self.y
            # disallow additional figure row dim
            assert self.rows == None
            # melt y values into one df column
            id_vars = [self.x]+series_keys+[self.cols]
            id_vars = [idv for idv in id_vars if idv is not None]
            if id_vars:
                # Melt data with error bars
                df = pd.melt(self.df, id_vars=id_vars, value_vars=self.err_bars.keys(),
                    var_name="variable_y", value_name="value")
                df_err = pd.melt(self.df, id_vars=id_vars, value_vars=self.err_bars.values(),
                    var_name="variable_y", value_name="_error")
                df_werr = pd.concat((df, df_err[["_error"]]), axis=1)

                # Melt data without error bars
                nerr_names = [a for a in self.y if a not in self.err_bars.keys()]
                df_nerr = pd.melt(self.df, id_vars=id_vars, value_vars=nerr_names,
                    var_name="variable_y", value_name="value")

                # Combine data
                self.df = pd.concat((df_werr, df_nerr), ignore_index=False)

                # Reset y and rows
                self.y = "value"
                self.rows = "variable_y"
                self.err_bars = "_error"

        if isinstance(self.x, list):
            self.col_order = self.x
            # disallow additional figure column dim
            assert self.cols == None
            # melt x values into one df column
            id_vars = [self.x]+series_keys+[self.cols]
            id_vars = [idv for idv in id_vars if idv is not None]
            if id_vars:
                self.df = pd.melt(self.df, id_vars=id_vars, value_vars=self.y,
                    var_name="variable_x", value_name="value")
                # Reset x and cols
                self.x = "value"
                self.cols = "variable_x"

        # Determine variables not bound to axis values
        if self.rows not in ["variable_y", None]:
            self.id_vars += [self.rows]
        if self.cols not in ["variable_x", None]:
            self.id_vars += [self.cols]

        if self.broadcast_nans:
            # Check if rows y-values
            for idv in self.id_vars:
                # Select data with and without nans.
                missing = self.df[self.df[idv].isna()].copy()
                not_missing = self.df.loc[~self.df[idv].isna()].copy()
                col_names = not_missing[idv].unique()
                self.df = not_missing.copy()
                for col in col_names:
                    missing[idv] = col
                    self.df = pd.concat((self.df, missing), ignore_index=True)

    def plot(self):
        df = self.df
        # Determine number of panels
        nrows = ncols = 1
        if self.rows:
            nrows = len(df[self.rows].unique())
        if self.cols:
            ncols = len(df[self.cols].unique())
        # TODO: autogenerate figsize more accurately
        self.fig, self.axs = plt.subplots(nrows, ncols, figsize=(5*nrows,
            5*ncols), squeeze=False)
        self.inset_axs = [[[] for j in range(ncols)] for i in range(nrows)]

        if self.multi_legend:
            multi_leg_data = self.build_multi_legend(df)
        else:
            # Set a default cycle for unstyled series
            defc = itertools.cycle(DCC)
            for s in df.series.unique():
                if s in self.series_styles:
                    if "c" not in self.series_styles[s]:
                        self.series_styles[s]["c"] = next(defc)
                else:
                    self.series_styles[s] = {"c": next(defc)}
        # Loop though panels and series
        for i, (row_name, rdf) in enumerate(groupby_iter(df, self.rows, order=self.row_order)):
            for j, (col_name, cdf) in enumerate(groupby_iter(rdf, self.cols, order=self.col_order)):
                # Select panel
                ax = self.axs[i, j]
                if self.inset_bounds is not None:
                    self.add_inset(i, j, ax)
                if self.multi_legend:
                    for handle, label, style in zip(*multi_leg_data):
                        sdf = cdf[require(cdf, handle)].copy()
                        self.plot_simple(i, j, ax, sdf, row_name, col_name, label, style)
                else:
                    for s_name, sdf in groupby_iter(cdf, self.series):
                        s_str = self.series_label_join.join(s_name)
                        # Style for this series
                        style = copy.deepcopy(self.default_style)
                        if s_str in self.series_styles:
                            style.update(self.series_styles[s_str])
                        label = dict_f(self.label_map, s_str)
                        self.plot_simple(i, j, ax, sdf, row_name, col_name, label, style)

        if self.suptitle:
            self.fig.suptitle(self.suptitle)

        if self.remove_ext:
            remove_ext_labels(self.axs)

        fs = 18
        if self.multi_legend:
            fs = 11
        if self.legend_style == "figure bottom":
            handles, labels = figure_legend(self.axs.flatten(), order=self.legend_ordering)
            self.figure_legend = self.fig.legend(handles, labels, loc="lower center",
                    ncol=1, prop={"size": fs, "family": "monospace"})
        if self.legend_style == "figure right":
            handles, labels = figure_legend(self.axs.flatten(), order=self.legend_ordering)
            self.figure_legend = self.fig.legend(handles, labels, loc="center right", ncol=1,
                    prop={"size": fs, "family": "monospace"})
            self.fig.subplots_adjust(right=0.6)

        if self.gui_configure:
            # Try to load a config if it exists.
            if os.path.exists(self.cfile) and self.load_spacing_config:
                self.read_spacing_config(self.cfile)

            plt.show()

            self.snap_lims()

            if self.save_spacing_config:
                self.write_spacing_config(self.cfile)

        if self.load_spacing_config:
            self.read_spacing_config(self.cfile)

        return self.fig, self.axs

    def plot_simple(self, i, j, ax, dat, row_name, col_name, label, style):
        """Plot a single set of data with a fixed style"""
        # if dat[self.x].isna().all():
            # return
        # Order data
        dat = dat.sort_values(self.x).reset_index(drop=True)
        # Plot everything
        if dat.isna().all().all():
            line = ax.plot([], [], label=label, **style)
        else:
            line = ax.plot(dat[self.x], dat[self.y], label=label, **style)

        # Optionally add err bars.
        if self.err_bars and not dat[self.err_bars].isna().all():
            ax.errorbar(dat[self.x], dat[self.y], dat[self.err_bars],
                    ls="", ecolor=line[0].get_color(),
                    elinewidth=1, capsize=2)

        # Do the same for insets, without forming legend
        if self.inset_axs[i][j] and not dat.isna().all().all():
            self.inset_axs[i][j].plot(dat[self.x], dat[self.y], **style)
            if self.err_bars and not dat[self.err_bars].isna().all():
                ax.errorbar(dat[self.x], dat[self.y], dat[self.err_bars],
                        ls="", ecolor=line[0].get_color(), elinewidth=1, capsize=2)
        # Label axis
        ax.set_xlabel(dict_f(self.label_map, self.x))
        ax.set_ylabel(dict_f(self.label_map, self.y))
        # Title panels
        titles = []
        if self.row_label_style=="title" and self.rows != "variable_y":
            titles.append(row_name)
        if self.col_label_style=="title" and self.cols != "variable_x":
            titles.append(col_name)
        ax.set_title(self.title_label_join.join(list(map(str,
            titles))), fontsize=16)
        if self.x == "value":
            ax.set_xlabel(dict_f(self.label_map, col_name))
        if self.y == "value":
            ax.set_ylabel(dict_f(self.label_map, row_name))
        # Configure axes scales
        ax.set_xscale(self.xscale)
        ax.set_yscale(self.yscale)
        ax.grid(which='major', color='gray', linestyle='-', alpha=0.5)
        ax.grid(which='minor', color='gray', linestyle='-', alpha=0.3)
        if self.inset_axs[i][j]:
            self.inset_axs[i][j].set_xscale(self.xscale)
            # self.inset_axs[-1].set_yscale(self.yscale)
            self.inset_axs[i][j].grid(which='major', color='gray', linestyle='-', alpha=0.5)
            self.inset_axs[i][j].grid(which='minor', color='gray', linestyle='-', alpha=0.3)
        plt.minorticks_on()
        if self.legend_style is None:
            # Create legend
            ax.legend(fontsize=14, ncol=2)

    def add_inset(self, i, j, ax):
        if self.inset_bounds[i][j] is None:
            self.inset_axs[i][j] = None
        else:
            inset_ax = plt.axes([0, 0, 1, 1])
            inset_ax.set_axes_locator(InsetPosition(ax, self.inset_bounds[i][j]))
            self.inset_axs[i][j] = inset_ax

    def write_spacing_config(self, cfile):
        fs = tuple(self.fig.get_size_inches())
        lims = [[ax.get_xlim(), ax.get_ylim()] for ax in self.axs.flatten()]
        spps = {p: getattr(self.fig.subplotpars, p) for p in
                ["left", "right", "top", "bottom", "wspace", "hspace"]}
        dat = {"fs": fs, "lims": lims, "spps": spps}
        with open(cfile, "w") as f:
            f.write(json.dumps(dat))

    def read_spacing_config(self, cfile):
        with open(cfile, "r") as f:
            d = json.load(f)
        self.fig.set_size_inches(d["fs"])
        self.fig.subplots_adjust(**d["spps"])
        for ax, lims in zip(self.axs.flatten(), d["lims"]):
            ax.set_xlim(lims[0])
            ax.set_ylim(lims[1])

    def snap_lims(self):
        lims = np.array([[list(ax.get_xlim())+list(ax.get_ylim())
                          for ax in axrow] for axrow in self.axs])
        if self.snapx:
            lims[:,:,0] = np.min(lims[:,:,0], axis=0)
            lims[:,:,1] = np.max(lims[:,:,1], axis=0)
        if self.snapy:
            mins = np.min(lims[:,:,2], axis=1).reshape(-1, 1)
            maxs = np.max(lims[:,:,3], axis=1).reshape(-1, 1)
            lims[:,:,2] = np.repeat(mins, 2, axis=1)
            lims[:,:,3] = np.repeat(maxs, 2, axis=1)

        for axrow, limrow in zip(self.axs, lims):
            for ax, lims in zip(axrow, limrow):
                ax.set_xlim(lims[:2])
                ax.set_ylim(lims[2:])

    def save(self, fdir):
        self.fig.savefig(pjoin(fdir, self.name+".pdf"))

    def build_multi_legend(self, df):
        """Create legend data from DataFrame of raw data, legend table format,
        column label maps, value label maps, and style maps."""
        # find preset keys that should not be included in legend
        pset_keys = [k for k in df if k in styles.presets and k not in [k_[0] for
                     k_ in self.series]]
        # Create style cycles, if they aren't cycles already
        for k in self.series_styles:
            if not isinstance(self.series_styles[k][2], itertools.cycle):
                self.series_styles[k][2] = itertools.cycle(self.series_styles[k][2])
        dfs = []
        # Isolate descriptor variables
        for key, col_d in self.series:
            for val, cols in col_d.items():
                dfs.append(df[df[key] == val][[key]+cols+pset_keys].drop_duplicates(
                    ).reset_index(drop=True))
            # Add everything else
            df_sep = df[~df[key].isin(col_d.keys())][[key]+pset_keys].drop_duplicates(
                ).reset_index(drop=True)
            # If there is anything else
            if len(df_sep):
                dfs.append(df_sep)
        # Store variable handles
        hands = []
        # Store legend strings
        legs = []
        # Store the styles that correspond to each legend strings
        stys = []
        # Loop through tables of legend
        for df_ in dfs:
            df_ = df_.sort_values(list(df_.columns)).reset_index(drop=True)
            # Make dataframe with string data types
            sdf = df_.copy()
            for col in sdf:
                sdf[col] = sdf[col].map(lambda v: self.value_map(col)(v)).astype("str")
            sdf.drop(pset_keys, axis=1, inplace=True)
            # Convert dfs to strings
            formatters = {}
            for col in sdf.columns:
                len_max = sdf[col].str.len().max()
                # Ignore special math operators that won't show up on the plot
                len_max = max(len_max, len(str(self.label_map(col).replace("$", "").replace("^", ""))))
                formatters[col] = len_max
            # Add legend header and blank legend header linestyle
            hands += [{}]
            legs += [" ".join(f"{self.label_map(s):^{mx}s}" for s, mx in formatters.items())]
            stys += [dict(ls="", marker="")] # Add empty style
            for (_, row), (_, srow) in zip(df_.iterrows(), sdf.iterrows()):
                # Add style for this line using actual data
                stys += [{}]
                # First check if there are specified styles from `self.series_styles`
                for col, (vs, styk, styvs) in self.series_styles.items():
                    if vs == None:
                        vs = [row[col]]
                        styvs = styles.DCYC[styk]
                    for v, styv in zip(vs, styvs):
                        if col in row and row[col] == v:
                            stys[-1].update({styk: styv})

                # Then check if there are any preset styles for this row
                for pskey in styles.presets:
                    for psval, pssty in styles.presets[pskey].items():
                        if pskey in row.index and row[pskey] == psval:
                            stys[-1].update(
                                {k: v for k, v in next(pssty).items()
                                if k not in stys[-1]}
                            )
                # Drop pset keys from row
                row = row[[k for k in row.index if k not in pset_keys]].copy()
                # Add legend text to legend list using string data.
                legs += [" ".join(f"{s:<{mx}s}" for col, s, mx in
                    zip(srow.index, srow.values, formatters.values()))]
                # Add data from `row` into handles list
                hands += [row.to_dict()]

        return hands, legs, stys

# Util
def get_hashable(h):
    """Return a hashable string that uniquely
    describes the plotting style of a handle."""
    props = [h.get_ls(), h.get_lw(), h.get_c(),
        h.get_ms(), h.get_marker()]
    s = "".join(list(map(str, props)))
    return s

def figure_legend(axs, order=None):
    """Take several axis and merge their legends into one
    figure legend."""
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
    return [handles[i] for i in order], [labels[i] for i in order]

def dict_f(d, x, default=None):
    if isinstance(d, dict):
        if x in d:
            return d[x]
        elif default is not None:
            return default
        return x
    return d(x)

def groupby_iter(df, groups, order=None):
    if groups == None:
        yield "", df
    elif order == None:
        for name, sub_df in df.groupby(groups):
            yield name, sub_df
    else:
        d = {name: sub_df for name, sub_df in df.groupby(groups)}
        for name in order:
            yield name, d[name]

def remove_ext_labels(axs_):
    # Determine shape of the figure.
    ncols = axs_.shape[1]
    axs = axs_.flatten()
    n_last_row = ncols
    nrows = len(axs) // ncols
    if len(axs) % ncols:
        n_last_row = len(axs) % ncols
        nrows = ((len(axs) - n_last_row) // ncols) + 1
    def split_set(n):
        # Return the number(s) half way to n
        if not n % 2:
            return {n // 2, n // 2 + 1}
        return {(n // 2) + 1}
    axl = [axs[n*ncols: (n+1)*ncols] for n in range(nrows - 1)] + [axs[-n_last_row:]]
    for i, row in enumerate(axl):
        for j, ax in enumerate(row):
            if not i+1 in split_set(nrows) or j:
                axl[i][j].set_ylabel("")
            if not j+1 in split_set(n_last_row) or nrows-i-1:
                axl[i][j].set_xlabel("")

def require(df, d, other=[]):
    """Return a mask where all the columns have values passed by the dictionary."""
    if d == {}:
        return pd.DataFrame(None, columns=df.columns)
    crit = [df[k] == v for k, v in d.items()]
    return reduce(land, crit+other)

