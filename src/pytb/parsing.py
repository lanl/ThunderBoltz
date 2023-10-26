"""Tools for parsing various formats of cross section and
transport parameter data."""

import io
import itertools
import os
from os import listdir as ls
from os.path import join as pjoin
import re
import time
import warnings

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

def read_LXCat_cs(efile):
    """Convert LXCat cross section text files into
    python pandas DataFrames."""
    # Process types
    ptypes = ["ELASTIC", "EFFECTIVE", "ROTATIONAL",
              "EXCITATION", "IONIZATION", "ATTACHMENT"]
    # Some conditions indicating the beginning end of sections in the file
    # End of file
    EOF = lambda l: l == ""
    # Data base block start and exit
    db_start = lambda l: "DATABASE" in l
    db_end = lambda l: "x"*10 in l
    # Target block start and exit
    tg_start = lambda l: "*"*58 + " " in l
    tg_end = lambda l: "*"*80 in l
    # Process block start
    p_start = lambda l: l.rstrip() in ptypes
    def goto(l, *conditions):
        """Skip lines in string stream until a condition is met."""
        while not any(c(l) for c in conditions):
            l = f.readline()
        return l

    def read_L3(d, l):
        """Interpret the special third line of the LXCat format."""
        s = l.split()
        if len(s) > 1 and ": " not in l:
            # This must be an excitation with a population ratio
            d["threshold"] = float(s[0])
            d["weight_ratio"] = float(s[1])
            l = f.readline()
        elif d["process_type"] == "ATTACHMENT":
            # 3rd line is omitted
            pass
        elif len(s) > 1 and ": " in l:
            # 3rd line is missing for some reason
            pass
        elif d["process_type"] in ["ELASTIC", "EFFECTIVE"]:
            # 3rd line is the ratio between the electron and target mass
            d["target_mass_ratio"] = float(l.rstrip())
            l = f.readline()
        else:
            # 3rd line is the threshold value
            d["threshold"] = float(l.rstrip())
            l = f.readline()
        return l

    def scan_process(l):
        """Scan one block representing a single cross section table.
        Assume `l` begins at a `process_type` line in the LXCat file.
        """
        d = {"process_type": l.rstrip()}
        l = f.readline()
        d["process_string"] = l.rstrip()
        l = f.readline()
        l = read_L3(d, l)
        while "----" not in l:
            meta_data = l.split(": ")
            key = meta_data[0]
            msg = ": ".join(meta_data[1:])
            d[key.lower()] = msg.strip()
            l = f.readline()
        if "columns" not in d:
            d["columns"] = ["Energy (eV)", "Cross section (m2)"]
        else:
            d["columns"] = d["columns"].split(" | ")
        l = f.readline()
        es = []
        css = []
        while "-"*5 not in l:
            e, cs = map(float, l.split())
            es.append(e)
            css.append(cs)
            l = f.readline()
        # Add metadata
        df = pd.DataFrame(np.array([es, css]).T, columns=d["columns"])
        for k, v in d.items():
            if k != "columns":
                df[k] = v
        return df, l

    df = pd.DataFrame()
    with open(efile, "r") as f:
        l = f.readline()
        while not EOF(l):
            l = goto(l, db_start)
            db = l.split()[1]
            l = goto(l, tg_start)
            while not db_end(l):
                target = l.split()[1]
                l = goto(l, p_start)
                while not (tg_end(l) or db_end(l)):
                    pdf, l = scan_process(l)
                    # Update target and db info, add to df
                    pdf["target"] = target
                    pdf["db"] = db
                    df = pd.concat((df, pdf), ignore_index=True)
                    l = goto(l, p_start, tg_end, db_end)
                l = goto(l, tg_start, db_end)
            l = goto(l, db_start, EOF)
    return df

def get_cs_defaults(infile):
    infile = "thunderboltz/input.in"
    cols = ["fid", "r1", "r2", "rtype", "B", "p1", "p1"]
    with open(infile, "r") as f:
        df = pd.DataFrame()
        for line in f.readlines():
            if "CS" not in line or "CC" in line:
                continue
            l = line.split()[1:]
            fid = l[0].split("/")[-1].split(".")[0]
            dat = [fid, int(l[1]), int(l[2]), l[3], float(l[4]),
                   int(l[5]), int(l[6])]
            df = df.append({k: v for k, v in zip(cols, dat)},
                ignore_index=True)

def read_LXCat_swarmp(efile):
    """Read an LXCat swarm parameter text file and return a DataFrame."""
    df = pd.DataFrame()
    with open(efile, "r") as f:
        # For every db
        l = f.readline()
        while l:
            # Find first db
            maxout = 0
            while "DATABASE:" not in l and maxout < 1000:
                maxout += 1
                l = f.readline()
            if maxout > 999:
                break
            db = l.split()[1].upper()
            # For every author
            while "xxxxxx" not in l:
                while "SPECIES" not in l:
                    l = f.readline()
                conditions = T = None
                sp = " ".join(l.split()[1:])
                while "COMMENT" not in l:
                    l = f.readline()
                auth = " ".join(l.split()[1:])
                while "COLUMNS" not in l:
                    l = f.readline()
                x_col = l.split(":")[1].split("|")[0].strip()
                y_col = l.split(":")[1].split("|")[1].strip()
                l = f.readline()
                l = f.readline()
                x, y = [], []
                while "---" not in l:
                    if "UPDATED" in l:
                        continue
                    x_,  y_ = list(map(float, l.split()))
                    x.append(x_)
                    y.append(y_)
                    l = f.readline()
                n = len(x)
                dbs = [db]*n
                auths = [auth]*n
                conds = [conditions]*n
                Ts = [T]*n
                df2= pd.DataFrame({"db": dbs, "author": auths,
                    x_col: x, y_col: y, "conditions": conds,
                    "T": Ts})
                df = df.append(df2,ignore_index=True, sort=True)
                while ("xxxxxx" not in l and "SPECIES"not in l):
                    l = f.readline()
    return df


def extract_vdf(filename, sample_cap=500000, v=False):
    """Read a single particle dump file and return a DataFrame.

    Args:
        filename (str): The directory to search for particle dump files.
        sample_cap (int): Limit the number of samples read from the dump
            file for very large files. Default is 500000. If bool(sample_cap)
            evaluates to ``False``, then no cap will be imposed.
        v (int): Verbosity -- 0: silent, 1: print file paths before reading.

    Returns:
        :class:`pandas.DataFrame`: The velocity distributions function sample
            data.
    """
    if v: print(f"Reading {filename}", flush=True)

    # Read data
    with open(filename) as f:
        if sample_cap:
            # Read the first `sample_cap` lines if available
            txt = "".join(itertools.islice(f, int(sample_cap)))
        else:
            txt = "".join(f.readlines())

    try:
        vdf_arr = np.genfromtxt(
            io.StringIO(txt), delimiter=",", autostrip=True)
    except ValueError as E:
        warnings.warn(f"Skipping {efile}: {str(E)}", RuntimeWarning)
        return

    # Build up dataframe
    cols = ["vx", "vy", "vz"]
    vdf_df = pd.DataFrame(vdf_arr, columns=cols)
    return vdf_df

def extract_vdfs(calc_dir, steps="all", sample_cap=500000,
        particle_type=0, v=0):
    """
    Pull all electron velocity distribution data from a calculation directory
    into a DataFrame.

    Args:
        calc_dir (str): The directory to search for particle dump files.
        steps (list[int] or str): What time steps to extract from the dump files.
            If a list is passed, the files with corresponding step integers
            contained both within the list and the calculation directory
            will be extracted.
        sample_cap (int): Limit the number of samples read from the dump
            file for very large files. Default is 500000. If bool(sample_cap)
            evaluates to ``False``, then no cap will be imposed.
        v (int): Verbosity -- 0: silent, 1: print file paths before reading.

    Returns:
        :class:`pandas.DataFrame`: The velocity distributions function sample
            data.
    """
    # Initialize empty data frame
    df = pd.DataFrame({k: pd.Series([], dtype=dtype)
        for k, dtype in zip(["ptype", "step", "vx", "vy", "vz"],
                            [np.int64]*2 + [np.float64]*3)})

    if steps == []: return df # Return nothing, don't read anything

    if isinstance(particle_type, int):
        particle_type = [particle_type]

    for efile in ls(calc_dir):
        if "Dump" in efile: # Only read velocity dump files

            # Extract particle type and time step out of each filename.
            toks = efile.split("_")
            ptype = int(toks[1])
            step = int(toks[-1].split(".")[0])

            # Skip unrequested data
            if steps != "all" and step not in steps:
                continue
            if particle_type != "all" and ptype not in particle_type:
                continue

            vdf_df = extract_vdf(
                pjoin(calc_dir, efile), sample_cap=sample_cap, v=v)

            if vdf_df is None: continue

            vdf_df["step"] = step
            vdf_df["ptype"] = ptype
            df = pd.concat((df, vdf_df), ignore_index=True)

    df = df.sort_values("step").reset_index(drop=True)
    return df
