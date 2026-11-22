import re

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from .tools import set_chrom_categorical


def use_seaborn_style():
    """The seaborn default look (darkgrid, notebook context, deep palette) that main set with `sns.set()`, from the
    style sheet matplotlib ships, so seaborn (and its scipy dependency) is not needed."""
    plt.style.use("seaborn-v0_8")


def sort_df(df):
    """
    helper for sorting dfs for chromosomes (natural order)
    """
    return set_chrom_categorical(df).sort_values(["Chr", "FullExonPos"], kind="stable")


def get_chrom_df(df):

    # observed=True: only chromosomes present in the data
    chrom_df = df.groupby("Chr", observed=True)["FullExonPos"].agg(["mean", "min", "max"]).dropna()
    chrom_df.index = pd.Index(chrom_df.index.astype(str), name="Chr")
    cols = list(chrom_df.columns)
    chrom_df["sum"] = chrom_df["max"] - chrom_df["min"]
    chrom_df["cummin"] = chrom_df["sum"].cumsum()
    chrom_df["dif"] = (chrom_df["max"] - chrom_df["cummin"]).astype(int)
    for col in cols:
        chrom_df[col] = (chrom_df[col] - chrom_df["dif"]).astype(int)
    cols.append("dif")
    return chrom_df.loc[:, cols]


def make_color_chroms(
    ax, chrom_df, color_chroms, ylimits=(-10, 10), colormap="coolwarm_r"
):

    # set the cmap from provided argument
    cmap = matplotlib.colormaps[colormap].resampled(23)
    # build one rect per chromosome present (main built 23 rects but coloured only the present ones)
    rects = []
    # set the height and ymin beyond the ylimits so borders are not seen
    ymin = ylimits[0] * 1.1
    height = (ylimits[1] - ymin) * 1.1

    for _, row in chrom_df.iterrows():
        rects.append(Rectangle((row["min"], ymin), width=row["max"] - row["min"], height=height))
    if color_chroms:
        rect_kwargs = dict(alpha=0.4, ec="none")
    else:
        rect_kwargs = dict(alpha=1, fc="none", ec="darkgray", lw=1, ls="-")
    # set the rectangle collection with colormap
    rect_collection = PatchCollection(rects, cmap=cmap, **rect_kwargs)
    # colour index from the chromosome number (X = 23, Y = 24, others 25)
    names = [re.sub(r"^chr", "", str(c)) for c in chrom_df.index]
    rect_collection.set_array(
        np.array([int(n) if n.isdigit() else {"X": 23, "Y": 24}.get(n, 25) for n in names])
    )

    # setting clim allows fixing the color for individual chroms
    # https://stackoverflow.com/questions/6028675/setting-color-range-in-matplotlib-patchcollection
    rect_collection.set_clim([0, 23])
    return ax.add_collection(rect_collection)


def add_chrom_labels(ax, chrom_df, ylimits=(-10, 10)):

    # YOFFSET is the upper-relative y-position
    YOFFSET = 1
    # get the min_chrom_fraction from minimum chrom_size / whole stretch
    min_chrom_frac = (chrom_df["max"] - chrom_df["min"]).min() / chrom_df["max"].max()
    chrom_size = min(25, max(15, 200 * min_chrom_frac))
    style = dict(size=chrom_size, color="#2f3832")
    # set the height and ymin beyond the ylimits so borders are not seen
    ypos = ylimits[0] + YOFFSET * (ylimits[1] - ylimits[0])
    for chrom, row in chrom_df.iterrows():
        if len(chrom_df.index) > 12:
            chrom = chrom.replace("chr", "")
        ax.text(row["mean"], ypos, chrom, ha="center", **style)


def make_nice(position):
    """
    takes position and returns closest multiple of 1, 2, 5 or 10
    """
    # set nice values
    nice_positions = np.array([1, 2, 2.5, 5, 10])
    # get the 10s
    power10 = np.power(10, np.floor(np.log10(position)))
    # reduce to value between 1 and 10
    num = position / power10
    # find the closest nice position
    base = nice_positions[np.argmin(np.abs(nice_positions / num - 1))]
    return base * power10


def get_tick_pos(tick_dist, chrom_df):
    """
    return from chrom_df the evenly-spread (tick_dist) positions per chrom
    """
    return [
        pos
        for _, row in chrom_df.iterrows()
        for pos in range(row["min"] + tick_dist, row["max"], tick_dist)
    ]


def str_pos(pos, df, precision=1):
    """
    returns string representation of base position
    on genomic coords
    """
    pos = df.iloc[np.argmin(np.abs(df["PlotPos"] - pos))]["Pos"]
    # get the closest base power
    power10 = int(np.round(np.log10(pos) / 3) * 3)
    # get the base fraction
    base = pos / np.power(10, power10)
    if power10 == 9:
        base = base * 1000
        power10 = 6
    if power10 == 6:
        suff = "Mb"
    elif power10 == 3:
        suff = "kb"
    base = re.sub(r"\.0$", "", str(round(base, precision)))
    if power10 == 0:
        suff = "b"
        base = int(base)
    # only print the suff if print_suff
    return f"{base}{suff}"


def get_precision(pos_list):
    """
    get the major tick precision from the range of chrom values
    """
    # get the range of positions
    prange = max(pos_list) - min(pos_list)
    # get the closest 10base
    power10 = int(np.round(np.log10(prange)))
    precision = max(7 - power10, 0)
    return precision


def set_ticks(ax, df, chrom_df, ticks=20, label_size=12):
    """
    for a given tick number, set nicely spread ticks
    """

    ## determine optimale tick distance
    # get the chrom_number
    chrom_count = len(chrom_df.index)
    # get the number of bases
    stretch = chrom_df["max"].iloc[-1]
    # set the number of desired ticks
    major_tick_dist = int(stretch / (ticks + 1))
    minor_tick_dist = int(stretch / ((ticks * 2) + 1))

    # feed tick distance into chrom_df to get chrom-specific coords
    major_pos = get_tick_pos(major_tick_dist, chrom_df)
    minor_pos = [pos - minor_tick_dist for pos in major_pos]

    ax.xaxis.set_major_locator(plt.FixedLocator(major_pos))
    # only print the genomic coords below a certain base total
    if stretch < 2e7:
        precision = get_precision(major_pos)
        major_labels = [
            str_pos(pos, df, precision=precision) for pos in major_pos  ###############
        ]
        ax.xaxis.set_major_formatter(plt.FixedFormatter(major_labels))
        # set the axis labels
        ax.set_xlabel("genomic coords", fontsize=1.25 * label_size)
    else:
        ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_minor_locator(plt.FixedLocator(minor_pos))
    ax.xaxis.grid(which="major", linestyle="-", linewidth=2)
    ax.xaxis.grid(which="minor", linestyle="--", linewidth=1)
    ax.xaxis.set_tick_params(which="major", length=20, labelsize=label_size)
    ax.yaxis.set_tick_params(which="major", length=20, labelsize=label_size)
    # set the tick labels
    for tick in ax.xaxis.get_majorticklabels():
        tick.set_verticalalignment("bottom")
    return ax


def extract_pos(region):
    def convert(pos):
        if pos.endswith("Mb"):
            pos = int(float(pos.replace("Mb", "")) * 1e6)
        elif pos.endswith("kb"):
            pos = int(float(pos.replace("kb", "")) * 1000)
        else:
            pos = int(pos)
        return pos

    split = region.split(":")
    chrom = split[0]

    # if start and are used
    if len(split) > 1 and "-" in split[1]:
        se = split[1].split("-")
        start = convert(se[0])
        end = convert(se[1])
    else:
        start = 0
        end = 1e10
    return chrom, start, end
