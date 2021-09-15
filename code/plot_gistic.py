import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import seaborn as sns

sns.set()


def add_chrom_blocks(
    ax, chrom_df, colormap="coolwarm_r", ylimits=(-10, 10), yoffset=5, label_size=1
):
    # set the cmap from provided argument
    cmap = plt.cm.get_cmap("coolwarm_r", 23)
    rects = []
    ylimits = ax.get_ylim()
    ymin = ylimits[0] * 1.1
    height = (ylimits[1] - ymin) * 1.1
    for _, row in chrom_df.iterrows():
        rect = Rectangle(
            (row["FullStart"], ymin),
            width=row["FullEnd"] - row["FullStart"],
            height=height,
        )
        rects.append(rect)
    rect_kwargs = dict(alpha=1, fc="none", ec="black", lw=1, ls="dotted")
    rect_collection = PatchCollection(rects, **rect_kwargs)
    # rect_collection.set_array(
    #    chrom_df['Chr'].str.replace("chr", "").str.replace("X", "23").astype(int)
    # )
    _ = ax.add_collection(rect_collection)

    # add the labels
    # get the min_chrom_fraction from minimum chrom_size / whole stretch
    min_chrom_frac = (chrom_df["FullEnd"] - chrom_df["FullStart"]).min() / chrom_df[
        "FullEnd"
    ].max()
    chrom_size = min(20, max(15, label_size * 200 * min_chrom_frac))

    label_style = dict(size=chrom_size, color="#2f3832")
    # set the height and ymin beyond the ylimits so borders are not seen
    ypos = ylimits[0] - yoffset
    for _, row in chrom_df.iterrows():
        mean = (row["FullEnd"] + row["FullStart"]) / 2
        chrom = row["Chr"].replace("chr", "")
        ax.text(mean, ypos, chrom, ha="center", **label_style)

    return ax


def add_cytoband_blocks(
    ax,
    band_df,
    yoffset=5,
    label_size=1,
    band_colors={"p": "lightgray", "q": "darkgray"},
):
    # set the cmap from provided argument
    # cmap = plt.cm.get_cmap("coolwarm_r", 23)

    ylimits = ax.get_ylim()
    ymin = ylimits[0] * 1.1
    height = (ylimits[1] - ymin) * 1.1
    ypos = ylimits[0] + yoffset

    label_style = dict(size=10 * label_size, color="#2f3832")
    for band in ["p", "q"]:
        rects = []
        for _, row in band_df.query("BAND == @band").iterrows():
            rect = Rectangle(
                (row["FullStart"], ymin),
                width=row["FullEnd"] - row["FullStart"],
                height=height,
            )
            rects.append(rect)
            # add the label
            mean = (row["FullEnd"] + row["FullStart"]) / 2
            ax.text(mean, ypos, band, ha="center", **label_style)

        # each rect collection has separate rect_kwargs
        rect_kwargs = dict(alpha=0.4, fc=band_colors[band], ec="none", lw=0.3, ls="-")
        rect_collection = PatchCollection(rects, **rect_kwargs)
        _ = ax.add_collection(rect_collection)

    return ax


def get_tick_range(df, col):
    _max = df[col].max()
    # positive values
    if _max > 0:
        step = 10 if _max >= 20 else 5

        return [s + step for s in range(0, step * math.ceil(_max / step), step)]
    else:
        _min = df[col].min()
        step = 10 if _min <= -20 else 5
        return [s for s in range(step * math.floor(_min / step), 0, step)]


def set_y_ticks(ax, df, label_size=12, margin=2):
    """
    for a given tick number, set nicely spread ticks
    """

    # get the ticks from the gains and losses range
    y_ticks = get_tick_range(df, "losses") + get_tick_range(df, "gains")

    # set the ylim from y_ticks range
    _ = ax.set_ylim(y_ticks[0] - margin, y_ticks[-1] + margin)

    ax.yaxis.set_major_locator(plt.FixedLocator(y_ticks))
    ax.yaxis.set_major_formatter(plt.FixedFormatter(np.abs(y_ticks)))

    ax.yaxis.grid(which="major", linestyle="dotted", linewidth=1)

    ax.yaxis.set_tick_params(which="major", length=20, labelsize=label_size)
    # set the tick labels
    for tick in ax.yaxis.get_majorticklabels():
        tick.set_horizontalalignment("left")
    return ax


def make_CNV_plot(
    df,
    chrom_df,
    band_df,
    figsize=(10, 4),
    y_margin=2,
    label_size=10,
    cnv_fill_colors={},
    chrom_color_scheme="coolwarm_r",
    chrom_label_yoffset=5,
    chrom_label_size=1,
    band_label_yoffset=5,
    band_label_size=1,
    band_colors={"p": "lightgray", "q": "darkgray"},
    y_label_size=10,
):
    """"""

    ######### data ############
    # get the FullPos from chrom_df into df
    df["FullPos"] = df["Start"] + df.merge(chrom_df, on="Chr")["FullStart"] - 1

    # get the FullPos from chrom_df into band_df
    band_df = band_df.merge(chrom_df, on="Chr")
    band_df["FullEnd"] = band_df["End"] + band_df["FullStart"] - 1
    band_df["FullStart"] = band_df["Start"] + band_df["FullStart"]
    band_df = band_df.loc[:, ["Chr", "BAND", "FullStart", "FullEnd"]]

    # create the figure
    fig, ax = plt.subplots(figsize=figsize)

    # remove seaborn grid
    ax.grid(False)
    # set the x_lims
    _ = ax.set_xlim(0, chrom_df["FullEnd"].max())

    _ = set_y_ticks(ax, df, label_size=y_label_size, margin=y_margin)

    # plot the graphs
    plot = ax.stackplot(
        df["FullPos"], df["gains"], alpha=0.9, color=cnv_fill_colors["gains"]
    )
    plot = ax.stackplot(
        df["FullPos"], df["losses"], alpha=0.9, color=cnv_fill_colors["losses"]
    )

    _ = add_chrom_blocks(
        ax,
        chrom_df,
        colormap=chrom_color_scheme,
        ylimits=ax.get_ylim(),
        yoffset=chrom_label_yoffset,
        label_size=chrom_label_size,
    )

    _ = add_cytoband_blocks(
        ax,
        band_df,
        yoffset=band_label_yoffset,
        label_size=band_label_size,
        band_colors=band_colors,
    )

    # remove the x-plot ticks
    plt.tick_params(
        axis="x",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False,
    )

    return fig, ax