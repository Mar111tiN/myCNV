import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# use seaborn plotting defaults
import seaborn as sns

from plot_helpers import (
    get_chrom_df,
    extract_pos,
    make_color_chroms,
    add_chrom_labels,
    set_ticks,
    sort_df,
)

sns.set()


def plot_cov(
    df,
    plots,
    chroms="all",
    color_chroms=True,
    colormap="coolwarm_r",
    region="",
    figsize=(20, 4),
    ylim=(-1, 1),
    label_size=12,
):

    # ### DATA MANGLING ##########
    # get cols for rearranging
    org_cols = list(df.columns)

    # sort the df
    df = sort_df(df)
    # reduce the df to the selected chromosomes
    if region:
        chrom, start, end = extract_pos(region)
        df = df.query("Chr == @chrom and @start <= Pos <= @end")
    elif chroms != "all":
        df = df.query("Chr in @chroms")

    # get the chrom_df for collapsing the
    chrom_df = get_chrom_df(df)

    df = df.merge(chrom_df.loc[:, "dif"], on="Chr")
    df["PlotPos"] = df["FullExonPos"] - df["dif"]

    # rearrange the df as return value
    new_cols = org_cols[:4] + ["PlotPos"] + org_cols[4:]
    df = df.loc[:, new_cols]

    # ####### PLOTTING #######
    # plot the figure
    fig, ax = plt.subplots(figsize=figsize)

    # set the x-axis limits
    _ = ax.set_xlim(0, df["PlotPos"].max())

    # plot the graphs #######
    for plot in plots:
        if plot["plot_type"] == "line":
            plot = ax.plot(df["PlotPos"], df[plot["data"]], **plot["plot_args"])
        elif plot["plot_type"] == "scatter":
            plot = ax.scatter(df["PlotPos"], df[plot["data"]], **plot["plot_args"])

    _ = ax.set_ylim(ylim)
    # add the color chroms
    _ = make_color_chroms(
        ax, chrom_df, color_chroms, ylimits=ax.get_ylim(), colormap=colormap
    )

    # ####### LABELS ###################

    # quick fix for one y-label
    _ = ax.set_ylabel(
        " / ".join([plot["title"] for plot in plots]), fontsize=1.25 * label_size
    )

    # ####### CHROM LABELS #############
    add_chrom_labels(ax, chrom_df, ax.get_ylim())

    # set chrom borders
    for m in chrom_df["min"][1:]:
        _ = ax.axvline(x=m, c="k", lw=0.5, alpha=0.5, ls="-")
    # ###### X-AXIS ####################
    # set major ticks and grid for chrom
    ax = set_ticks(ax, df, chrom_df, label_size=label_size)

    # return fig and ax for further plotting and return edited dataframe
    return fig, ax, df, chrom_df


def plot_snp(
    df,
    plots=[],
    chroms="all",
    color_chroms=True,
    colormap="coolwarm_r",
    region="",
    label_size=12,
    figsize=(20, 4),
    ylim=(-1, 1),
):

    # ### DATA MANGELING ##########
    # get cols for rearranging
    org_cols = list(df.columns)

    # sort the df
    df = sort_df(df)
    # reduce the df to the selected chromosomes
    if region:
        chrom, start, end = extract_pos(region)
        df = df.query("Chr == @chrom and @start <= Pos <= @end")
    elif chroms != "all":
        df = df.query("Chr in @chroms")

    # get the chrom_df for collapsing
    chrom_df = get_chrom_df(df)
    df = df.merge(chrom_df.loc[:, "dif"], on="Chr")
    df["PlotPos"] = df["FullExonPos"] - df["dif"]

    # rearrange the df as return value
    new_cols = org_cols[:4] + ["PlotPos"] + org_cols[4:]
    df = df.loc[:, new_cols]

    # ########################
    # ####### PLOTTING #######
    # plot the figure
    fig, ax = plt.subplots(figsize=figsize)

    # set the x-axis limits
    _ = ax.set_xlim(0, df["PlotPos"].max())

    # ####### plot the SNP graphs #######
    for plot in plots:
        if plot["plot_type"] == "line":
            plot = ax.plot(df["PlotPos"], df[plot["data"]], **plot["plot_args"])
        elif plot["plot_type"] == "scatter":
            # highjack plot_args with
            pa = plot["plot_args"]
            if "c" in pa:
                pa["c"] = df[pa["c"]]
            if "s" in pa:
                if isinstance(pa["s"], str):
                    pa["s"] = df[pa["s"]] * 20 + 1
            plot = ax.scatter(df["PlotPos"], df[plot["data"]], **pa)

    _ = ax.set_ylim(ylim)
    # add the color chroms
    _ = make_color_chroms(
        ax, chrom_df, color_chroms, ylimits=ax.get_ylim(), colormap=colormap
    )

    # ####### LABELS ###################

    # quick fix for one y-label
    _ = ax.set_ylabel(
        " / ".join([plot["title"] for plot in plots]), fontsize=1.25 * label_size
    )

    # ####### CHROM LABELS #############
    add_chrom_labels(ax, chrom_df, ax.get_ylim())

    # ###### X-AXIS ####################
    # set major ticks and grid for chrom

    ax = set_ticks(ax, df, chrom_df, label_size=label_size)
    # set helper lines
    #  _ = ax.axhline(y=1, c="k", lw=2, ls="-")
    _ = ax.axhline(y=0.5, c="k", lw=1.5, alpha=0.5, ls="--")

    # set chrom borders
    for m in chrom_df["min"][1:]:
        _ = ax.axvline(x=m, c="k", lw=0.5, alpha=0.5, ls="-")
    # return fig and ax for further plotting and return edited dataframe
    return fig, ax, df, chrom_df


def plot_CNV(
    df,
    snp_plots=[],
    cov_plots=[],
    chroms="all",
    color_chroms=True,
    colormap="coolwarm_r",
    region="",
    label_size=12,
    figsize=(20, 4),
    ylims=dict(cov=(-1, 2.5), snp=(0, 1)),
):

    # ### DATA MANGELING ##########
    # get cols for rearranging
    org_cols = list(df.columns)
    # sort the df
    df = sort_df(df)
    # reduce the df to the selected chromosomes
    if region:
        chrom, start, end = extract_pos(region)
        df = df.query("Chr == @chrom and @start <= Pos <= @end")
    elif chroms != "all":
        df = df.query("Chr in @chroms")

    # get the chrom_df for collapsing
    chrom_df = get_chrom_df(df)
    df = df.merge(chrom_df.loc[:, "dif"], on="Chr")
    df["PlotPos"] = df["FullExonPos"] - df["dif"]

    # rearrange the df as return value
    new_cols = org_cols[:4] + ["PlotPos"] + org_cols[4:]
    df = df.loc[:, new_cols]

    # ########################
    # ####### PLOTTING #######
    # plot the figure
    fig, axes = plt.subplots(2, figsize=figsize, gridspec_kw={"height_ratios": [1, 2]})

    for ax in axes:
        # set the x-axis limits
        _ = ax.set_xlim(0, df["PlotPos"].max())

    # ######## plot COVERAGE
    for plot in cov_plots:

        df[plot["data"]] = df[plot["data"]]
        if plot["plot_type"] == "line":
            _ = axes[0].plot(df["PlotPos"], df[plot["data"]], **plot["plot_args"])

        elif plot["plot_type"] == "scatter":
            # highjack plot_args
            pa = plot["plot_args"]
            if "c" in pa:
                pa["c"] = df[pa["c"]]
            if "s" in pa:
                if isinstance(pa["s"], str):
                    pa["s"] = df[pa["s"]] * 20 + 1
            _ = axes[0].scatter(df["PlotPos"], df[plot["data"]], **pa)

    _ = axes[0].set_ylim(ylims["cov"])

    # ####### plot the SNP graphs #######
    for plot in snp_plots:
        if plot["plot_type"] == "line":
            plot = axes[1].plot(df["PlotPos"], df[plot["data"]], **plot["plot_args"])
        elif plot["plot_type"] == "scatter":
            # highjack plot_args with
            pa = plot["plot_args"]
            if "c" in pa:
                pa["c"] = df[pa["c"]]
            if "s" in pa:
                if isinstance(pa["s"], str):
                    pa["s"] = df[pa["s"]] * 20 + 1
            plot = axes[1].scatter(df["PlotPos"], df[plot["data"]], **pa)

    _ = axes[1].set_ylim(ylims["snp"])

    # add the color chroms
    for ax in axes:
        _ = make_color_chroms(
            ax, chrom_df, color_chroms, ylimits=ax.get_ylim(), colormap=colormap
        )

    # ####### LABELS ###################

    # quick fix for one y-label
    _ = axes[0].set_ylabel("COV [log2r]", fontsize=1.25 * label_size)
    _ = axes[1].set_ylabel("BAF", fontsize=1.25 * label_size)

    # ###### X-AXIS ####################
    # chrom lables
    add_chrom_labels(axes[1], chrom_df, ax.get_ylim())
    # set major ticks and grid for chrom
    axes[1] = set_ticks(axes[1], df, chrom_df, label_size=label_size)
    # remove ticks for coverage plot
    axes[0].xaxis.set_tick_params(which="both", labelbottom=False)
    # set helper lines
    # cov_plot
    for line_pos in [-1, 0, 1]:
        _ = axes[0].axhline(y=line_pos, c="k", lw=1.5, alpha=0.5, ls="--")
    #  VAF plot
    _ = axes[1].axhline(y=0.5, c="k", lw=1.5, alpha=0.5, ls="--")

    # set chrom borders
    for m in chrom_df["min"][1:]:
        for ax in axes:
            _ = ax.axvline(x=m, c="k", lw=0.5, alpha=0.5, ls="-")
    # return fig and ax for further plotting and return edited dataframe
    return fig, axes, df, chrom_df


# ############# CLUSTER FIGURES ##########################################
def plot_2d(df, xcol, ycol, df2=pd.DataFrame(), figsize=(5, 5)):
    fig, ax = plt.subplots(figsize=figsize)
    _ = ax.scatter(df[xcol], df[ycol], s=0.1)
    if len(df2.index):
        _ = ax.scatter(df2[xcol], df2[ycol], s=1, alpha=0.5, color="red")
    _ = ax.set_xlabel(xcol, fontsize=10)
    _ = ax.set_ylabel(ycol, fontsize=10)

    def get_lims(col):
        if "log" in col:
            return (-1.5, 3)
        if "abs" in col:
            return (0, 1)
        if col == "deltaVAFvar":
            return (0, 0.2)
        if col == "deltaVAFstd":
            return (0, 1)
        if col == "VAF":
            return (0, 1)
        else:
            return (-1, 1)

    _ = ax.set_xlim(get_lims(xcol))
    _ = ax.set_ylim(get_lims(ycol))


def plot_3d(df, xcol, ycol, zcol, df2=pd.DataFrame(), figsize=(10, 10)):
    fig = plt.figure()
    ax = plt.axes(projection="3d")

    _ = ax.scatter3D(df[xcol], df[ycol], df[zcol], color="green", alpha=0.2, s=0.1)
    if len(df2.index):
        _ = ax.scatter3D(df2[xcol], df2[ycol], df2[zcol], s=5, color="red")
    # labels
    _ = ax.set_xlabel(xcol, fontsize=10)
    _ = ax.set_ylabel(ycol, fontsize=10)
    _ = ax.set_zlabel(zcol, fontsize=10)

    def get_lims(col):
        if "log" in col:
            return (-1.5, 3)
        if "abs" in col:
            return (0, 1)
        if col == "deltaVAFvar":
            return (0, 0.2)
        if col == "deltaVAFstd":
            return (0, 1)
        else:
            return (-1, 1)

    _ = ax.set_xlim(get_lims(xcol))
    _ = ax.set_ylim(get_lims(ycol))
    _ = ax.set_zlim(get_lims(zcol))
    return fig, ax


def make_GC_plot(cov_df, sample="", agg="mean", max_plots=99):
    """
    create GC plot for the coverages
    """
    cov_cols = [col for col in cov_df.columns if col.startswith("Cov")][:max_plots]
    # create the agg dictionary
    cov_agg = {col: agg for col in cov_cols}
    # make the agg
    df = (
        cov_df.loc[cov_df["map50_0"] > 0.5, :]
        .loc[cov_df["map30_0"] > 0.5, :]
        .loc[cov_df["map75_1"] > 0.5, :]
        .groupby("GCratio")
        .agg(cov_agg)
    )
    fig, ax = plt.subplots(figsize=(10, 10))
    for col in cov_cols:
        _ = ax.plot(df.index, df[col], alpha=0.4)
    _ = ax.set_xlabel("GCratio", fontsize=14)
    _ = ax.set_ylabel("Coverage", fontsize=14)
    if sample:
        _ = ax.set_title(f"Sample {sample} | GCratio vs coverage", fontsize=20)
    return fig, ax
