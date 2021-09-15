covfig_params = dict(
    figsize=(34, 4),
    colormap="coolwarm_r",
    color_chroms=True,
    ylim=(-1.5, 2.5),
    label_size=15,
)

snpfig_params = dict(
    figsize=(34, 4),
    colormap="coolwarm_r",
    color_chroms=True,
    ylim=(0, 1),
    label_size=15,
)

cnvfig_params = dict(
    figsize=(34, 6),
    colormap="coolwarm_r",
    color_chroms=True,
    label_size=13,
    ylims=dict(cov=(-1.4, 2.5), snp=(0, 1)),
)


vaf1 = dict(
    title="BAF",
    plot_type="scatter",  # ['line', 'scatter']
    data="VAF1",
    plot_args=dict(s=3, color="black", cmap="viridis", alpha=1),
)

vaf2 = dict(
    title="BAF",
    plot_type="scatter",  # ['line', 'scatter']
    data="VAF2",
    plot_args=dict(s=3, color="black", cmap="viridis", alpha=1),
)

log1 = dict(
    title="log2ratio",
    plot_type="scatter",  # ['line', 'scatter']
    data="log2ratio1",
    plot_args=dict(linewidth=0.3, color="black", s=0.5, alpha=0.7),
)

log2 = dict(
    title="log2ratio",
    plot_type="scatter",  # ['line', 'scatter']
    data="log2ratio2",
    plot_args=dict(linewidth=0.3, color="black", s=0.5, alpha=0.7),
)

L2Rmean1 = dict(
    title="rollinglog2ratio",
    plot_type="line",  # ['line', 'ascatter']
    data="log2ratio1_mean",
    plot_args=dict(linewidth=1, color="yellow", alpha=0.7),
)

L2Rmean2 = dict(
    title="rollinglog2ratio",
    plot_type="line",  # ['line', 'ascatter']
    data="log2ratio2_mean",
    plot_args=dict(linewidth=1, color="yellow", alpha=0.7),
)
