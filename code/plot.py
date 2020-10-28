import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib as mpl
import matplotlib.pyplot as plt
# use seaborn plotting defaults
import seaborn as sns
sns.set()


def sort_df(df):
    '''
    helper for sorting dfs for chromosomes
    '''
    df2 = df.copy()
    # make Chr column categorical for sorting .. and sort
    chrom_list = [f"chr{i+1}" for i in range(22)] + ['chrX']
    df2['Chr'] = pd.Categorical(df2['Chr'], chrom_list)
    return df2.sort_values(['Chr', 'FullExonPos'])


def get_chrom_df(df):

    # dropna is neccessary because grouping on categorical returns all categories
    chrom_df = df.groupby('Chr')['FullExonPos'].agg(
        ['mean', 'min', 'max']).dropna()
    cols = list(chrom_df.columns)
    chrom_df['sum'] = chrom_df['max'] - chrom_df['min']
    chrom_df['cummin'] = chrom_df['sum'].cumsum()
    chrom_df['dif'] = (chrom_df['max'] - chrom_df['cummin']).astype(int)
    for col in cols:
        chrom_df[col] = (chrom_df[col] - chrom_df['dif']).astype(int)
    cols.append('dif')
    return chrom_df.loc[:, cols]


def make_color_chroms(ax, chrom_df, color_chroms, ylimits=(-10, 10), colormap='coolwarm_r'):

    # set the cmap from provided argument
    cmap = plt.cm.get_cmap(colormap, 23)

    # build the rects
    rects = []
    # set the height and ymin beyond the ylimits so borders are not seen
    ymin = ylimits[0] * 1.1
    height = (ylimits[1] - ymin) * 1.1

    for chrom in chrom_df.index:
        row = chrom_df.loc[chrom]
        rect = Rectangle((row['min'], ymin),
                         width=row['max']-row['min'], height=height)
        rects.append(rect)

    if color_chroms:
        rect_kwargs = dict(
            alpha=.6,
            ec='none'
        )
    else:
        rect_kwargs = dict(
            alpha=1,
            fc='none',
            ec='darkgray',
            lw=1,
            ls='-'
        )
    # set the rectangle collection with colormap
    rect_collection = PatchCollection(rects, cmap=cmap, **rect_kwargs)
    # set the index for the color map from chrom integers
    rect_collection.set_array(chrom_df.index.str.replace(
        'chr', '').str.replace('X', '23').astype(int))
    return ax.add_collection(rect_collection)


def add_chrom_labels(ax, chrom_df, ylimits=(-10, 10)):

    # YOFFSET is the upper-relative y-position
    YOFFSET = 0.9
    # get the min_chrom_fraction from minimum chrom_size / whole stretch
    min_chrom_frac = (chrom_df['max'] - chrom_df['min']
                      ).min() / chrom_df['max'].max()
    chrom_size = min(50, max(15, 200*min_chrom_frac))
    style = dict(size=chrom_size, color='#2f3832')
    # set the height and ymin beyond the ylimits so borders are not seen
    ypos = ylimits[0] + YOFFSET * (ylimits[1] - ylimits[0])
    for chrom, row in chrom_df.iterrows():
        if (len(chrom_df.index) > 12):
            chrom = chrom.replace('chr', '')
        ax.text(row['mean'], ypos, chrom, ha='center', **style)


def make_nice(position):
    '''
    takes position and returns closest multiple of 1, 2, 5 or 10
    '''
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
    '''
    return from chrom_df the evenly-spread (tick_dist) positions per chrom
    '''
    return [pos for _, row in chrom_df.iterrows() for pos in range(row['min'] + tick_dist, row['max'], tick_dist)]


def str_pos(pos, df, precision=1):
    '''
    returns string representation of base position
    on genomic coords
    '''
    pos = df.iloc[np.argmin(np.abs(df['PlotPos'] - pos))]['Pos']
    # get the closest base power
    power10 = int(np.round(np.log10(pos) / 3) * 3)
    # get the base fraction
    base = pos / np.power(10, power10)
    if power10 == 9:
        base = base * 1000
        power10 = 6
    if power10 == 6:
        suff = 'Mb'
    elif power10 == 3:
        suff = 'kb'
    base = round(base, precision)
    if power10 == 0:
        suff = 'b'
        base = int(base)
    return f"{base}{suff}"


def set_ticks(ax, df, chrom_df, ticks=20):
    '''
    for a given tick number, set nicely spread ticks
    '''

    # determine optimale tick distance
    # get the chrom_number
    chrom_count = len(chrom_df.index)
    # get the number of bases
    stretch = chrom_df['max'][-1]
    # set the number of desired ticks
    major_tick_dist = int(stretch / (ticks + 1))
    minor_tick_dist = int(stretch / ((ticks * 2) + 1))

    # feed tick distance into chrom_df to get chrom-specific coords
    major_pos = get_tick_pos(major_tick_dist, chrom_df)
    minor_pos = [pos - minor_tick_dist for pos in major_pos]

    ax.xaxis.set_major_locator(plt.FixedLocator(major_pos))
    # only print the genomic coords below a certain base total
    if stretch < 2e7:
        major_labels = [str_pos(pos, df) for pos in major_pos]
        ax.xaxis.set_major_formatter(plt.FixedFormatter(major_labels))
        # set the axis labels
        _ = ax.set_xlabel('genomic coords')
    else:
        ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_minor_locator(plt.FixedLocator(minor_pos))
    ax.xaxis.grid(which='major', linestyle='-', linewidth=2)
    ax.xaxis.grid(which='minor', linestyle='--', linewidth=1)
    ax.xaxis.set_tick_params(which='major', length=20, )

    # set the tick labels
    for tick in ax.xaxis.get_majorticklabels():
        tick.set_verticalalignment("bottom")
    return ax


def extract_pos(region):

    def convert(pos):
        if pos.endswith('Mb'):
            pos = int(pos.replace('Mb', '')) * 1e6
        elif pos.endswith('kb'):
            pos = int(pos.replace('kb', '')) * 1000
        else:
            pos = int(pos)
        return pos

    split = region.split(':')
    chrom = split[0]

    # if start and are used
    if len(split) > 1 and '-' in split[1]:
        se = split[1].split('-')
        start = convert(se[0])
        end = convert(se[1])
    else:
        start = 0
        end = 1e10
    return chrom, start, end


def plot_genomic(df, plots, chroms='all', color_chroms=True, colormap='coolwarm_r', region='', figsize=(20, 4), ylim=(-1, 1), ):

    #### DATA MANGELING ##########
    # get cols for rearranging
    org_cols = list(df.columns)

    # sort the df
    df = sort_df(df)
    # reduce the df to the selected chromosomes
    if region:
        chrom, start, end = extract_pos(region)
        df = df.query('Chr == @chrom and @start <= Pos <= @end')
    elif chroms != 'all':
        df = df.query('Chr in @chroms')

    # get the chrom_df for collapsing the
    chrom_df = get_chrom_df(df)

    df = df.merge(chrom_df.loc[:, 'dif'], on='Chr')
    df['PlotPos'] = df['FullExonPos'] - df['dif']

    # rearrange the df as return value
    new_cols = org_cols[:4] + ['PlotPos'] + org_cols[4:]
    df = df.loc[:, new_cols]

    ######## PLOTTING #######
    # plot the figure
    fig, ax = plt.subplots(figsize=figsize)
    # set the x-axis limits
    _ = ax.set_xlim(0, df['PlotPos'].max())

    # plot the graphs #######
    # !!!!
    # print('Datapoints', len(df.index))
    # suggestions for the best scatter
    #

    for plot in plots:
        if plot['plot_type'] == 'line':
            plot = ax.plot(df['PlotPos'], df[plot['data']],
                           **plot['plot_args'])
        elif plot['plot_type'] == 'scatter':
            plot = ax.scatter(
                df['PlotPos'], df[plot['data']], **plot['plot_args'])

    _ = ax.set_ylim(ylim)
    # add the color chroms
    _ = make_color_chroms(ax, chrom_df, color_chroms,
                          ylimits=ax.get_ylim(), colormap=colormap)

    ######## LABELS ###################
    # quick fix for one y-label
    label = ' / '.join([plot['title'] for plot in plots])
    if len(label) > 80:
        label = label.replace(' / ', ' /\n')
    _ = ax.set_ylabel(label)

    ######## CHROM LABELS #############
    add_chrom_labels(ax, chrom_df, ax.get_ylim())

    ####### X-AXIS ####################
    # set major ticks and grid for chrom
    ax = set_ticks(ax, df, chrom_df)

    # return fig and ax for further plotting and return edited dataframe
    return fig, ax, df, chrom_df
