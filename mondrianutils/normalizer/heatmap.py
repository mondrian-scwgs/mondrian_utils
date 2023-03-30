from itertools import groupby

import csverve.api as csverve
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.colors import rgb2hex
from mondrianutils import helpers
from mpl_toolkits.axes_grid1 import make_axes_locatable


def generate_colormap_heatmap(df):
    """generating a custom heatmap 2:gray 0: blue 2+: reds
    :param maxval highest value in the data
    :returns listedcolormap
    """
    vmax = df.max().max()
    vmin = df.min().min()

    maxval = 20

    color_reference = {0: '#3498DB', 1: '#85C1E9', 2: '#D3D3D3'}

    low_max = 3 + int((maxval - 3) / 2) + 1
    hi_max = maxval
    low_states = np.arange(3, low_max)
    hi_states = np.arange(low_max, hi_max)

    low_cmap = matplotlib.cm.get_cmap('OrRd', low_max + 1)
    hi_cmap = matplotlib.cm.get_cmap('RdPu', hi_max + 1)

    for cn_level in low_states:
        rgb = low_cmap(int(cn_level))[:3]
        color_reference[cn_level] = rgb2hex(rgb)

    for cn_level in hi_states:
        rgb = hi_cmap(int(cn_level))[:3]
        color_reference[cn_level] = rgb2hex(rgb)

    color_reference[maxval] = "#000000"

    colors = [color_reference[cnlevel]
              for cnlevel in np.arange(vmin, vmax + 1)]

    cmap = ListedColormap(colors)

    return cmap


def get_chr_idxs(bins):
    """
    :param bins: sorted bins used for the plot
    :return chr_idxs: list with the index where chromosome changes
    returns the index where the chromosome changes
    used for marking chr boundaries on the plot
    """
    # chr 1 starts at beginning
    chr_idxs = [0]

    chrom = bins[0][0]
    for i, bin_v in enumerate(bins):
        if bin_v[0] != chrom:
            chr_idxs.append(i)
            chrom = bin_v[0]

    return chr_idxs


def load_metrics(metrics_file):
    metrics = csverve.read_csv(metrics_file)

    metrics = metrics[~ metrics['aneuploidy_score'].isna()]

    metrics = metrics.set_index('cell_id')
    metrics = metrics.sort_values(by=['aneuploidy_score'])

    return metrics


def hmmcopy_heatmap(reads, metrics, row_column=None):
    cmap = generate_colormap_heatmap(reads)

    if row_column is not None:
        lut = dict(zip(metrics[row_column].unique(), "rbg"))
        row_colors = metrics[row_column].map(lut)
    else:
        row_colors = None

    cluster = sns.clustermap(
        reads,
        col_cluster=False,
        row_cluster=False,
        cmap=cmap,
        row_colors=row_colors,
        yticklabels=False
    )

    chr_idxs = get_chr_idxs(reads.columns.values)
    cluster.ax_heatmap.set_xticks(chr_idxs)
    cluster.ax_heatmap.set_xticklabels([str(v) for v in range(1, 23)] + ['X', 'Y'])

    for val in chr_idxs:
        cluster.ax_heatmap.plot(
            [val, val], [len(reads.columns), 0], ':', linewidth=0.5, color='black', )

    if row_column is not None:
        num_normal_cells = len(metrics[metrics['is_normal'] == 'Normal'])
        total_cells = len(metrics)
        num_tumor_cells = total_cells - num_normal_cells

        borders = np.cumsum([0] + [sum(1 for i in g) for k, g in groupby(row_colors)])
        labels = [f'Normal ({num_normal_cells}/{total_cells})', f'Tumor ({num_tumor_cells}/{total_cells})']

        for b0, b1, label in zip(borders[:-1], borders[1:], labels):
            cluster.ax_row_colors.text(
                -0.06, (b0 + b1) / 2,
                label,
                color='black',
                ha='right',
                va='center',
                rotation=90,
                transform=cluster.ax_row_colors.get_yaxis_transform()
            )

    return cluster


def aneuploidy_heatmap(reads, metrics, output, aneuploidy_score=0.005):
    reads = helpers.load_and_pivot_reads_data(reads, 'state')

    metrics = load_metrics(metrics)

    cell_order = list(metrics.index)

    reads = reads[reads.index.isin(metrics.index)]
    reads = reads.reindex(cell_order)

    cluster = hmmcopy_heatmap(reads, metrics, row_column='is_normal')

    divider = make_axes_locatable(cluster.ax_heatmap)
    ax = divider.append_axes("right", size="50%", pad=0.1)

    line_data = np.array(metrics.aneuploidy_score)
    ax.plot(line_data, range(len(line_data)))
    ax.plot([aneuploidy_score] * len(line_data), range(len(line_data)), linestyle='dashed', alpha=0.25)

    sns.despine(ax=ax, left=True, right=False, trim=True)

    ticks = np.arange(1, len(metrics), 100)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f'{v}/{len(metrics)}' for v in ticks])
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.set_ylim(0.5, len(reads.index) - .5)
    ax.invert_yaxis()

    ax.set_xscale('symlog', linthresh=0.0001)
    ax.set_xlim(0, 1)
    ax.set_xticks(ax.get_xticks(), [float(v) for v in ax.get_xticks()], rotation='vertical')

    sns.set_style('whitegrid')
    ax.xaxis.grid(True, alpha=0.25)
    ax.yaxis.grid(False)

    plt.tight_layout()
    plt.savefig(output)
