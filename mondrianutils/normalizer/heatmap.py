import csverve.api as csverve
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.colors import rgb2hex
from matplotlib.patches import Patch
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

    metrics = metrics[~ metrics['relative_aneuploidy'].isna()]

    metrics = metrics.set_index('cell_id')
    metrics = metrics.sort_values(by=['relative_aneuploidy'])

    return metrics


def hmmcopy_heatmap(reads, metrics, row_column=None):
    cmap = generate_colormap_heatmap(reads)

    if row_column is not None:
        lut = dict(zip(metrics[row_column].unique(), "rbg"))
        row_colors = metrics[row_column].map(lut)
    else:
        row_colors = None
        lut = None

    cluster = sns.clustermap(
        reads,
        col_cluster=False,
        row_cluster=False,
        cmap=cmap,
        row_colors=row_colors
    )

    chr_idxs = get_chr_idxs(reads.columns.values)
    cluster.ax_heatmap.set_xticks(chr_idxs)
    cluster.ax_heatmap.set_xticklabels([str(v) for v in range(1, 23)] + ['X', 'Y'])

    for val in chr_idxs:
        cluster.ax_heatmap.plot(
            [val, val], [len(reads.columns), 0], ':', linewidth=0.5, color='black', )

    if row_column is not None:
        handles = [Patch(facecolor=lut[name]) for name in lut]
        plt.legend(
            handles,
            lut,
            title='Species',
            bbox_to_anchor=(1, 1),
            bbox_transform=plt.gcf().transFigure,
            loc='upper right'
        )

    return cluster


def normal_heatmap(reads, metrics, output):
    reads = helpers.load_and_pivot_reads_data(reads, 'state')

    metrics = load_metrics(metrics)

    cell_order = list(metrics.index)

    reads = reads[reads.index.isin(metrics.index)]
    reads = reads.reindex(cell_order)

    hmmcopy_heatmap(reads, metrics, row_column='is_normal')

    plt.savefig(output)


def aneuploidy_heatmap(reads, metrics, output):
    reads = helpers.load_and_pivot_reads_data(reads, 'state')

    metrics = load_metrics(metrics)

    cell_order = list(metrics.index)

    reads = reads[reads.index.isin(metrics.index)]
    reads = reads.reindex(cell_order)

    cluster = hmmcopy_heatmap(reads, metrics)

    divider = make_axes_locatable(cluster.ax_heatmap)
    ax = divider.append_axes("left", size="50%", pad=0.1)

    line_data = np.array(metrics.relative_aneuploidy)
    ax.plot(line_data, range(len(line_data)))
    ax.plot([0.005] * len(line_data), range(len(line_data)), linestyle='dashed', alpha=0.25)

    sns.despine(ax=ax, left=True, trim=True)

    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.set_ylim(-0.5, len(reads.index) - .5)
    ax.invert_yaxis()

    ax.set_xscale('symlog', linthresh=0.0001)
    ax.set_xlim(1, 0)
    ax.set_xticks(ax.get_xticks(), [float(v) for v in ax.get_xticks()], rotation='vertical')

    sns.set_style('whitegrid')
    ax.xaxis.grid(True, alpha=0.25)
    ax.yaxis.grid(False)

    plt.savefig(output)
