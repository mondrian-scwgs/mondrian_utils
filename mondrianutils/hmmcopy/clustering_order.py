import csverve.api as csverve
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as hc
from mondrianutils.dtypes import hmmcopy_metrics
import scipy.spatial as sp


def sort_bins(bins, chromosomes):
    bins = bins.drop_duplicates()

    if not chromosomes:
        chromosomes = list(map(str, range(1, 23))) + ['X', 'Y']

    bins["chr"] = pd.Categorical(bins["chr"], chromosomes)

    bins = bins.sort_values(['start', ])

    bins = [tuple(v) for v in bins.values.tolist()]

    return bins


def get_hierarchical_clustering_order(
        reads_filename, chromosomes=None):
    data = []
    chunksize = 10 ** 5
    for chunk in csverve.read_csv(
            reads_filename, chunksize=chunksize):
        chunk["bin"] = list(zip(chunk.chr, chunk.start, chunk.end))

        # for some reason pivot doesnt like an Int64 state col
        chunk['state'] = chunk['state'].astype('float')

        chunk = chunk.pivot(index='cell_id', columns='bin', values='state')

        data.append(chunk)

    # merge chunks, sum cells that get split across chunks
    table = pd.concat(data)
    table = table.groupby(table.index).sum()

    bins = pd.DataFrame(
        table.columns.values.tolist(),
        columns=[
            'chr',
            'start',
            'end'])

    bins['chr'] = bins['chr'].astype(str)

    bins = sort_bins(bins, chromosomes)

    table = table.sort_values(bins, axis=0)

    data_mat = np.array(table.values)

    data_mat[np.isnan(data_mat)] = -1

    row_linkage = hc.linkage(sp.distance.pdist(data_mat, 'cityblock'),
                             method='ward')

    order = hc.leaves_list(row_linkage)

    samps = table.index
    order = [samps[i] for i in order]
    order = {v: i for i, v in enumerate(order)}

    return order


def add_clustering_order(
        reads, metrics, output, chromosomes=None):
    """
    adds sample information to metrics in place
    """

    order = get_hierarchical_clustering_order(
        reads, chromosomes=chromosomes
    )

    annotation_df = []
    for cellid, rank in order.items():
        annotation_df.append({'cell_id': cellid, 'clustering_order': rank})
    annotation_df = pd.DataFrame(annotation_df)

    csverve.annotate_csv(metrics, annotation_df, output, hmmcopy_metrics.dtypes)
