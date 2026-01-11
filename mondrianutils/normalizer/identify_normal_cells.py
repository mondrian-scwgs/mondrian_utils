import warnings

import csverve.api as csverve
import numpy as np
import pandas as pd
import yaml
from mondrianutils import helpers
from mondrianutils.normalizer.heatmap import aneuploidy_heatmap
from scipy.spatial.distance import cdist
from scipy.stats import mode


def get_relative_aneuploidy(cn_data):
    # ignore sex chromosomes for this calculation
    my_columns = np.where([a[0] not in ['X', 'Y']
                           for a in cn_data.columns])[0]
    cn_data = cn_data.iloc[:, my_columns]

    ampdel = np.zeros(cn_data.shape, np.int8)
    modes = mode(cn_data, axis=1, keepdims=True)[0]
    ampdel[cn_data > modes] = 1
    ampdel[cn_data < modes] = -1

    rel_aneuploidy = np.sum(ampdel != 0, axis=1) / ampdel.shape[1]

    return rel_aneuploidy


def load_metrics(
        metrics_file,
        min_reads: int = 500000,
        min_quality: float = 0.85,
):
    metrics = csverve.read_csv(metrics_file)

    metrics = metrics.query('~is_control')
    metrics = metrics.query('is_s_phase == False')

    metrics = metrics.query(f'total_reads >= {min_reads}')
    metrics = metrics.query(f'quality >= {min_quality}')

    return metrics


def get_supernormal_reads_data(
        reads, observations,
        relative_aneuploidy_threshold=0.05,
        ploidy_threshold=2.5
):
    supernormal = observations.query(
        f'rel_aneuploidy < {relative_aneuploidy_threshold}'
    ).query(
        f'ploidy < {ploidy_threshold}'
    )
    return reads[reads.index.isin(supernormal.index)]


def remove_blacklist_bins(bins, blacklist_file):

    blacklist = pd.read_csv(blacklist_file, dtype='str')
    blacklist_bins = blacklist['chr'] + '-' + blacklist['start'] + ':' + blacklist['end']
    return [v for v in bins if v not in blacklist_bins]


def get_mean_copy_by_chromosome(reads_data, chromosome):
    copies = reads_data[[v for v in reads_data if v[0] in (chromosome, 'chr' + chromosome)]]
    copies = np.mean(copies.to_numpy()).round()

    return copies


def get_coverage(metrics, normal_cells):
    normal_metrics = metrics[metrics['cell_id'].isin(normal_cells)]
    readcount = sum(normal_metrics['total_reads'])

    coverage = (readcount / 3e9) * 150

    return readcount, coverage


def identify_normal_cells(
        hmmcopy_reads_path: str,
        metrics_data_path: str,
        normal_copy_path: str,
        output_yaml: str,
        output_csv: str,
        output_plot: str,
        min_reads: int = 500000,
        min_quality: float = 0.85,
        aneuploidy_score_threshold: float = 0.005,
        ploidy_threshold: float = 2.5
):
    """
    Identify normal cells based on the provided HMMcopy reads and metrics data.
    
    Parameters
    ----------
    hmmcopy_reads_path : str
        Path to the HMMcopy reads data.
    metrics_data_path : str
        Path to the metrics data.
    normal_copy_path : str
        Path to the normal copy-number profile.
    output_yaml : str
        Path to the output YAML file.
    output_csv : str
        Path to the output CSV file.
    output_plot : str
        Path to the output aneuploidy heatmap file.
    min_reads : int, optional
        Minimum number of reads to consider a cell.
    min_quality : float, optional
        Minimum quality to consider a cell.
    aneuploidy_score_threshold : float, optional
        Relative aneuploidy threshold score for calling normal cells.
    ploidy_threshold : float, optional
        Ploidy threshold for calling normal cells.
    """

    cn_data = pd.read_csv(
        hmmcopy_reads_path,
        dtype={'chr': 'str', 'cell_id': 'str'},
        usecols=[
            'cell_id',
            'chr',
            'start',
            'end',
            'state',
        ])

    # Filter cells based on metrics
    metrics = load_metrics(metrics_data_path, min_reads, min_quality)
    cn_data = cn_data[cn_data['cell_id'].isin(metrics['cell_id'].values)]

    normal_cn = pd.read_csv(
        normal_copy_path,
        dtype={'chr': 'str'},
        usecols=[
            'chr',
            'start',
            'end',
            'state',
        ])
    normal_cn = normal_cn.set_index(['chr', 'start', 'end'])['state']

    # Create a matrix of cn data
    cn_matrix = cn_data.pivot_table(
        index=('chr', 'start', 'end'),
        columns='cell_id',
        values='state',
    )

    # Subset cn_matrix to bins in normal profile
    cn_matrix = cn_matrix.loc[normal_cn.index].copy()

    scores = pd.DataFrame(index=cn_matrix.columns)
    scores['ploidy'] = cn_matrix.mean(axis=0)
    scores['aneuploidy_score'] = cdist(cn_matrix.values.T, normal_cn.values[:, None].T, metric='hamming')[:, 0]
    assert not scores['aneuploidy_score'].isna().any()

    scores['is_normal'] = scores['aneuploidy_score'] <= aneuploidy_score_threshold

    normal_cells = scores.query('is_normal').index.tolist()

    num_reads_normal, coverage_normal = get_coverage(metrics, normal_cells)

    aneuploidy_heatmap(cn_matrix, scores, output_plot, aneuploidy_score_threshold=aneuploidy_score_threshold)

    with open(output_yaml, 'wt') as writer:
        yamldata = {
            'cells': list(normal_cells),
            'params': {
                'min_reads': min_reads,
                'min_quality': min_quality,
                'aneuploidy_score_threshold': aneuploidy_score_threshold,
                'ploidy_threshold': ploidy_threshold
            },
            'normal_metrics': {
                'total_read_count': num_reads_normal,
                'coverage': coverage_normal,
                'total_eligible_cells': len(metrics),
                'normal_cells': len(normal_cells)
            }
        }
        yaml.dump(yamldata, writer, default_flow_style=False)

    scores.to_csv(output_csv, index=True)
