import warnings

import csverve.api as csverve
import numpy as np
import pandas as pd
import yaml
from mondrianutils import helpers
from mondrianutils.dtypes.hmmcopy import dtypes as hmmcopy_dtypes
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


def load_reads_data(reads_data, cells):
    df = helpers.load_and_pivot_reads_data(reads_data, 'state')
    df = df[df.index.isin(cells)]
    return df


def get_overlapping_bin(value, bins):
    if np.isnan(value):
        return 'nan'

    contains = bins.contains(value)

    idx = np.where(contains == True)

    assert len(idx) == 1

    return str(bins[idx][0])


def annotate_metrics(input_metrics, output_metrics, normal_cells, aneuploidy_scores):
    normal_cells = set(normal_cells)
    df = csverve.read_csv(input_metrics)

    df['is_normal'] = df['cell_id'].apply(lambda x: 'Normal' if x in normal_cells else 'Tumor')

    df['aneuploidy_score'] = df['cell_id'].apply(lambda x: aneuploidy_scores.get(x, float('nan')))

    organisms = [v for v in df.columns.values if v.startswith('fastqscreen_')]
    organisms = sorted(set([v.split('_')[1] for v in organisms]))
    organisms = [v for v in organisms if v not in ['nohit', 'total']]

    csverve.write_dataframe_to_csv_and_yaml(
        df, output_metrics,
        hmmcopy_dtypes(fastqscreen_genomes=organisms)['metrics']
    )


def get_coverage(metrics, normal_cells):
    normal_metrics = metrics[metrics['cell_id'].isin(normal_cells)]
    readcount = sum(normal_metrics['total_reads'])

    coverage = (readcount / 3e9) * 150

    return readcount, coverage

def identify_normal_cells(
        hmmcopy_reads_path: str,
        metrics_data_path: str,
        output_yaml: str,
        output_csv: str,
        blacklist_file: str = None,
        min_reads: int = 500000,
        min_quality: float = 0.85,
        # TODO: change allowed_aneuploidy_score to "normal_allowed_aneuploidy_score" for clarity
        allowed_aneuploidy_score: float = 0,
        # TODO: change relative_aneuploidy_threshold to "supernormal_aneuploidy_threshold" for clarity
        relative_aneuploidy_threshold: float = 0.005,
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
    output_yaml : str
        Path to the output YAML file.
    output_csv : str
        Path to the output CSV file.
    blacklist_file : str, optional
        Path to the blacklist file.
    min_reads : int, optional
        Minimum number of reads to consider a cell.
    min_quality : float, optional
        Minimum quality to consider a cell.
    supernormal_aneuploidy_threshold : float, optional
        Maximum allowed aneuploidy score for a cell to be considered for normal X copy estimation. 
    normal_aneuploidy_score_threshold : float, optional
        Relative aneuploidy threshold for a score 
    """
    metrics = load_metrics(metrics_data_path, min_reads, min_quality)

    cn_data = load_reads_data(hmmcopy_reads_path, metrics['cell_id'])

    # drop Y-chromosome immediately -- unreliable for CN inference
    my_columns = np.where([a[0] not in ['Y']
                           for a in cn_data.columns])[0]
    cn_data = cn_data.iloc[:, my_columns]

    # compute ploidy and relative aneuploidy
    observations = pd.DataFrame(index=cn_data.index)
    observations['ploidy'] = cn_data.apply(lambda x: np.nanmean(x), axis=1)
    observations['rel_aneuploidy'] = get_relative_aneuploidy(cn_data)

    # use obviously normal ("supernormal") cells to estimate normal X copy number
    normal_reads = get_supernormal_reads_data(
        cn_data, observations,
        relative_aneuploidy_threshold=relative_aneuploidy_threshold,
        ploidy_threshold=ploidy_threshold
    )
    xcopies = get_mean_copy_by_chromosome(normal_reads, 'X')
    if not (xcopies == 1 or xcopies == 2):
        warnings.warn(f"Found abnormal sex chromosome copies: chrX={xcopies}")

    # remove blacklist bins
    if blacklist_file is not None:
        non_blacklist_bins = remove_blacklist_bins(cn_data.columns, blacklist_file)
    else:
        non_blacklist_bins = cn_data.columns
    observed = cn_data[non_blacklist_bins]
    observed = observed.to_numpy()

    # determine expected normal copy-number profile
    expected = pd.DataFrame(columns=non_blacklist_bins)
    expected.loc[0] = 2
    expected[[v for v in expected if v[0] == 'X']] = xcopies
    expected = expected.to_numpy()

    # compute aneuploidy score
    aneu = cdist(observed, expected, metric='hamming')
    observations['aneuploidy_score'] = aneu
    assert not observations.aneuploidy_score.isna().any()

    # filter cells
    normal_cells = observations.query(f'aneuploidy_score <= {allowed_aneuploidy_score}').index
    aneuploidy_scores = observations['aneuploidy_score'].to_dict()

    num_reads_normal, coverage_normal = get_coverage(metrics, normal_cells)

    with open(output_yaml, 'wt') as writer:
        yamldata = {
            'cells': list(normal_cells),
            'params': {
                'min_reads': min_reads,
                'min_quality': min_quality,
                'allowed_aneuploidy_score': allowed_aneuploidy_score,
                'relative_aneuploidy_threshold': relative_aneuploidy_threshold,
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

    annotate_metrics(metrics_data_path, output_csv, normal_cells, aneuploidy_scores)
