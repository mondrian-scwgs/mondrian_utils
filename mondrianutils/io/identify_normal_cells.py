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


def remove_blacklist_bins(bins, reference_name):
    if reference_name.lower() == 'grch37':
        return [v for v in bins if not (v[0] == '9' and v[1] >= 38500001 and v[1] <= 69500001)]
    else:
        raise NotImplementedError('Only GRCh37 is supported at this time')


def get_mean_copy_by_chromosome(reads_data, chromosome):
    copies = reads_data[[v for v in reads_data if v[0] in (chromosome, 'chr' + chromosome)]]
    copies = np.mean(copies.to_numpy()).round()

    return copies


def load_reads_data(reads_data, cells):
    df = helpers.load_and_pivot_reads_data(reads_data, 'state')
    df = df[df.index.isin(cells)]
    return df


def annotate_metrics(input_metrics, output_metrics, normal_cells):
    normal_cells = set(normal_cells)
    df = csverve.read_csv(input_metrics)

    df['is_normal'] = df['cell_id'].apply(lambda x: 'Normal' if x in normal_cells else 'Tumor')

    organisms = [v for v in df.columns.values if v.startswith('fastqscreen_')]
    organisms = sorted(set([v.split('_')[1] for v in organisms]))
    organisms = [v for v in organisms if v not in ['nohit', 'total']]

    csverve.write_dataframe_to_csv_and_yaml(
        df, output_metrics,
        hmmcopy_dtypes(fastqscreen_genomes=organisms)['metrics']
    )


def identify_normal_cells(
        hmmcopy_reads_path: str,
        metrics_data_path: str,
        output_yaml: str,
        output_csv: str,
        reference_name: str,
        min_reads: int = 500000,
        min_quality: float = 0.85,
        allowed_aneuploidy_score: float = 0,
        relative_aneuploidy_threshold: float = 0.05,
        ploidy_threshold: float = 2.5
):
    metrics = load_metrics(metrics_data_path, min_reads, min_quality)

    cn_data = load_reads_data(hmmcopy_reads_path, metrics['cell_id'])

    # drop Y-chromosome immediately
    my_columns = np.where([a[0] not in ['Y']
                           for a in cn_data.columns])[0]
    cn_data = cn_data.iloc[:, my_columns]

    observations = pd.DataFrame(index=cn_data.index)
    observations['ploidy'] = cn_data.apply(lambda x: np.nanmean(x), axis=1)
    observations['rel_aneuploidy'] = get_relative_aneuploidy(cn_data)

    normal_reads = get_supernormal_reads_data(
        cn_data, observations,
        relative_aneuploidy_threshold=relative_aneuploidy_threshold,
        ploidy_threshold=ploidy_threshold
    )

    xcopies = get_mean_copy_by_chromosome(normal_reads, 'X')

    if not (xcopies == 1 or xcopies == 2):
        warnings.warn(f"Found abnormal sex chromosome copies: chrX={xcopies}")

    non_blacklist_bins = remove_blacklist_bins(cn_data.columns, reference_name)

    observed = cn_data[non_blacklist_bins]
    observed = observed.to_numpy()

    expected = pd.DataFrame(columns=non_blacklist_bins)
    expected.loc[0] = 2
    expected[[v for v in expected if v[0] == 'X']] = xcopies
    expected = expected.to_numpy()

    aneu = cdist(observed, expected, metric='hamming')
    aneu_norm = aneu / np.nanmax(aneu)

    observations['aneuploidy_score'] = aneu_norm

    normal_cells = observations.query(f'aneuploidy_score <= {allowed_aneuploidy_score}').index

    with open(output_yaml, 'wt') as writer:
        yaml.dump(list(normal_cells), writer, default_flow_style=False)

    annotate_metrics(metrics_data_path, output_csv, normal_cells)
