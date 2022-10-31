import csverve.api as csverve
import numpy as np
import pandas as pd
import yaml
from mondrianutils import helpers
from scipy.spatial.distance import cdist
from scipy.stats import mode


def get_relative_aneuploidy(cn_data):
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


def get_supernormal_reads_data(reads, observations):
    supernormal = observations.query('rel_aneuploidy < 0.05').query('ploidy < 2.5')
    return reads[reads.index.isin(supernormal.index)]


def remove_blacklist_bins(bins, reference_name):

    if reference_name.lower() == 'grch37':
        return [v for v in bins if not (v[0] == '9' and v[1] >= 38500001 and v[1] <= 69500001)]
    else:
        raise NotImplementedError('Only GRCh37 is supported at this time')

def get_mean_copy_by_chromosome(reads_data, chromosome):
    copies = reads_data[[v for v in reads_data if v[0] in (chromosome, 'chr'+chromosome)]]
    copies = np.mean(copies.to_numpy()).round()

    return copies


def load_reads_data(reads_data, cells):
    df = helpers.load_and_pivot_reads_data(reads_data, 'state')
    df = df[df.index.isin(cells)]
    return df


def identify_normal_cells(
        hmmcopy_reads_path: str,
        metrics_data_path: str,
        output_yaml: str,
        reference_name: str,
        min_reads: int = 500000,
        min_quality: float = 0.85,
        allowed_aneuploidy_score: float = 0,
):
    metrics = load_metrics(metrics_data_path, min_reads, min_quality)

    cn_data = load_reads_data(hmmcopy_reads_path, metrics['cell_id'])

    observations = pd.DataFrame(index=cn_data.index)
    observations['ploidy'] = cn_data.apply(lambda x: np.nanmean(x), axis=1)
    observations['rel_aneuploidy'] = get_relative_aneuploidy(cn_data)

    normal_reads = get_supernormal_reads_data(cn_data, observations)

    xcopies = get_mean_copy_by_chromosome(normal_reads, 'X')
    ycopies = get_mean_copy_by_chromosome(normal_reads, 'Y')

    assert (xcopies == 2 and ycopies == 0) or \
           (xcopies == 1 and ycopies == 1), \
        f"Found abnormal sex chromosome copies: chrX={xcopies}, chrY={ycopies}"

    non_blacklist_bins = remove_blacklist_bins(cn_data.columns, reference_name)

    observed = cn_data[non_blacklist_bins]
    observed = observed.to_numpy()

    expected = pd.DataFrame(columns=non_blacklist_bins)
    expected.loc[0] = 2
    expected[[v for v in expected if v[0] == 'X']] = xcopies
    expected[[v for v in expected if v[0] == 'Y']] = ycopies
    expected = expected.to_numpy()

    aneu = cdist(observed, expected, metric='hamming')
    aneu_norm = aneu / np.nanmax(aneu)

    observations['aneuploidy_score'] = aneu_norm

    normal_cells = observations.query(f'aneuploidy_score <= {allowed_aneuploidy_score}').index

    with open(output_yaml, 'wt') as writer:
        yaml.dump(list(normal_cells), writer, default_flow_style=False)
