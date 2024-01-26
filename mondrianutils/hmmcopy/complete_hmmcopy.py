import os

import csverve.api as csverve
import mondrianutils.helpers as helpers
import pandas as pd

from mondrianutils.dtypes.hmmcopy import dtypes as hmmcopy_dtypes
from mondrianutils.hmmcopy.correct_read_count import CorrectReadCount
from mondrianutils.hmmcopy.plot_hmmcopy import plot_hmmcopy

from mondrianutils.hmmcopy import add_mappability
from mondrianutils.hmmcopy import add_quality





def run_hmmcopy(
        corrected_reads,
        tempdir,
        metrics,
        params,
        reads,
        segments,
        output_tarball,
        multipliers=tuple(range(1, 7)),
        strength=1000,
        e=0.999999,
        mu=tuple(range(12)),
        lambda_p=20,
        nu=2.1,
        kappa=(100, 100, 700, 100, 25, 25, 25, 25, 25, 25, 25, 25),
        m=tuple(range(12)),
        eta=50000,
        g=3,
        s=1
):
    df = pd.read_csv(corrected_reads)
    cell_id = list(df['cell_id'].unique())
    assert len(cell_id) == 1
    cell_id = cell_id[0]

    scripts_directory = os.path.realpath(os.path.dirname(__file__))
    run_hmmcopy_rscript = os.path.join(scripts_directory, 'hmmcopy_single_cell.R')
    cmd = [run_hmmcopy_rscript]

    # run hmmcopy
    cmd += ['--corrected_data=' + corrected_reads,
            '--outdir=' + tempdir,
            '--sample_id=' + cell_id]

    cmd.append('--param_str=' + str(strength))
    cmd.append('--param_e=' + str(e))
    cmd.append('--param_mu=' + ','.join(map(str, mu)))
    cmd.append('--param_l=' + str(lambda_p))
    cmd.append('--param_nu=' + str(nu))
    cmd.append('--param_k=' + ','.join(map(str, kappa)))
    cmd.append('--param_m=' + ','.join(map(str, m)))
    cmd.append('--param_eta=' + str(eta))
    cmd.append('--param_g=' + str(g))
    cmd.append('--param_s=' + str(s))
    cmd.append('--param_multiplier=' + ','.join(map(str, multipliers)))

    helpers.run_cmd(cmd)

    csverve.rewrite_csv_file(f'{tempdir}/0/reads.csv', reads, dtypes=hmmcopy_dtypes()['reads'])
    csverve.rewrite_csv_file(f'{tempdir}/0/params.csv', params, dtypes=hmmcopy_dtypes()['params'])
    csverve.rewrite_csv_file(f'{tempdir}/0/segs.csv', segments, dtypes=hmmcopy_dtypes()['segs'])
    csverve.rewrite_csv_file(f'{tempdir}/0/metrics.csv', metrics, dtypes=hmmcopy_dtypes()['metrics'])

    helpers.make_tarfile(output_tarball, tempdir)


def complete_hmmcopy(
        readcount_wig, gc_wig_file, map_wig_file, alignment_metrics,
        metrics, params, reads, segments,
        output_tarball, reference, segments_output, bias_output, cell_id, tempdir,
        quality_classifier_training_data, quality_classifier_model=None,
        mappability_cutoff=0.9
):
    helpers.makedirs(tempdir)

    corrected_reads = os.path.join(tempdir, 'corrected_reads.csv')
    correction = CorrectReadCount(
        gc_wig_file, map_wig_file, readcount_wig, corrected_reads,
        cell_id, mappability=mappability_cutoff
    )
    correction.main()

    hmmcopy_tempdir = os.path.join(tempdir, 'hmmcopy')
    helpers.makedirs(hmmcopy_tempdir)
    raw_reads = os.path.join(hmmcopy_tempdir, 'reads.csv.gz')
    raw_metrics = os.path.join(hmmcopy_tempdir, 'metrics.csv.gz')
    run_hmmcopy(
        corrected_reads, hmmcopy_tempdir, raw_metrics,
        params, raw_reads, segments, output_tarball
    )
    add_mappability(raw_reads, reads)

    plot_hmmcopy(
        reads, segments, params, raw_metrics, reference,
        segments_output, bias_output
    )

    merged_metrics = os.path.join(hmmcopy_tempdir, 'alignment_merged_metrics.csv.gz')
    csverve.merge_csv(
        [raw_metrics, alignment_metrics], merged_metrics, how='inner', on='cell_id', write_header=True, lenient=True
    )

    add_quality(merged_metrics, metrics, quality_classifier_training_data,
                joblib_model=quality_classifier_model)
