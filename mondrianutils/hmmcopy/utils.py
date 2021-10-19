import os

import argparse
import csverve.api as csverve
import mondrianutils.dtypes.hmmcopy_reads as reads_dtypes
import mondrianutils.helpers as helpers
import mondrianutils.hmmcopy.classify as classify
import pandas as pd
from mondrianutils.dtypes import hmmcopy_metrics
from mondrianutils.dtypes import hmmcopy_params
from mondrianutils.dtypes import hmmcopy_reads
from mondrianutils.dtypes import hmmcopy_segs
from mondrianutils.hmmcopy.correct_read_count import CorrectReadCount
from mondrianutils.hmmcopy.plot_hmmcopy import GenHmmPlots
from mondrianutils.hmmcopy.readcounter import ReadCounter


def plot_hmmcopy(
        reads, segments, params, metrics, ref_genome, segs_out,
        bias_out, num_states=12,
        annotation_cols=None, sample_info=None, max_cn=None
):
    if not annotation_cols:
        annotation_cols = ['cell_call', 'experimental_condition', 'sample_type',
                           'mad_neutral_state', 'MSRSI_non_integerness',
                           'total_mapped_reads_hmmcopy']

    with GenHmmPlots(reads, segments, params, metrics, ref_genome, segs_out,
                     bias_out, num_states=num_states,
                     annotation_cols=annotation_cols,
                     sample_info=sample_info, max_cn=max_cn) as plot:
        plot.main()


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

    csverve.rewrite_csv_file('output/0/reads.csv', reads, dtypes=hmmcopy_reads.dtypes)
    csverve.rewrite_csv_file('output/0/params.csv', params, dtypes=hmmcopy_params.dtypes)
    csverve.rewrite_csv_file('output/0/segs.csv', segments, dtypes=hmmcopy_segs.dtypes)
    csverve.rewrite_csv_file('output/0/metrics.csv', metrics, dtypes=hmmcopy_metrics.dtypes)

    helpers.make_tarfile(output_tarball, tempdir)


def add_mappability(reads, annotated_reads):
    reads = csverve.read_csv_and_yaml(reads, chunksize=100)

    alldata = []
    for read_data in reads:
        read_data['is_low_mappability'] = (read_data['map'] <= 0.9)
        alldata.append(read_data)

    alldata = pd.concat(alldata)

    csverve.write_dataframe_to_csv_and_yaml(
        alldata, annotated_reads, reads_dtypes.dtypes, write_header=True
    )


def add_quality(hmmcopy_metrics_csv, alignment_metrics, tempdir, output, training_data):
    helpers.makedirs(tempdir)
    tempout = os.path.join(tempdir, 'added_quality.csv')

    model = classify.train_classifier(training_data)

    feature_names = model.feature_names_

    data = classify.load_data(hmmcopy_metrics_csv, alignment_metrics,
                              feature_names)

    predictions = classify.classify(model, data)

    classify.write_to_output(
        hmmcopy_metrics_csv,
        tempout,
        predictions)

    csverve.rewrite_csv_file(tempout, output, dtypes=hmmcopy_metrics.dtypes)


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    readcounter = subparsers.add_parser('readcounter')
    readcounter.set_defaults(which='readcounter')
    readcounter.add_argument(
        '--infile'
    )
    readcounter.add_argument(
        '--outdir',
    )

    readcounter.add_argument(
        '--chromosomes',
        nargs='*',
        default=[str(v) for v in range(1, 23)] + ['X', 'Y'],
        help='specify target chromosomes'
    )
    readcounter.add_argument(
        '-w', '--window_size',
        type=int,
        default=1000,
        help='specify bin size'
    )
    readcounter.add_argument(
        '-m', '--mapping_quality_threshold',
        type=int,
        default=0,
        help='threshold for the mapping quality, reads ' \
             'with quality lower than threshold will be ignored'
    )

    readcounter.add_argument(
        '--exclude_list',
        default=None,
        help='regions to skip'
    )

    correct_readcount = subparsers.add_parser('correct_readcount')
    correct_readcount.set_defaults(which='correct_readcount')
    correct_readcount.add_argument(
        '--infile'
    )
    correct_readcount.add_argument(
        '--gc_wig_file'
    )
    correct_readcount.add_argument(
        '--map_wig_file'
    )
    correct_readcount.add_argument(
        '--outfile'
    )
    correct_readcount.add_argument(
        '--cell_id'
    )
    correct_readcount.add_argument(
        '--map_cutoff',
        default=0.9,
        type=float,
    )

    plot_hmmcopy = subparsers.add_parser('plot_hmmcopy')
    plot_hmmcopy.set_defaults(which='plot_hmmcopy')
    plot_hmmcopy.add_argument(
        '--reads'
    )
    plot_hmmcopy.add_argument(
        '--segments'
    )
    plot_hmmcopy.add_argument(
        '--params'
    )
    plot_hmmcopy.add_argument(
        '--metrics'
    )
    plot_hmmcopy.add_argument(
        '--reference'
    )
    plot_hmmcopy.add_argument(
        '--segments_output'
    )
    plot_hmmcopy.add_argument(
        '--bias_output'
    )
    plot_hmmcopy.add_argument(
        '--cell_id'
    )

    run_hmmcopy = subparsers.add_parser('run_hmmcopy')
    run_hmmcopy.set_defaults(which='run_hmmcopy')
    run_hmmcopy.add_argument(
        '--corrected_reads'
    )
    run_hmmcopy.add_argument(
        '--tempdir'
    )
    run_hmmcopy.add_argument(
        '--metrics'
    )
    run_hmmcopy.add_argument(
        '--params'
    )
    run_hmmcopy.add_argument(
        '--reads'
    )
    run_hmmcopy.add_argument(
        '--segments'
    )
    run_hmmcopy.add_argument(
        '--output_tarball'
    )

    add_mappability = subparsers.add_parser('add_mappability')
    add_mappability.set_defaults(which='add_mappability')
    add_mappability.add_argument(
        '--infile'
    )
    add_mappability.add_argument(
        '--outfile'
    )

    add_quality = subparsers.add_parser('add_quality')
    add_quality.set_defaults(which='add_quality')
    add_quality.add_argument(
        '--hmmcopy_metrics'
    )
    add_quality.add_argument(
        '--alignment_metrics'
    )
    add_quality.add_argument(
        '--training_data'
    )
    add_quality.add_argument(
        '--output'
    )
    add_quality.add_argument(
        '--tempdir'
    )

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'readcounter':
        with ReadCounter(args['infile'], args['outdir'], args['window_size'],
                         args['chromosomes'], args['mapping_quality_threshold'],
                         excluded=args['exclude_list']) as rcount:
            rcount.main()
    elif args['which'] == 'correct_readcount':
        CorrectReadCount(args["gc_wig_file"],
                         args['map_wig_file'],
                         args['infile'],
                         args['outfile'],
                         args['cell_id'],
                         mappability=args['map_cutoff']).main()
    elif args['which'] == 'plot_hmmcopy':
        plot_hmmcopy(
            args['reads'], args['segments'], args['params'], args['metrics'],
            args['reference'], args['segments_output'], args['bias_output']
        )
    elif args['which'] == 'run_hmmcopy':
        run_hmmcopy(
            args['corrected_reads'], args['tempdir'], args['metrics'], args['params'],
            args['reads'], args['segments'], args['output_tarball']
        )
    elif args['which'] == 'add_mappability':
        add_mappability(args['infile'], args['outfile'])
    elif args['which'] == 'add_quality':
        add_quality(args['hmmcopy_metrics'], args['alignment_metrics'], args['tempdir'], args['output'],
                    args['training_data'])
    else:
        raise Exception()


utils()
