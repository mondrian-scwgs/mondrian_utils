import os
import shutil

import argparse
import csverve.api as csverve
import mondrianutils.helpers as helpers
import mondrianutils.hmmcopy.classify as classify
import pandas as pd
import pysam
import yaml
from mondrianutils import __version__
from mondrianutils.dtypes.hmmcopy import dtypes as hmmcopy_dtypes
from mondrianutils.hmmcopy.clustering_order import add_clustering_order
from mondrianutils.hmmcopy.complete_hmmcopy import complete_hmmcopy
from mondrianutils.hmmcopy.generate_qc_html import generate_html_report
from mondrianutils.hmmcopy.plot_heatmap import PlotPcolor
from mondrianutils.hmmcopy.readcounter import ReadCounter


def _readcounter_command(infile, outdir, chromosome, exclude_list=None, mapping_quality_threshold=0, window_size=1000):
    return [
        'hmmcopy_utils', 'readcounter',
        '--infile', infile,
        '--outdir', outdir,
        '-w', window_size,
        '--chromosomes', chromosome,
        '-m', mapping_quality_threshold,
        '--exclude_list', exclude_list,
        '--ncores', 1
    ]


def _merge_wig(infiles, outfile, cell):
    with open(outfile, 'wt') as writer:
        writer.write(f'track type=wiggle_0 name={cell}\n')
        for infile in infiles:
            with open(infile, 'rt') as reader:
                for line in reader:
                    if line.startswith('track type'):
                        continue
                    writer.write(line)


def readcounter(
        infile,
        outdir,
        tempdir,
        chromosomes,
        exclude_list=None,
        ncores=16,
        mapping_quality_threshold=20,
        window_size=500000
):
    if ncores == 1:
        with ReadCounter(infile, outdir, window_size,
                         chromosomes, mapping_quality_threshold,
                         excluded=exclude_list) as rcount:
            rcount.main()
        return

    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)
    if os.path.exists(outdir):
        shutil.rmtree(outdir)

    os.makedirs(tempdir)
    os.makedirs(outdir)

    cells = helpers.get_cells(pysam.AlignmentFile(infile, 'rb'))

    commands = [
        _readcounter_command(infile, os.path.join(tempdir, chromosome), chromosome, exclude_list=exclude_list,
                             mapping_quality_threshold=mapping_quality_threshold, window_size=window_size)
        for chromosome in chromosomes
    ]

    scripts_tempdir = os.path.join(tempdir, 'scripts_split')

    helpers.run_in_gnu_parallel(commands, scripts_tempdir, ncores)

    for cell in cells:
        wig_files = [os.path.join(tempdir, chromosome, f'{cell}.wig') for chromosome in chromosomes]

        _merge_wig(wig_files, os.path.join(outdir, f'{cell}.wig'), cell)


def plot_heatmap(reads, metrics, chromosomes, output, sidebar_column='pick_met', disable_clustering=None):
    plot = PlotPcolor(
        reads, metrics, output,
        column_name='state',
        max_cn=12,
        scale_by_cells=False,
        mappability_threshold=0.9,
        chromosomes=chromosomes,
        sidebar_column=sidebar_column,
        disable_clustering=disable_clustering
    )
    plot.main()


def add_mappability(reads, annotated_reads):
    reads = csverve.read_csv(reads, chunksize=100)

    alldata = []
    for read_data in reads:
        read_data['is_low_mappability'] = (read_data['map'] <= 0.9)
        alldata.append(read_data)

    alldata = pd.concat(alldata)

    csverve.write_dataframe_to_csv_and_yaml(
        alldata, annotated_reads, hmmcopy_dtypes()['reads'], skip_header=False
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

    data = pd.read_csv(tempout)
    organisms = [v for v in data.columns.values if v.startswith('fastqscreen_')]
    organisms = sorted(set([v.split('_')[1] for v in organisms]))
    organisms = [v for v in organisms if v not in ['nohit', 'total']]

    csverve.rewrite_csv_file(tempout, output, dtypes=hmmcopy_dtypes(fastqscreen_genomes=organisms)['metrics'])


def create_segs_tar(segs_files, metrics, pass_tar, fail_tar, tempdir):
    helpers.makedirs(tempdir)

    metrics_data = csverve.read_csv(metrics)
    all_cells = metrics_data.cell_id.tolist()
    metrics_data = metrics_data[metrics_data['quality'] >= 0.75]
    metrics_data = metrics_data[metrics_data['is_contaminated'] == False]

    good_cells = metrics_data.cell_id.tolist()
    bad_cells = [cell for cell in all_cells if cell not in good_cells]

    pass_dir = os.path.join(tempdir, 'segs_pass')
    helpers.makedirs(pass_dir)
    fail_dir = os.path.join(tempdir, 'segs_fail')
    helpers.makedirs(fail_dir)

    for filepath in segs_files:
        cell_id = open(filepath + '.sample', 'rt').readlines()
        assert len(cell_id) == 1
        cell_id = cell_id[0].strip()

        if cell_id in good_cells:
            shutil.copyfile(filepath, os.path.join(pass_dir, '{}_segments.pdf'.format(cell_id)))
        elif cell_id in bad_cells:
            shutil.copyfile(filepath, os.path.join(fail_dir, '{}_segments.pdf'.format(cell_id)))
        else:
            raise Exception('cell_id {} for file {} not found in metrics file {}'.format(cell_id, filepath, metrics))

    helpers.make_tarfile(pass_tar, pass_dir)
    helpers.make_tarfile(fail_tar, fail_dir)


def generate_metadata(
        metrics, params, reads, segments, segments_tar_pass,
        segments_tar_fail, heatmap, metadata_input, metadata_output
):
    with open(metadata_input, 'rt') as reader:
        data = yaml.safe_load(reader)

    out_data = dict()
    out_data['meta'] = dict(
        type='hmmcopy',
        version=__version__,
        lane_ids=data['meta']['lane_ids'],
        sample_ids=data['meta']['sample_ids'],
        library_ids=data['meta']['library_ids'],
        cell_ids=data['meta']['cell_ids'],
    )

    files = {
        os.path.basename(metrics[0]): {
            'result_type': 'hmmcopy_metrics',
            'auxiliary': helpers.get_auxiliary_files(metrics[0])
        },
        os.path.basename(metrics[1]): {
            'result_type': 'hmmcopy_metrics',
            'auxiliary': helpers.get_auxiliary_files(metrics[1])
        },
        os.path.basename(reads[0]): {
            'result_type': 'hmmcopy_reads',
            'auxiliary': helpers.get_auxiliary_files(reads[0])
        },
        os.path.basename(reads[1]): {
            'result_type': 'hmmcopy_reads',
            'auxiliary': helpers.get_auxiliary_files(reads[1])
        },
        os.path.basename(params[0]): {
            'result_type': 'hmmcopy_params',
            'auxiliary': helpers.get_auxiliary_files(params[0])
        },
        os.path.basename(params[1]): {
            'result_type': 'hmmcopy_params',
            'auxiliary': helpers.get_auxiliary_files(params[1])
        },
        os.path.basename(segments[0]): {
            'result_type': 'hmmcopy_segments',
            'auxiliary': helpers.get_auxiliary_files(segments[0])
        },
        os.path.basename(segments[1]): {
            'result_type': 'hmmcopy_segments',
            'auxiliary': helpers.get_auxiliary_files(segments[1])
        },
        os.path.basename(segments_tar_pass): {
            'result_type': 'hmmcopy_segment_plots', 'filtering': 'quality_passed',
            'auxiliary': helpers.get_auxiliary_files(segments_tar_pass)
        },
        os.path.basename(segments_tar_fail): {
            'result_type': 'hmmcopy_segment_plots', 'filtering': 'quality_failed',
            'auxiliary': helpers.get_auxiliary_files(segments_tar_fail)
        },
        os.path.basename(heatmap): {
            'result_type': 'hmmcopy_heatmap_plots',
            'auxiliary': helpers.get_auxiliary_files(heatmap)
        }
    }

    out_data['files'] = files

    with open(metadata_output, 'wt') as writer:
        yaml.dump(out_data, writer, default_flow_style=False)


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    readcounter = subparsers.add_parser('readcounter')
    readcounter.set_defaults(which='readcounter')
    readcounter.add_argument(
        '--infile',
        required=True
    )
    readcounter.add_argument(
        '--outdir',
        required=True
    )
    readcounter.add_argument(
        '--tempdir',
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
        default=500000,
        help='specify bin size'
    )
    readcounter.add_argument(
        '-m', '--mapping_quality_threshold',
        type=int,
        default=20,
        help='threshold for the mapping quality, reads ' \
             'with quality lower than threshold will be ignored'
    )
    readcounter.add_argument(
        '--exclude_list',
        default=None,
        help='regions to skip'
    )
    readcounter.add_argument(
        '--ncores',
        default=12,
        type=int,
        help='regions to skip'
    )


    run_hmmcopy = subparsers.add_parser('hmmcopy')
    run_hmmcopy.set_defaults(which='hmmcopy')
    run_hmmcopy.add_argument('--readcount_wig', required=True)
    run_hmmcopy.add_argument('--gc_wig_file', required=True)
    run_hmmcopy.add_argument('--map_wig_file', required=True)
    run_hmmcopy.add_argument('--metrics', required=True)
    run_hmmcopy.add_argument('--params', required=True)
    run_hmmcopy.add_argument('--reads', required=True)
    run_hmmcopy.add_argument('--segments', required=True)
    run_hmmcopy.add_argument('--output_tarball', required=True)
    run_hmmcopy.add_argument('--reference', required=True)
    run_hmmcopy.add_argument('--segments_output', required=True)
    run_hmmcopy.add_argument('--bias_output', required=True)
    run_hmmcopy.add_argument('--cell_id', required=True)
    run_hmmcopy.add_argument('--tempdir', required=True)
    run_hmmcopy.add_argument(
        '--map_cutoff',
        default=0.9,
        type=float,
    )

    add_mappability = subparsers.add_parser('add_mappability')
    add_mappability.set_defaults(which='add_mappability')
    add_mappability.add_argument(
        '--infile'
    )
    add_mappability.add_argument(
        '--outfile'
    )

    heatmap = subparsers.add_parser('heatmap')
    heatmap.set_defaults(which='heatmap')
    heatmap.add_argument(
        '--reads'
    )
    heatmap.add_argument(
        '--metrics'
    )
    heatmap.add_argument(
        '--chromosomes',
        default=[str(v) for v in range(1, 23)] + ['X', 'Y'],
        nargs='*'
    )
    heatmap.add_argument(
        '--output'
    )
    heatmap.add_argument(
        '--sidebar_column',
        default='pick_met'
    )
    heatmap.add_argument(
        '--disable_clustering',
        default=False,
        action='store_true'
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

    create_segs_tar = subparsers.add_parser('create_segs_tar')
    create_segs_tar.set_defaults(which='create_segs_tar')
    create_segs_tar.add_argument(
        '--segs_pdf',
        nargs='*'
    )
    create_segs_tar.add_argument(
        '--metrics'
    )
    create_segs_tar.add_argument(
        '--pass_output'
    )
    create_segs_tar.add_argument(
        '--fail_output'
    )
    create_segs_tar.add_argument(
        '--tempdir'
    )

    generate_html_report = subparsers.add_parser('generate_html_report')
    generate_html_report.set_defaults(which='generate_html_report')
    generate_html_report.add_argument(
        '--tempdir',
    )
    generate_html_report.add_argument(
        '--html'
    )
    generate_html_report.add_argument(
        '--metrics'
    )
    generate_html_report.add_argument(
        '--gc_metrics'
    )

    add_clustering_order = subparsers.add_parser('add_clustering_order')
    add_clustering_order.set_defaults(which='add_clustering_order')
    add_clustering_order.add_argument(
        '--reads',
    )
    add_clustering_order.add_argument(
        '--metrics'
    )
    add_clustering_order.add_argument(
        '--output'
    )
    add_clustering_order.add_argument(
        '--chromosomes',
        default=[str(v) for v in range(1, 23)] + ['X', 'Y'],
        nargs='*'
    )

    generate_metadata = subparsers.add_parser('generate_metadata')
    generate_metadata.set_defaults(which='generate_metadata')
    generate_metadata.add_argument(
        '--metrics', nargs=2
    )
    generate_metadata.add_argument(
        '--params', nargs=2
    )
    generate_metadata.add_argument(
        '--reads', nargs=2
    )
    generate_metadata.add_argument(
        '--segments', nargs=2
    )
    generate_metadata.add_argument(
        '--segments_tar_pass'
    )
    generate_metadata.add_argument(
        '--segments_tar_fail'
    )
    generate_metadata.add_argument(
        '--heatmap'
    )
    generate_metadata.add_argument(
        '--metadata_input'
    )
    generate_metadata.add_argument(
        '--metadata_output'
    )

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'readcounter':
        readcounter(
            args['infile'],
            args['outdir'],
            args['tempdir'],
            args['chromosomes'],
            exclude_list=args['exclude_list'],
            ncores=args['ncores'],
            mapping_quality_threshold=args['mapping_quality_threshold'],
            window_size=args['window_size']
        )

    elif args['which'] == 'hmmcopy':
        complete_hmmcopy(
            args['readcount_wig'], args["gc_wig_file"], args['map_wig_file'],
            args['metrics'], args['params'],
            args['reads'], args['segments'], args['output_tarball'],
            args['reference'], args['segments_output'], args['bias_output'],
            args['cell_id'], args['tempdir'],
            mappability_cutoff=args['map_cutoff']
        )
    elif args['which'] == 'add_mappability':
        add_mappability(args['infile'], args['outfile'])
    elif args['which'] == 'add_quality':
        add_quality(args['hmmcopy_metrics'], args['alignment_metrics'], args['tempdir'], args['output'],
                    args['training_data'])
    elif args['which'] == 'create_segs_tar':
        create_segs_tar(
            args['segs_pdf'], args['metrics'],
            args['pass_output'], args['fail_output'],
            args['tempdir']
        )
    elif args['which'] == 'heatmap':
        plot_heatmap(
            args['reads'], args['metrics'], args['chromosomes'], args['output'],
            sidebar_column=args['sidebar_column'], disable_clustering=args['disable_clustering']
        )
    elif args['which'] == 'generate_html_report':
        generate_html_report(
            args['tempdir'], args['html'], args['metrics'], args['gc_metrics']
        )
    elif args['which'] == 'add_clustering_order':
        add_clustering_order(args['reads'], args['metrics'], args['output'], chromosomes=args['chromosomes'])

    elif args['which'] == 'generate_metadata':
        generate_metadata(
            args['metrics'], args['params'], args['reads'], args['segments'],
            args['segments_tar_pass'], args['segments_tar_fail'],
            args['heatmap'], args['metadata_input'], args['metadata_output']
        )
    else:
        raise Exception()


utils()
