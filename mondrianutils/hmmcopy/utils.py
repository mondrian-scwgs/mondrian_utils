import os
import shutil

import csverve.api as csverve
import mondrianutils.helpers as helpers
import mondrianutils.hmmcopy.classify as classify
import pandas as pd
import pysam
import yaml
from mondrianutils import __version__
from mondrianutils.dtypes.hmmcopy import dtypes as hmmcopy_dtypes
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
    reads = csverve.read_csv(reads, chunksize=1e6)

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
