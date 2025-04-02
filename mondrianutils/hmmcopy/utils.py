import os
import shutil

import csverve.api as csverve
import mondrianutils.helpers as helpers
import pandas as pd
import pysam
import yaml
from multiprocessing import Pool
from mondrianutils import __version__
from mondrianutils.dtypes.hmmcopy import dtypes as hmmcopy_dtypes
from mondrianutils.hmmcopy.plot_heatmap import PlotPcolor
from mondrianutils.hmmcopy.readcounter import ReadCounter, merge_cell_read_counts
import cell_cycle_classifier.api as cell_cycle


def get_cells_from_bam(bam_filename, tag_name):
    with pysam.AlignmentFile(bam_filename, 'rb') as bam:
        if isinstance(bam.header, dict):
            comments = [v.replace(tag_name+':', '') for v in bam.header['CO']]
            return comments

        cells = []
        for line in str(bam.header).split('\n'):
            if not line.startswith("@CO"):
                continue
            line = line.strip().split()
            cb = line[1]
            cell = cb.split(':')[1]
            cells.append(cell)

        return cells


def readcounter(
        infile,
        outdir,
        chromosomes,
        exclude_list=None,
        ncores=16,
        mapping_quality_threshold=20,
        window_size=500000,
        tag_name='CB',
        discover_cells=False,
        ignore_missing_tags=False,
        tabular=False
):
    readcounter = ReadCounter(
        infile, window_size, mapping_quality_threshold,
        excluded=exclude_list, tag_name=tag_name,
        ignore_missing_tags=ignore_missing_tags)

    if ncores == 1:
        counts = []
        for chromosome in chromosomes:
            counts.append(readcounter.get_data(chromosome))

    else:
        with Pool(processes=ncores) as pool:
            counts = pool.map(readcounter.get_data, chromosomes)

    counts = merge_cell_read_counts(counts, chromosomes)

    cells = None
    if not discover_cells:
        cells = get_cells_from_bam(infile, tag_name)

    if tabular:
        counts.write_dataframe(outdir, cells=cells)
    else:
        counts.write_wig_files(outdir, cells=cells)


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


def create_segs_tar(segs_files, metrics, pass_tar, fail_tar, tempdir, segs_samples=None):
    helpers.makedirs(tempdir)
    if segs_samples is not None:
        segs_samples = {k: v for k, v in zip(segs_files, segs_samples)}
    else:
        segs_samples = {}

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
        cell_id = segs_samples.get(filepath, None)
        if cell_id is None:
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


def cell_cycle_classifier(reads, alignment_metrics, hmmcopy_metrics, tempdir, output):
    helpers.makedirs(tempdir)
    cell_cycle_preds = os.path.join(tempdir, 'cellcycle.csv.gz')

    reads_df = csverve.read_csv(reads)
    reads_df['chr'] = reads_df['chr'].astype(str)

    align_metrics_df = csverve.read_csv(alignment_metrics)
    hmm_metrics_df = csverve.read_csv(hmmcopy_metrics)

    predictions = cell_cycle.train_classify(
        reads_df, hmm_metrics_df, align_metrics_df, figures_prefix=None
    )
    csverve.write_dataframe_to_csv_and_yaml(predictions, cell_cycle_preds, dtypes=hmmcopy_dtypes()['metrics'])
    csverve.merge_csv([hmmcopy_metrics, cell_cycle_preds], output, on='cell_id', how='outer')
