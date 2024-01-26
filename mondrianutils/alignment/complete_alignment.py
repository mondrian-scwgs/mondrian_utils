import os
import pathlib
import shutil
import csverve

import mondrianutils.helpers as helpers
import yaml

from mondrianutils.alignment.coverage_metrics import get_coverage_metrics
from mondrianutils.alignment.fastqscreen import merge_fastq_screen_counts
from mondrianutils.alignment.fastqscreen import organism_filter
from mondrianutils.alignment.trim_galore import trim_galore

from mondrianutils.alignment import utils as alignment_utils

from mondrianutils.alignment import mark_duplicates
from mondrianutils.alignment import insert_metrics
from mondrianutils.alignment import wgs_metrics
from mondrianutils.alignment import flagstat
from mondrianutils.alignment import gc_metrics
from mondrianutils.alignment import tss_enrichment
import pysam
import subprocess


def load_metadata(metadata_yaml, lane_id, flowcell_id, cell_id):
    with open(metadata_yaml, 'rt') as reader:
        data = yaml.safe_load(reader)

    lane_data = data['meta']['lanes'][flowcell_id][lane_id]

    center = lane_data['sequencing_centre']

    cell_data = data['meta']['cells']
    library_id = cell_data[cell_id]['library_id']
    sample_id = cell_data[cell_id]['sample_id']

    return center, library_id, sample_id


def bwa_align(
        fastq1, fastq2, reference, reference_name, metadata_yaml,
        output, lane_id, flowcell_id, cell_id, num_threads, tempdir
):
    center, library_id, sample_id = load_metadata(metadata_yaml, lane_id, flowcell_id, cell_id)

    script_path = pathlib.Path(__file__).parent.resolve()
    script_path = os.path.join(script_path, 'bwa_align.sh')

    if reference_name.lower() == 'grch38':
        script_path = os.path.join(script_path, 'bwa_align.sh')

    cmd = [
        script_path, fastq1, fastq2, reference, output,
        sample_id, library_id, cell_id, lane_id,
        flowcell_id, center, num_threads, tempdir
    ]
    helpers.run_cmd(cmd)


def bam_sort(bam_filename, sorted_bam_filename, tempdir, num_threads=1, mem="2G"):
    if not os.path.exists(tempdir):
        helpers.makedirs(tempdir)

    cmd = ['picard', '-Xmx' + mem, '-Xms' + mem]
    if num_threads == 1:
        cmd.append('-XX:ParallelGCThreads=1')
    cmd.extend([
        'SortSam',
        'INPUT=' + bam_filename,
        'OUTPUT=' + sorted_bam_filename,
        'SORT_ORDER=coordinate',
        'VALIDATION_STRINGENCY=LENIENT',
        'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        'QUIET=true'
    ])
    helpers.run_cmd(cmd)


def bam_index(infile):
    helpers.run_cmd([
        'samtools', 'index',
        infile,
        infile + '.bai'
    ])


def add_cell_id_to_header(bamfile, cellid, outputbam, tempdir):
    helpers.makedirs(tempdir)
    new_header = os.path.join(tempdir, 'header.sam')

    subprocess.run(['samtools', 'view', '-H', bamfile, '-o', new_header])
    with open(new_header, 'at') as header:
        header.write('@CO\tCB:{}\n'.format(cellid))

    subprocess.run(
        ['picard', 'ReplaceSamHeader', 'I={}'.format(bamfile),
         'HEADER={}'.format(new_header), 'O={}'.format(outputbam)
         ]
    )


def alignment(
        fastq_pairs, metadata_yaml, reference,
        supplementary_references,
        tempdir, adapter1, adapter2, cell_id, wgs_metrics_mqual, wgs_metrics_bqual, wgs_metrics_count_unpaired,
        bam_output, metrics_output, metrics_gc_output,
        tar_output, num_threads, run_fastqc=False
):
    reference_name, reference_version, reference = reference.split(',')
    supplementary_reference_files = [v.split(',')[2] for v in supplementary_references]
    supplementary_reference_names = [v.split(',')[0] for v in supplementary_references]

    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)

    final_lane_bams = []
    all_detailed_counts = []
    all_summary_counts = []
    for lane in fastq_pairs:
        lane = lane.split(',')
        assert len(lane) == 4, lane
        r1 = lane[2]
        r2 = lane[3]
        lane_id = lane[0]
        flowcell_id = lane[1]

        print("Processing Lane {} Flowcell {}".format(lane_id, flowcell_id))

        lane_tempdir = os.path.join(tempdir, lane_id, flowcell_id)
        assert not os.path.exists(lane_tempdir)

        helpers.makedirs(os.path.join(lane_tempdir, 'fastqscreen'))
        fastqscreen_r1 = os.path.join(lane_tempdir, 'fastqscreen', 'r1.fastq.gz')
        fastqscreen_r2 = os.path.join(lane_tempdir, 'fastqscreen', 'r2.fastq.gz')
        detailed_metrics = os.path.join(lane_tempdir, 'fastqscreen', 'detailed.csv.gz')
        summary_metrics = os.path.join(lane_tempdir, 'fastqscreen', 'summary.csv.gz')
        fastqscreen_temp = os.path.join(lane_tempdir, 'fastqscreen')

        print("Starting FastqScreen")
        organism_filter(
            r1, r2, fastqscreen_r1, fastqscreen_r2,
            detailed_metrics, summary_metrics,
            fastqscreen_temp, cell_id,
            reference, reference_name,
            supplementary_reference_files, supplementary_reference_names,
            num_threads
        )
        all_detailed_counts.append(detailed_metrics)
        all_summary_counts.append(summary_metrics)

        print("Starting Trim Galore")
        helpers.makedirs(os.path.join(lane_tempdir, 'trim_galore'))
        trim_galore_r1 = os.path.join(lane_tempdir, 'trim_galore', 'r1.fastq.gz')
        trim_galore_r2 = os.path.join(lane_tempdir, 'trim_galore', 'r2.fastq.gz')
        trim_galore_temp = os.path.join(lane_tempdir, 'trim_galore')
        trim_galore(
            fastqscreen_r1, fastqscreen_r2, trim_galore_r1, trim_galore_r2,
            adapter1, adapter2, trim_galore_temp, num_threads, run_fastqc=run_fastqc
        )

        print("Starting Alignment")
        helpers.makedirs(os.path.join(lane_tempdir, 'bwa_mem'))
        helpers.makedirs(os.path.join(lane_tempdir, 'bwa_mem', 'alignment_temp_files'))
        lane_aligned_bam = os.path.join(lane_tempdir, 'bwa_mem', 'aligned.bam')
        bwa_align(
            trim_galore_r1, trim_galore_r2, reference, reference_name, metadata_yaml, lane_aligned_bam,
            lane_id, flowcell_id, cell_id, num_threads, os.path.join(lane_tempdir, 'bwa_mem', 'alignment_temp_files')
        )

        print("Tagging Bam with cell id")
        helpers.makedirs(os.path.join(lane_tempdir, 'tagging'))
        lane_tagged_bam = os.path.join(lane_tempdir, 'tagging', 'tagged.bam')
        alignment_utils.tag_bam_with_cellid(
            lane_aligned_bam, lane_tagged_bam, cell_id
        )

        print("Starting Bam Sort")
        helpers.makedirs(os.path.join(lane_tempdir, 'sorting'))
        lane_sorted_bam = os.path.join(lane_tempdir, 'sorting', 'sorted.bam')
        lane_sorted_temp = os.path.join(lane_tempdir, 'sorting')
        bam_sort(lane_tagged_bam, lane_sorted_bam, lane_sorted_temp, num_threads=num_threads, mem='4G')

        final_lane_bams.append(lane_sorted_bam)

    print("Merging all Lanes")
    merge_tempdir = os.path.join(tempdir, cell_id, 'merge')
    bam_merged = os.path.join(tempdir, cell_id, 'merge', 'merged.bam')
    helpers.merge_bams(final_lane_bams, bam_merged, merge_tempdir, num_threads)

    print("reheader")
    reheader_tempdir = os.path.join(tempdir, cell_id, 'reheader')
    reheader_bam = os.path.join(tempdir, cell_id, 'reheader', 'reheader.bam')
    add_cell_id_to_header(bam_merged, cell_id, reheader_bam, reheader_tempdir)

    print("Marking Duplicates")
    tempdir_markdups = os.path.join(tempdir, cell_id, 'markdups')
    metrics_markdups = os.path.join(tempdir, cell_id, 'markdups', 'metrics.txt')
    parsed_markdups = os.path.join(tempdir, cell_id, 'markdups', 'metrics.csv.gz')
    mark_duplicates(
        reheader_bam, bam_output, metrics_markdups, parsed_markdups, tempdir_markdups, cell_id,
        num_threads=num_threads, mem='4G'
    )
    bam_index(bam_output)

    print("Starting Flagstat")
    helpers.makedirs(os.path.join(tempdir, cell_id, 'flagstat'))
    metrics_flagstat = os.path.join(tempdir, cell_id, 'flagstat', 'flagstat.txt')
    parsed_flagstat = os.path.join(tempdir, cell_id, 'flagstat', 'flagstat.csv.gz')
    flagstat(bam_output, metrics_flagstat, parsed_flagstat, cell_id)

    print("Starting Insert Metrics")
    metrics_insert = os.path.join(tempdir, cell_id, 'insert_metrics', 'metrics.txt')
    histogram_insert = os.path.join(tempdir, cell_id, 'insert_metrics', 'histogram.pdf')
    tempdir_insert = os.path.join(tempdir, cell_id, 'insert_metrics')
    parsed_insert = os.path.join(tempdir, cell_id, 'insert_metrics', 'metrics.csv.gz')
    insert_metrics(
        bam_output, metrics_flagstat, metrics_insert, histogram_insert, parsed_insert,
        tempdir_insert, cell_id, num_threads=num_threads
    )

    print("Starting Picard WGS metrics")
    tempdir_wgs = os.path.join(tempdir, cell_id, 'wgs_metrics')
    metrics_wgs = os.path.join(tempdir, cell_id, 'wgs_metrics', 'metrics.txt')
    parsed_wgs = os.path.join(tempdir, cell_id, 'wgs_metrics', 'metrics.csv.gz')
    wgs_metrics(
        bam_output, reference, wgs_metrics_mqual,
        wgs_metrics_bqual, wgs_metrics_count_unpaired,
        metrics_wgs, parsed_wgs, tempdir_wgs, cell_id,
        num_threads=num_threads, mem='4G'
    )

    print("Collecting read attrition metrics")
    helpers.makedirs(os.path.join(tempdir, cell_id, 'read_attrition'))
    read_attrition_metrics = os.path.join(tempdir, cell_id, 'read_attrition', 'metrics.csv.gz')
    get_coverage_metrics(bam_output, cell_id, read_attrition_metrics)

    print('calculating TSS score')
    tss_tempdir = os.path.join(tempdir, cell_id, 'tss_enrichment')
    temp_tss_metrics = os.path.join(tempdir, cell_id, 'tss_enrichment', 'temp_metrics.csv.gz')
    tss_enrichment(
        bam_output, temp_tss_metrics, reference_version, cell_id, tss_tempdir
    )

    print("merging fastqscreen counts")
    detailed_metrics = os.path.join(tempdir, cell_id, 'fastqscreen', 'detailed.csv.gz')
    summary_metrics = os.path.join(tempdir, cell_id, 'fastqscreen', 'summary.csv.gz')
    merge_fastq_screen_counts(
        all_detailed_counts, all_summary_counts, detailed_metrics,
        summary_metrics
    )

    print("merging all bam metrics")
    merged_metrics = os.path.join(tempdir, cell_id, 'fastqscreen', 'merged_metrics.csv.gz')
    csverve.merge_csv(
        [parsed_markdups, parsed_flagstat, parsed_insert, parsed_wgs,
         read_attrition_metrics, temp_tss_metrics, summary_metrics],
        merged_metrics,
        how='outer', on='cell_id', skip_header=False
    )

    is_contaminated_metrics = os.path.join(tempdir, cell_id, 'fastqscreen', 'contaminated_metrics.csv.gz')
    alignment_utils.add_contamination_status(merged_metrics, is_contaminated_metrics, reference_name)
    alignment_utils.add_metadata(is_contaminated_metrics, metadata_yaml, metrics_output)

    print("Collecting GC Metrics")
    tempdir_gc = os.path.join(tempdir, cell_id, 'gc_metrics')
    metrics_gc = os.path.join(tempdir, cell_id, 'gc_metrics', 'metrics.txt')
    summary_gc = os.path.join(tempdir, cell_id, 'gc_metrics', 'summary.txt')
    chart_gc = os.path.join(tempdir, cell_id, 'gc_metrics', 'chart.pdf')
    gc_metrics(
        bam_output, reference, metrics_gc,
        summary_gc, chart_gc, metrics_gc_output,
        tempdir_gc, cell_id, num_threads=num_threads, mem="4G"
    )

    print("building tar file of supplementary metrics")
    tar_dir = os.path.join(tempdir, '{}_metrics'.format(cell_id))
    helpers.makedirs(tar_dir)
    shutil.copy(metrics_markdups, tar_dir)
    shutil.copy(metrics_gc, tar_dir)
    shutil.copy(summary_gc, tar_dir)
    shutil.copy(chart_gc, tar_dir)
    shutil.copy(metrics_flagstat, tar_dir)
    shutil.copy(metrics_insert, tar_dir)
    shutil.copy(histogram_insert, tar_dir)
    shutil.copy(metrics_wgs, tar_dir)
    shutil.copy(detailed_metrics, tar_dir)

    helpers.make_tarfile(tar_output, tar_dir)
