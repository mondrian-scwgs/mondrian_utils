import json
import os
import pathlib
import shutil

import mondrianutils.helpers as helpers
import pysam
import yaml
from mondrianutils.alignment.collect_gc_metrics import collect_gc_metrics
from mondrianutils.alignment.collect_metrics import collect_metrics
from mondrianutils.alignment.coverage_metrics import get_coverage_metrics
from mondrianutils.alignment.fastqscreen import merge_fastq_screen_counts
from mondrianutils.alignment.fastqscreen import organism_filter
from mondrianutils.alignment.trim_galore import trim_galore


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
        fastq1, fastq2, reference, metadata_yaml,
        output, lane_id, flowcell_id, cell_id
):
    center, library_id, sample_id = load_metadata(metadata_yaml, lane_id, flowcell_id, cell_id)

    script_path = pathlib.Path(__file__).parent.resolve()
    script_path = os.path.join(script_path, 'bwa_align.sh')

    cmd = [
        script_path, fastq1, fastq2, reference, output,
        sample_id, library_id, cell_id, lane_id,
        flowcell_id, center
    ]
    helpers.run_cmd(cmd)


def tag_bam_with_cellid(infile, outfile, cell_id):
    infile = pysam.AlignmentFile(infile, "rb")
    outfile = pysam.AlignmentFile(outfile, "wb", template=infile)

    iter = infile.fetch(until_eof=True)
    for read in iter:
        read.set_tag('CB', cell_id, replace=False)
        outfile.write(read)
    infile.close()
    outfile.close()


def bam_sort(bam_filename, sorted_bam_filename, tempdir, mem="2G"):
    if not os.path.exists(tempdir):
        helpers.makedirs(tempdir)

    cmd = [
        'picard', '-Xmx' + mem, '-Xms' + mem,
        '-XX:ParallelGCThreads=1',
        'SortSam',
                  'INPUT=' + bam_filename,
                  'OUTPUT=' + sorted_bam_filename,
        'SORT_ORDER=coordinate',
        'VALIDATION_STRINGENCY=LENIENT',
                  'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        'QUIET=true'
    ]
    helpers.run_cmd(cmd)


def merge_bams(inputs, output, mem="2G"):
    if isinstance(inputs, dict):
        inputs = inputs.values()

    cmd = ['picard', '-Xmx' + mem, '-Xms' + mem,
           '-XX:ParallelGCThreads=1',
           'MergeSamFiles',
           'OUTPUT=' + output,
           'SORT_ORDER=coordinate',
           'ASSUME_SORTED=true',
           'VALIDATION_STRINGENCY=LENIENT',
           'MAX_RECORDS_IN_RAM=150000',
           'QUIET=true'
           ]

    for bamfile in inputs:
        cmd.append('I=' + os.path.abspath(bamfile))

    helpers.run_cmd(cmd)


def bam_markdups(bam_filename, markduped_bam_filename, metrics_filename,
                 tempdir, mem="2G"):
    if not os.path.exists(tempdir):
        helpers.makedirs(tempdir)

    helpers.run_cmd([
        'picard', '-Xmx' + mem, '-Xms' + mem,
        '-XX:ParallelGCThreads=1',
        'MarkDuplicates',
                  'INPUT=' + bam_filename,
                  'OUTPUT=' + markduped_bam_filename,
                  'METRICS_FILE=' + metrics_filename,
        'REMOVE_DUPLICATES=False',
        'ASSUME_SORTED=True',
        'VALIDATION_STRINGENCY=LENIENT',
                  'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        'QUIET=true'
    ])


def bam_collect_gc_metrics(bam_filename, ref_genome, metrics_filename,
                           summary_filename, chart_filename, tempdir,
                           mem="2G"):
    if not os.path.exists(tempdir):
        helpers.makedirs(tempdir)

    helpers.run_cmd([
        'picard', '-Xmx' + mem, '-Xms' + mem,
        '-XX:ParallelGCThreads=1',
        'CollectGcBiasMetrics',
                  'INPUT=' + bam_filename,
                  'OUTPUT=' + metrics_filename,
                  'REFERENCE_SEQUENCE=' + ref_genome,
                  'S=' + summary_filename,
                  'CHART_OUTPUT=' + chart_filename,
        'VALIDATION_STRINGENCY=LENIENT',
                  'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        'QUIET=true'
    ])


def bam_collect_insert_metrics(bam_filename, flagstat_metrics_filename,
                               metrics_filename, histogram_filename, tempdir,
                               mem="2G"):
    # Check if any paired reads exist
    has_paired = None
    with open(flagstat_metrics_filename) as f:
        for line in f:
            if 'properly paired' in line:
                if line.startswith('0 '):
                    has_paired = False
                else:
                    has_paired = True

    if has_paired is None:
        raise Exception('Unable to determine number of properly paired reads from {}'.format(
            flagstat_metrics_filename))

    if not has_paired:
        with open(metrics_filename, 'w') as f:
            f.write('## FAILED: No properly paired reads\n')
        with open(histogram_filename, 'w'):
            pass
        return

    if not os.path.exists(tempdir):
        helpers.makedirs(tempdir)

    helpers.run_cmd([
        'picard', '-Xmx' + mem, '-Xms' + mem,
        '-XX:ParallelGCThreads=1',
        'CollectInsertSizeMetrics',
                  'INPUT=' + bam_filename,
                  'OUTPUT=' + metrics_filename,
                  'HISTOGRAM_FILE=' + histogram_filename,
        'ASSUME_SORTED=True',
        'VALIDATION_STRINGENCY=LENIENT',
                  'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        'QUIET=true'
    ])


def bam_collect_wgs_metrics(bam_filename, ref_genome, metrics_filename,
                            mqual, bqual, count_unpaired, tempdir, mem="2G"):
    if not os.path.exists(tempdir):
        helpers.makedirs(tempdir)

    helpers.run_cmd([
        'picard', '-Xmx' + mem, '-Xms' + mem,
        '-XX:ParallelGCThreads=1',
        'CollectWgsMetrics',
                  'INPUT=' + bam_filename,
                  'OUTPUT=' + metrics_filename,
                  'REFERENCE_SEQUENCE=' + ref_genome,
                  'MINIMUM_BASE_QUALITY=' +
                  bqual,
                  'MINIMUM_MAPPING_QUALITY=' +
                  mqual,
        'COVERAGE_CAP=500',
        'VALIDATION_STRINGENCY=LENIENT',
                  'COUNT_UNPAIRED=' +
                  ('True' if count_unpaired else 'False'),
                  'TMP_DIR=' + tempdir,
        'MAX_RECORDS_IN_RAM=150000',
        'QUIET=true'
    ])


def bam_flagstat(bam, metrics):
    script_path = pathlib.Path(__file__).parent.resolve()
    script_path = os.path.join(script_path, 'flagstat.sh')
    helpers.run_cmd([
        script_path,
        bam,
        metrics,
    ])


def bam_index(infile):
    helpers.run_cmd([
        'samtools', 'index',
        infile,
        infile + '.bai'
    ])


def alignment(
        fastq_files, metadata_yaml, reference, reference_name, supplementary_references, tempdir,
        adapter1, adapter2, cell_id, wgs_metrics_mqual, wgs_metrics_bqual, wgs_metrics_count_unpaired,
        bam_output, metrics_output, metrics_gc_output, fastqscreen_detailed_output, fastqscreen_summary_output,
        tar_output
):

    with open(supplementary_references, 'rt') as reader:
        supplementary_references = json.load(reader)

    with open(fastq_files, 'rt') as reader:
        fastqdata = json.load(reader)

    final_lane_bams = []
    all_detailed_counts = []
    all_summary_counts = []
    for lane in fastqdata:
        r1 = lane['fastq1']
        r2 = lane['fastq2']
        lane_id = lane['lane_id']
        flowcell_id = lane['flowcell_id']

        helpers.makedirs(os.path.join(tempdir, lane_id, 'fastqscreen'))
        fastqscreen_r1 = os.path.join(tempdir, lane_id, 'fastqscreen', 'r1.fastq.gz')
        fastqscreen_r2 = os.path.join(tempdir, lane_id, 'fastqscreen', 'r2.fastq.gz')
        detailed_metrics = os.path.join(tempdir, lane_id, 'fastqscreen', 'detailed.csv.gz')
        summary_metrics = os.path.join(tempdir, lane_id, 'fastqscreen', 'summary.csv.gz')
        fastqscreen_temp = os.path.join(tempdir, lane_id, 'fastqscreen')

        organism_filter(
            r1, r2, fastqscreen_r1, fastqscreen_r2,
            detailed_metrics, summary_metrics, fastqscreen_temp,
            cell_id, reference, reference_name, supplementary_references
        )
        all_detailed_counts.append(detailed_metrics)
        all_summary_counts.append(summary_metrics)

        helpers.makedirs(os.path.join(tempdir, lane_id, 'trim_galore'))
        trim_galore_r1 = os.path.join(tempdir, lane_id, 'trim_galore', 'r1.fastq.gz')
        trim_galore_r2 = os.path.join(tempdir, lane_id, 'trim_galore', 'r2.fastq.gz')
        trim_galore_temp = os.path.join(tempdir, lane_id, 'trim_galore')
        trim_galore(
            fastqscreen_r1, fastqscreen_r2, trim_galore_r1, trim_galore_r2,
            adapter1, adapter2, trim_galore_temp
        )

        helpers.makedirs(os.path.join(tempdir, lane_id, 'bwa_mem'))
        lane_aligned_bam = os.path.join(tempdir, lane_id, 'bwa_mem', 'aligned.bam')
        bwa_align(
            trim_galore_r1, trim_galore_r2, reference, metadata_yaml, lane_aligned_bam,
            lane_id, flowcell_id, cell_id
        )

        helpers.makedirs(os.path.join(tempdir, lane_id, 'tagging'))
        lane_tagged_bam = os.path.join(tempdir, lane_id, 'tagging', 'tagged.bam')
        tag_bam_with_cellid(
            lane_aligned_bam, lane_tagged_bam, cell_id
        )

        helpers.makedirs(os.path.join(tempdir, lane_id, 'sorting'))
        lane_sorted_bam = os.path.join(tempdir, lane_id, 'sorting', 'sorted.bam')
        lane_sorted_temp = os.path.join(tempdir, lane_id, 'sorting')
        bam_sort(lane_tagged_bam, lane_sorted_bam, lane_sorted_temp, mem='4G')

        final_lane_bams.append(lane_sorted_bam)

    helpers.makedirs(os.path.join(tempdir, cell_id, 'merge'))
    bam_merged = os.path.join(tempdir, cell_id, 'merge', 'merged.bam')
    merge_bams(final_lane_bams, bam_merged, mem='4G')

    helpers.makedirs(os.path.join(tempdir, cell_id, 'markdups'))
    metrics_markdups = os.path.join(tempdir, cell_id, 'markdups', 'metrics.txt')
    tempdir_markdups = os.path.join(tempdir, cell_id, 'markdups')
    bam_markdups(bam_merged, bam_output, metrics_markdups, tempdir_markdups, mem='4G')
    bam_index(bam_output)

    helpers.makedirs(os.path.join(tempdir, cell_id, 'gc_metrics'))
    metrics_gc = os.path.join(tempdir, cell_id, 'gc_metrics', 'metrics.txt')
    summary_gc = os.path.join(tempdir, cell_id, 'gc_metrics', 'summary.txt')
    chart_gc = os.path.join(tempdir, cell_id, 'gc_metrics', 'chart.pdf')
    tempdir_gc = os.path.join(tempdir, cell_id, 'gc_metrics')
    bam_collect_gc_metrics(
        bam_output, reference, metrics_gc,
        summary_gc, chart_gc, tempdir_gc, mem="4G"
    )

    helpers.makedirs(os.path.join(tempdir, cell_id, 'flagstat'))
    metrics_flagstat = os.path.join(tempdir, cell_id, 'flagstat', 'flagstat.txt')
    bam_flagstat(bam_output, metrics_flagstat)

    helpers.makedirs(os.path.join(tempdir, cell_id, 'insert_metrics'))
    metrics_insert = os.path.join(tempdir, cell_id, 'insert_metrics', 'metrics.txt')
    histogram_insert = os.path.join(tempdir, cell_id, 'insert_metrics', 'histogram.pdf')
    tempdir_insert = os.path.join(tempdir, cell_id, 'insert_metrics')
    bam_collect_insert_metrics(
        bam_output, metrics_flagstat, metrics_insert, histogram_insert, tempdir_insert
    )

    helpers.makedirs(os.path.join(tempdir, cell_id, 'wgs_metrics'))
    metrics_wgs = os.path.join(tempdir, cell_id, 'wgs_metrics', 'metrics.txt')
    tempdir_wgs = os.path.join(tempdir, cell_id, 'wgs_metrics')
    bam_collect_wgs_metrics(
        bam_output, reference, metrics_wgs, wgs_metrics_mqual,
        wgs_metrics_bqual, wgs_metrics_count_unpaired, tempdir_wgs, mem='4G'
    )

    helpers.makedirs(os.path.join(tempdir, cell_id, 'read_attrition'))
    read_attrition_metrics = os.path.join(tempdir, cell_id, 'read_attrition', 'metrics.csv.gz')
    get_coverage_metrics(bam_output, read_attrition_metrics)

    collect_metrics(
        metrics_wgs, metrics_insert, metrics_flagstat, metrics_markdups,
        read_attrition_metrics, metrics_output, cell_id
    )

    collect_gc_metrics(metrics_gc, metrics_gc_output, cell_id)

    merge_fastq_screen_counts(
        all_detailed_counts, all_summary_counts, fastqscreen_detailed_output,
        fastqscreen_summary_output
    )

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

    helpers.make_tarfile(tar_output, tar_dir)
