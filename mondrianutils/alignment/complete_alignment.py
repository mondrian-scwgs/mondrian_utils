import json
import os
import pathlib
import shutil
import csverve

import mondrianutils.helpers as helpers
import pysam
import yaml
from mondrianutils.alignment.collect_gc_metrics import collect_gc_metrics
from mondrianutils.alignment.collect_metrics import collect_metrics
from mondrianutils.alignment.coverage_metrics import get_coverage_metrics
from mondrianutils.alignment.fastqscreen import merge_fastq_screen_counts
from mondrianutils.alignment.fastqscreen import organism_filter
from mondrianutils.alignment.trim_galore import trim_galore

from mondrianutils.dtypes.alignment import dtypes


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


def tag_bam_with_cellid(infile, outfile, cell_id):
    infile = pysam.AlignmentFile(infile, "rb")
    outfile = pysam.AlignmentFile(outfile, "wb", template=infile)

    iter = infile.fetch(until_eof=True)
    for read in iter:
        read.set_tag('CB', cell_id, replace=False)
        outfile.write(read)
    infile.close()
    outfile.close()


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


def merge_bams(inputs, output, num_threads=1, mem="2G"):
    if isinstance(inputs, dict):
        inputs = inputs.values()

    cmd = ['picard', '-Xmx' + mem, '-Xms' + mem]
    if num_threads == 1:
        cmd.append('-XX:ParallelGCThreads=1')
    cmd.extend([
        'MergeSamFiles',
        'OUTPUT=' + output,
        'SORT_ORDER=coordinate',
        'ASSUME_SORTED=true',
        'VALIDATION_STRINGENCY=LENIENT',
        'MAX_RECORDS_IN_RAM=150000',
        'QUIET=true'
    ])

    for bamfile in inputs:
        cmd.append('I=' + os.path.abspath(bamfile))

    if not num_threads == 1:
        cmd.append('USE_THREADING=true')

    helpers.run_cmd(cmd)


def bam_markdups(bam_filename, markduped_bam_filename, metrics_filename,
                 tempdir, num_threads=1, mem="2G"):
    if not os.path.exists(tempdir):
        helpers.makedirs(tempdir)

    cmd = ['picard', '-Xmx' + mem, '-Xms' + mem]
    if num_threads == 1:
        cmd.append('-XX:ParallelGCThreads=1')
    cmd.extend([
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
    helpers.run_cmd(cmd)


def bam_collect_gc_metrics(bam_filename, ref_genome, metrics_filename,
                           summary_filename, chart_filename, tempdir,
                           num_threads=1, mem="2G"):
    if not os.path.exists(tempdir):
        helpers.makedirs(tempdir)

    cmd = ['picard', '-Xmx' + mem, '-Xms' + mem]
    if num_threads == 1:
        cmd.append('-XX:ParallelGCThreads=1')
    cmd.extend([
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
    helpers.run_cmd(cmd)


def bam_collect_insert_metrics(
        bam_filename, flagstat_metrics_filename,
        metrics_filename, histogram_filename, tempdir,
        num_threads=1, mem="2G"
):
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

    cmd = ['picard', '-Xmx' + mem, '-Xms' + mem]
    if num_threads == 1:
        cmd.append('-XX:ParallelGCThreads=1')
    cmd.extend([
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
    helpers.run_cmd(cmd)


def bam_collect_wgs_metrics(bam_filename, ref_genome, metrics_filename,
                            mqual, bqual, count_unpaired, tempdir, num_threads=1, mem="2G"):
    if not os.path.exists(tempdir):
        helpers.makedirs(tempdir)

    cmd = ['picard', '-Xmx' + mem, '-Xms' + mem]
    if num_threads == 1:
        cmd.append('-XX:ParallelGCThreads=1')
    cmd.extend([
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
    helpers.run_cmd(cmd)


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

def is_valid_tss_error(stdout):
    if 'Can not get any signals' in stdout:
        return True
    if 'Can not get any proper mapped reads' in stdout:
        return True
    return False

def add_tss_enrichment(bamfile, metrics_file, annotated_metrics, genome_version, tempdir):
    genome_version = genome_version.lower()

    if genome_version not in ('grch37', 'grch38', 'hg18', 'hg19'):
        csverve.simple_annotate_csv(
            metrics_file,
            annotated_metrics,
            'tss_enrichment_score',
            'NA',
            dtypes()['metrics']['tss_enrichment_score']
        )
        return

    with pysam.AlignmentFile(bamfile, "rb") as bam_reader:
        bam_chromosomes = bam_reader.references

    chromosomes = [str(v) for v in range(1, 23)] + ['X']
    if bam_chromosomes[0].startswith('chr'):
        chromosomes = ['chr' + v for v in chromosomes]
    chromosomes = [v for v in chromosomes if v in bam_chromosomes]
    assert len(chromosomes) > 1

    helpers.makedirs(tempdir)

    tempoutput = os.path.join(tempdir, 'temp_tss_output.txt')

    scripts_directory = os.path.realpath(os.path.dirname(__file__))
    run_rscript = os.path.join(scripts_directory, 'tss_enrichment.R')
    cmd = [
        run_rscript,
        '--bamfile', bamfile,
        '--output', tempoutput,
        '--tempdir', tempdir,
        '--genome_version', genome_version,
        '--chromosomes', ','.join(chromosomes)
    ]

    stdout, stderr = helpers.run_cmd(cmd)

    if not os.path.exists(tempoutput):
        if is_valid_tss_error(stdout):
            tss_score = float('nan')
        else:
            raise Exception(stdout)
    else:
        tss_score = open(tempoutput, 'rt').readlines()
        assert len(tss_score) == 1
        tss_score = tss_score[0]

    csverve.simple_annotate_csv(
        metrics_file,
        annotated_metrics,
        'tss_enrichment_score',
        tss_score,
        dtypes()['metrics']['tss_enrichment_score']
    )


def alignment(
        fastq_files, metadata_yaml, reference, reference_name, reference_version, supplementary_references, tempdir,
        adapter1, adapter2, cell_id, wgs_metrics_mqual, wgs_metrics_bqual, wgs_metrics_count_unpaired,
        bam_output, metrics_output, metrics_gc_output, fastqscreen_detailed_output, fastqscreen_summary_output,
        tar_output, num_threads, run_fastqc=False
):
    with open(supplementary_references, 'rt') as reader:
        supplementary_references = json.load(reader)

    with open(fastq_files, 'rt') as reader:
        fastqdata = json.load(reader)

    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)

    final_lane_bams = []
    all_detailed_counts = []
    all_summary_counts = []
    for lane in fastqdata:
        r1 = lane['fastq1']
        r2 = lane['fastq2']
        lane_id = lane['lane_id']
        flowcell_id = lane['flowcell_id']

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
            detailed_metrics, summary_metrics, fastqscreen_temp,
            cell_id, reference, reference_name, supplementary_references, num_threads
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
        tag_bam_with_cellid(
            lane_aligned_bam, lane_tagged_bam, cell_id
        )

        print("Starting Bam Sort")
        helpers.makedirs(os.path.join(lane_tempdir, 'sorting'))
        lane_sorted_bam = os.path.join(lane_tempdir, 'sorting', 'sorted.bam')
        lane_sorted_temp = os.path.join(lane_tempdir, 'sorting')
        bam_sort(lane_tagged_bam, lane_sorted_bam, lane_sorted_temp, num_threads=num_threads, mem='4G')

        final_lane_bams.append(lane_sorted_bam)

    print("Merging all Lanes")
    helpers.makedirs(os.path.join(tempdir, cell_id, 'merge'))
    bam_merged = os.path.join(tempdir, cell_id, 'merge', 'merged.bam')
    merge_bams(final_lane_bams, bam_merged, num_threads=num_threads, mem='4G')

    print("Marking Duplicates")
    helpers.makedirs(os.path.join(tempdir, cell_id, 'markdups'))
    metrics_markdups = os.path.join(tempdir, cell_id, 'markdups', 'metrics.txt')
    tempdir_markdups = os.path.join(tempdir, cell_id, 'markdups')
    bam_markdups(bam_merged, bam_output, metrics_markdups, tempdir_markdups, num_threads=num_threads, mem='4G')
    bam_index(bam_output)

    print("Collecting GC Metrics")
    helpers.makedirs(os.path.join(tempdir, cell_id, 'gc_metrics'))
    metrics_gc = os.path.join(tempdir, cell_id, 'gc_metrics', 'metrics.txt')
    summary_gc = os.path.join(tempdir, cell_id, 'gc_metrics', 'summary.txt')
    chart_gc = os.path.join(tempdir, cell_id, 'gc_metrics', 'chart.pdf')
    tempdir_gc = os.path.join(tempdir, cell_id, 'gc_metrics')
    bam_collect_gc_metrics(
        bam_output, reference, metrics_gc,
        summary_gc, chart_gc, tempdir_gc, num_threads=num_threads, mem="4G"
    )

    print("Starting Flagstat")
    helpers.makedirs(os.path.join(tempdir, cell_id, 'flagstat'))
    metrics_flagstat = os.path.join(tempdir, cell_id, 'flagstat', 'flagstat.txt')
    bam_flagstat(bam_output, metrics_flagstat)

    print("Starting Insert Metrics")
    helpers.makedirs(os.path.join(tempdir, cell_id, 'insert_metrics'))
    metrics_insert = os.path.join(tempdir, cell_id, 'insert_metrics', 'metrics.txt')
    histogram_insert = os.path.join(tempdir, cell_id, 'insert_metrics', 'histogram.pdf')
    tempdir_insert = os.path.join(tempdir, cell_id, 'insert_metrics')
    bam_collect_insert_metrics(
        bam_output, metrics_flagstat, metrics_insert, histogram_insert, tempdir_insert, num_threads=num_threads
    )

    print("Starting Picard WGS metrics")
    helpers.makedirs(os.path.join(tempdir, cell_id, 'wgs_metrics'))
    metrics_wgs = os.path.join(tempdir, cell_id, 'wgs_metrics', 'metrics.txt')
    tempdir_wgs = os.path.join(tempdir, cell_id, 'wgs_metrics')
    bam_collect_wgs_metrics(
        bam_output, reference, metrics_wgs, wgs_metrics_mqual,
        wgs_metrics_bqual, wgs_metrics_count_unpaired, tempdir_wgs, num_threads=num_threads, mem='4G'
    )

    print("Collecting read attrition metrics")
    helpers.makedirs(os.path.join(tempdir, cell_id, 'read_attrition'))
    read_attrition_metrics = os.path.join(tempdir, cell_id, 'read_attrition', 'metrics.csv.gz')
    get_coverage_metrics(bam_output, read_attrition_metrics)

    print("parsing all bam metrics")
    helpers.makedirs(os.path.join(tempdir, cell_id, 'metrics'))
    temp_collect_metrics = os.path.join(tempdir, cell_id, 'metrics', 'temp_metrics.csv.gz')
    collect_metrics(
        metrics_wgs, metrics_insert, metrics_flagstat, metrics_markdups,
        read_attrition_metrics, temp_collect_metrics, cell_id
    )
    helpers.makedirs(os.path.join(tempdir, cell_id, 'tss_enrichment'))
    add_tss_enrichment(
        bam_output, temp_collect_metrics, metrics_output, reference_version,
        os.path.join(tempdir, cell_id, 'tss_enrichment')
    )

    print("parsing GC metrics")
    collect_gc_metrics(metrics_gc, metrics_gc_output, cell_id)

    print("merging fastqscreen counts")
    merge_fastq_screen_counts(
        all_detailed_counts, all_summary_counts, fastqscreen_detailed_output,
        fastqscreen_summary_output
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

    helpers.make_tarfile(tar_output, tar_dir)
