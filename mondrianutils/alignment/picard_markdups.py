import os
import mondrianutils.helpers as helpers
import pandas as pd
import csverve
from mondrianutils.dtypes.alignment import dtypes


def run_picard_markdups(
        bam_filename, markduped_bam_filename,
        metrics_filename, tempdir,
        num_threads=1, mem="2G"
):
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


def extract_duplication_metrics(markdups_metrics, cell_id, parsed_metrics):
    """
    extract from markdups
    """

    mfile = open(markdups_metrics)

    targetlines = []

    line = mfile.readline()

    while line != '':
        if line.startswith('## METRICS CLASS'):
            targetlines.append(mfile.readline().strip('\n').split('\t'))
            targetlines.append(mfile.readline().strip('\n').split('\t'))
            break
        line = mfile.readline()

    mfile.close()

    header, data = targetlines

    header = [v.lower() for v in header]
    header = {v: i for i, v in enumerate(header)}

    unprd_mpd_rds = int(data[header['unpaired_reads_examined']])
    prd_mpd_rds = int(data[header['read_pairs_examined']])
    unprd_dup_rds = int(data[header['unpaired_read_duplicates']])
    prd_dup_rds = int(data[header['read_pair_duplicates']])
    unmpd_rds = data[header['unmapped_reads']]
    est_lib_size = data[header['estimated_library_size']]

    rd_pair_opt_dup = int(data[header['read_pair_optical_duplicates']])

    try:
        perc_dup_reads = (unprd_dup_rds + ((prd_dup_rds + rd_pair_opt_dup) * 2)) / (
                unprd_mpd_rds + (prd_mpd_rds * 2))
    except ZeroDivisionError:
        perc_dup_reads = 0

    outdata = {
        'cell_id': cell_id,
        'unpaired_mapped_reads': unprd_mpd_rds,
        'paired_mapped_reads': prd_mpd_rds,
        'unpaired_duplicate_reads': unprd_dup_rds,
        'paired_duplicate_reads': prd_dup_rds,
        'unmapped_reads': unmpd_rds,
        'percent_duplicate_reads': perc_dup_reads,
        'estimated_library_size': est_lib_size,
    }

    outdata = {k: 0 if v == '' else v for k, v in outdata.items()}

    outdata = pd.DataFrame.from_dict(outdata, orient='index').T

    csverve.write_dataframe_to_csv_and_yaml(
        outdata, parsed_metrics, dtypes()['metrics'], skip_header=False)

def mark_duplicates(
        bam_filename, markdups_bam,
        markdups_metrics, parsed_metrics, tempdir, cell_id,
        num_threads=1, mem="2G"
):
    run_picard_markdups(
        bam_filename, markdups_bam,
        markdups_metrics, tempdir,
        num_threads=num_threads, mem=mem
    )

    extract_duplication_metrics(markdups_metrics, cell_id, parsed_metrics)
