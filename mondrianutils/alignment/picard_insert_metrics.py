import os
import mondrianutils.helpers as helpers
import pandas as pd
import csverve
from mondrianutils.dtypes.alignment import dtypes


def bam_collect_insert_metrics(
        bam_filename, flagstat_metrics_filename,
        metrics_filename, histogram_filename, tempdir,
        num_threads=1, mem="2G"
):
    helpers.makedirs(tempdir)
    helpers.makedirs(os.path.basename(metrics_filename))
    helpers.makedirs(os.path.basename(histogram_filename))

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


def extract_insert_metrics(insert_metrics):
    ''' Extract median and mean insert size '''

    # picardtools insertmetrics completes with code 0 and doesn't generate metrics file
    # if inputs don't have sufficient read count
    if not os.path.isfile(insert_metrics):
        return 0, 0, 0

    # if the insert metrics fails due to low coverage
    if open(insert_metrics).readline().startswith("## FAILED"):
        return 0, 0, 0

    mfile = open(insert_metrics)

    targetlines = []

    line = mfile.readline()

    while line != '':
        if line.startswith('## METRICS CLASS'):
            targetlines.append(mfile.readline().strip().split('\t'))
            targetlines.append(mfile.readline().strip().split('\t'))
            break
        line = mfile.readline()

    mfile.close()

    header, data = targetlines

    header = [v.lower() for v in header]
    header = {v: i for i, v in enumerate(header)}

    median_ins_size = data[header['median_insert_size']]
    mean_ins_size = data[header['mean_insert_size']]
    std_dev_ins_size = data[header['standard_deviation']]

    median_ins_size = 0 if median_ins_size == '?' else median_ins_size
    mean_ins_size = 0 if mean_ins_size == '?' else mean_ins_size
    std_dev_ins_size = 0 if std_dev_ins_size == '?' else std_dev_ins_size

    return median_ins_size, mean_ins_size, std_dev_ins_size


def parse_insert_metrics(insert_metrics, cell_id, parsed_metrics):
    median_ins_size, mean_ins_size, std_dev_ins_size = extract_insert_metrics(insert_metrics)

    outdata = {
        'cell_id': cell_id,
        'median_insert_size': median_ins_size,
        'mean_insert_size': mean_ins_size,
        'standard_deviation_insert_size': std_dev_ins_size,
    }

    outdata = pd.DataFrame.from_dict(outdata, orient='index').T

    csverve.write_dataframe_to_csv_and_yaml(
        outdata, parsed_metrics, dtypes()['metrics'], skip_header=False)


def insert_metrics(
        bam, flagstat_metrics,
        insert_metrics, histogram, parsed_metrics,
        tempdir, cell_id, num_threads=1, mem="2G"
):
    bam_collect_insert_metrics(
        bam, flagstat_metrics, insert_metrics, histogram, tempdir,
        num_threads=num_threads, mem=mem
    )

    parse_insert_metrics(insert_metrics, cell_id, parsed_metrics)
