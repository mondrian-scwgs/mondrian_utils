'''
Created on Sep 8, 2015

@author: dgrewal
'''
import logging
import csverve.api as csverve
import pandas as pd

from mondrianutils.dtypes.alignment import dtypes

import os
import mondrianutils.helpers as helpers


def bam_collect_gc_metrics(
        bam_filename, ref_genome, metrics_filename,
        summary_filename, chart_filename, tempdir,
        num_threads=1, mem="2G"
):
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


def collect_gc_metrics(bias_file, output, sample_id):
    """
    parses the gcbias data
    """

    data = open(bias_file).readlines()
    skiprows = [i for i, v in enumerate(data) if v[0] == '#' or v == '\n']

    # If the file is empty (only header no data) then return 0s (dummy data)
    try:
        data = pd.read_csv(bias_file, sep='\t', skiprows=skiprows)
        data = data[['NORMALIZED_COVERAGE', 'WINDOWS', 'GC']]
        data = data.set_index('GC')
        data = data.rename(columns={'NORMALIZED_COVERAGE': sample_id, 'WINDOWS': 'reference'})
        data = data.T
        data['cell_id'] = data.index
        data.columns = data.columns.astype(str)
        data = data.replace('?', 0)
    except pd.errors.EmptyDataError:
        logging.getLogger("single_cell.align.gcbias").warn(
            'No data in the GCBias output')
        # If the file is empty (only header no data) then return 0s (dummy data)

        dummy_data = [[0 for v in range(101)] + [sample_id]]
        columns = [str(v) for v in range(101)] + ['cell_id']
        data = pd.DataFrame(dummy_data, columns=columns)

    csverve.write_dataframe_to_csv_and_yaml(data, output, dtypes()['gc'], skip_header=False)


def gc_metrics(
        bam_filename, ref_genome, metrics_filename,
        summary_filename, chart_filename, parsed_metrics,
        tempdir, cell_id, num_threads=1, mem="2G"
):
    bam_collect_gc_metrics(
        bam_filename, ref_genome, metrics_filename,
        summary_filename, chart_filename, tempdir,
        num_threads=num_threads, mem=mem
    )

    collect_gc_metrics(metrics_filename, parsed_metrics, cell_id)
