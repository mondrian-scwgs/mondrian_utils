import csverve.api as csverve
import pandas as pd

from mondrianutils.dtypes.alignment import dtypes

import os
import mondrianutils.helpers as helpers


def bam_collect_wgs_metrics(
        bam_filename, ref_genome, metrics_filename,
        mqual, bqual, count_unpaired, tempdir,
        num_threads=1, mem="2G"
):
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


def extract_wgs_metrics(wgs_metrics, cell_id, parsed_output):
    """
    get the coverage_depth (mean_coverage column)
    get the coverage_breadth (count/genome_territory)
    """

    mfile = open(wgs_metrics)

    metrics = []
    hist = {}

    addmetrics = False
    addhist = False

    for line in mfile:
        if line.strip() == '':
            continue
        if line.startswith('## METRICS CLASS'):
            addmetrics = True
            addhist = False
            continue

        if line.startswith('## HISTOGRAM'):
            addhist = True
            addmetrics = False
            continue

        if addmetrics:
            metrics.append(line.strip().split('\t'))
        if addhist:
            line = line.strip().split('\t')
            if line[0] == 'coverage':
                continue
            hist[int(line[0])] = int(line[1])

    mfile.close()
    header, data = metrics

    header = [v.lower() for v in header]
    header = {v: i for i, v in enumerate(header)}

    gen_territory = int(data[header['genome_territory']])
    cov_depth = float(data[header['mean_coverage']])
    count = int(hist[0])
    cov_breadth = (gen_territory - count) / gen_territory

    outdata = {
        'cell_id': cell_id,
        'coverage_depth': cov_depth,
        'coverage_breadth': cov_breadth,
    }

    outdata = pd.DataFrame.from_dict(outdata, orient='index').T

    csverve.write_dataframe_to_csv_and_yaml(
        outdata, parsed_output, dtypes()['metrics'], skip_header=False

    )


def wgs_metrics(
        bam_filename, ref_genome,
        mqual, bqual, count_unpaired,
        metrics_filename, parsed_metrics,
        tempdir, cell_id,
        num_threads=1, mem="2G"

):
    bam_collect_wgs_metrics(
        bam_filename, ref_genome, metrics_filename,
        mqual, bqual, count_unpaired, tempdir,
        num_threads=num_threads, mem=mem
    )

    extract_wgs_metrics(metrics_filename, cell_id, parsed_metrics)
