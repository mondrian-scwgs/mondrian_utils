import csverve.api as csverve
import pandas as pd

from mondrianutils.dtypes.alignment import dtypes

import os
import mondrianutils.helpers as helpers
import pysam


def is_valid_tss_error(stdout):
    if 'Can not get any signals' in stdout:
        return True
    if 'Can not get any proper mapped reads' in stdout:
        return True
    if 'Only single end reads' in stdout:
        return True
    return False


def get_tss_score(bamfile, genome_version, tempdir):
    genome_version = genome_version.lower()

    if genome_version not in ('grch37', 'grch38', 'hg18', 'hg19'):
        return float('nan')

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
        if is_valid_tss_error(stdout + stderr):
            tss_score = float('nan')
        else:
            raise Exception(stdout)
    else:
        tss_score = open(tempoutput, 'rt').readlines()
        assert len(tss_score) == 1
        tss_score = tss_score[0]

    return tss_score


def tss_enrichment(bamfile, tss_metrics, genome_version, cell_id, tempdir):
    helpers.makedirs(tempdir)

    score = get_tss_score(bamfile, genome_version, tempdir)

    data = {
        'cell_id': cell_id,
        'tss_enrichment_score': score
    }

    data = pd.DataFrame.from_dict(data, orient='index').T

    csverve.write_dataframe_to_csv_and_yaml(
        data, tss_metrics, dtypes()['metrics'], skip_header=False
    )
