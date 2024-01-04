import os
import mondrianutils.helpers as helpers
import pandas as pd
import csverve
from mondrianutils.dtypes.alignment import dtypes
import pathlib


def bam_flagstat(bam, metrics):
    script_path = pathlib.Path(__file__).parent.resolve()
    script_path = os.path.join(script_path, 'flagstat.sh')
    helpers.run_cmd([
        script_path,
        bam,
        metrics,
    ])


def extract_flagstat_metrics(flagstat_metrics, cell_id, parsed_metrics):
    """
    extract from flagstat
    """

    df = pd.read_csv(
        flagstat_metrics,
        sep=r'\s\+\s0\s',
        header=None,
        names=['value', 'type'],
        engine='python'
    )

    tot_reads = df[df['type'] == 'in total (QC-passed reads + QC-failed reads)']['value']
    tot_mpd_reads = df[
        (df['type'].str.contains('mapped') == True) &
        (df['type'].str.contains('primary mapped') == False) &
        (df['type'].str.contains('mate mapped') == False)
        ]
    tot_dup_reads = df[df['type'] == 'duplicates']['value']
    tot_prop_paired = df[df['type'].str.contains('properly paired')]

    assert len(tot_reads) == 1
    assert len(tot_mpd_reads) == 1
    assert len(tot_dup_reads) == 1
    assert len(tot_prop_paired) == 1

    tot_reads = tot_reads.iloc[0]
    tot_mpd_reads = tot_mpd_reads['value'].iloc[0]
    tot_dup_reads = tot_dup_reads.iloc[0]
    tot_prop_paired = tot_prop_paired['value'].iloc[0]

    outdata = {
        'cell_id': cell_id,
        'total_reads': tot_reads,
        'total_mapped_reads': tot_mpd_reads,
        'total_duplicate_reads': tot_dup_reads,
        'total_properly_paired': tot_prop_paired
    }

    outdata = pd.DataFrame.from_dict(outdata, orient='index').T

    csverve.write_dataframe_to_csv_and_yaml(
        outdata, parsed_metrics, dtypes()['metrics'], skip_header=False
    )


def flagstat(bam, flagstat_metrics, parsed_metrics, cell_id):
    bam_flagstat(bam, flagstat_metrics)
    extract_flagstat_metrics(flagstat_metrics, cell_id, parsed_metrics)
