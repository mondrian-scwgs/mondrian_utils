import os

import csverve.api as csverve
import pandas as pd
from mondrianutils.dtypes.breakpoint import dtypes

from .breakpoint_db import BreakpointDatabase
from .vcf_sv_parser import SvVcfData


def parse_region(region):
    if region is None:
        return None, None, None

    if ':' not in region:
        return region, None, None

    chrom, coords = region.split(':')

    assert '-' in coords

    beg, end = coords.split('-')

    return chrom, int(beg), int(end)


def read_destruct(destruct_calls, region=None):
    df = pd.read_csv(destruct_calls, sep='\t', dtype={'chromosome_1': str, 'chromosome_2': str})

    df = df[
        ['prediction_id', 'chromosome_1', 'position_1', 'strand_1', 'chromosome_2', 'position_2', 'strand_2', 'type']]
    df['breakpoint_id'] = df.prediction_id
    del df['prediction_id']
    df['caller'] = 'destruct'

    df['breakpoint_id'] = df['breakpoint_id'].astype(str) + '_' + df['caller']

    if region is not None:
        chromosome, start, end = parse_region(region)

        if chromosome:
            query_str = f'(chromosome_1 == "{chromosome}") | (chromosome_2 == "{chromosome}")'
            df = df.query(query_str)

        if start:
            query_str = f'({start} <= position_1 <= {end}) | ({start} <= position_2 <= {end})'
            df = df.query(query_str)

    return df


def check_common(x, df_db, calls):
    val = df_db.query(x, extend=500)

    val = sorted(val)

    if len(val) == 1:
        return

    if val[0] not in calls:
        calls[val[0]] = set()

    for v in val[1:]:
        calls[val[0]].add(v)


def get_common_calls(df, df_db):
    calls = {}

    for i, row in df.iterrows():
        check_common(row, df_db, calls)

    new_groups = {}
    for i, (key, vals) in enumerate(calls.items()):
        new_groups[key] = i
        for val in vals:
            new_groups[val] = i

    return new_groups


def consensus(destruct_calls, lumpy_calls, svaba_calls, gridss_calls, consensus_calls, sample_id, tempdir,
              region=None):
    temp_consensus_output = os.path.join(tempdir, 'consensus.csv')
    allcalls = [
        read_destruct(destruct_calls, region=region),
        SvVcfData(lumpy_calls, region=region).as_data_frame(),
        SvVcfData(svaba_calls, region=region).as_data_frame(),
        SvVcfData(gridss_calls, region=region).as_data_frame()
    ]

    allcalls = pd.concat(allcalls)

    allcalls_db = BreakpointDatabase(allcalls)

    groups = get_common_calls(allcalls, allcalls_db)

    allcalls['grouped_breakpoint_id'] = allcalls['breakpoint_id'].apply(lambda x: groups.get(x, float("nan")))

    allcalls = allcalls[~ pd.isnull(allcalls.grouped_breakpoint_id)]

    allcalls = allcalls.groupby('grouped_breakpoint_id')

    outdata = []
    for _, brkgrp in allcalls:

        # filter multiple calls by same tool in the window
        # without confirmation from another tool
        if len(brkgrp.caller.unique()) == 1:
            continue

        brkgrp['caller'] = ','.join(list(brkgrp['caller']))
        brkgrp = brkgrp[:1]

        outdata.append(brkgrp)

    if outdata:
        outdata = pd.concat(outdata)
    else:
        columns = [
            'breakpoint_id', 'caller', 'chromosome_1', 'chromosome_2', 'position_1',
            'position_2', 'strand_1', 'strand_2', 'type', 'grouped_breakpoint_id'
        ]
        outdata = pd.DataFrame(columns=columns)

    outdata['sample_id'] = sample_id

    outdata.to_csv(temp_consensus_output, index=False)

    csverve.rewrite_csv_file(temp_consensus_output, consensus_calls, dtypes=dtypes()['consensus'])
