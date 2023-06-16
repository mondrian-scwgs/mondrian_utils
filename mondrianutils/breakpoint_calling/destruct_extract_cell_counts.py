import pysam
import pandas as pd
import csverve
from mondrianutils import helpers
from mondrianutils.dtypes.breakpoint import dtypes


def build_read_id_to_cb_table(sam_filename, reads_of_interest, region):
    read_id_to_cb = {}

    chrom, start, end = helpers.parse_region(region)

    with pysam.AlignmentFile(sam_filename, 'r', check_sq=False) as samfile:

        if start is not None:
            fetcher = samfile.fetch(chrom, start, end)
        elif chrom is not None:
            fetcher = samfile.fetch(chrom)
        else:
            fetcher = samfile.fetch()

        for read in fetcher:
            read_id = read.query_name

            if read_id not in reads_of_interest:
                continue

            cb_tag = read.get_tag('CB')

            if read_id and cb_tag:
                read_id_to_cb[read_id] = cb_tag

    return read_id_to_cb


def get_counts(breakpoint_read_filename, bam_filename, output_csv, region=None):
    destruct_reads_cols = [
        'prediction_id',
        'library_id',
        'fragment_id',
        'read_end',
        'seq',
        'qual',
        'comment',
        'filtered',
    ]
    bkp_read_table = pd.read_csv(
        breakpoint_read_filename,
        header=None, sep='\t',
        names=destruct_reads_cols,
        usecols=['prediction_id', 'comment', 'filtered'])

    bkp_read_table['read_id'] = bkp_read_table['comment'].str.lstrip('+')
    bkp_read_table = bkp_read_table.drop(columns=['comment']).drop_duplicates()

    reads = set(bkp_read_table['read_id'].to_list())

    read_id_to_cb = build_read_id_to_cb_table(bam_filename, reads, region=region)

    read_id_to_cb = pd.Series(read_id_to_cb).rename('cell_id').rename_axis('read_id').reset_index()

    breakpoint_cell_reads = bkp_read_table.merge(read_id_to_cb, on='read_id', how='left')

    breakpoint_cell_counts = breakpoint_cell_reads[breakpoint_cell_reads['filtered'] == False]
    breakpoint_cell_counts = breakpoint_cell_counts[['prediction_id', 'read_id', 'cell_id']].drop_duplicates()
    breakpoint_cell_counts = breakpoint_cell_counts.groupby(['prediction_id', 'cell_id']).size().rename(
        'read_count').reset_index()

    csverve.write_dataframe_to_csv_and_yaml(
        breakpoint_cell_counts,
        output_csv,
        dtypes()['genotyping'],
        write_header=True
    )


def merge_counts_files(infiles, outfile):
    on = ['prediction_id', 'cell_id']
    output_df = csverve.read_csv(infiles[0])

    for filepath in infiles[1:]:
        new_df = pd.read_csv(filepath)
        output_df = pd.concat([output_df, new_df]).groupby(on).sum().reset_index()

    csverve.write_dataframe_to_csv_and_yaml(
        output_df, outfile, dtypes()['genotyping'], write_header=True
    )
