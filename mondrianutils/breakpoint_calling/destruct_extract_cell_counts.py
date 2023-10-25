import pandas as pd
import csverve
from mondrianutils.dtypes.breakpoint import dtypes


def get_counts(breakpoint_read_filename, output_csv):
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

    bkp_read_table['comment'] = bkp_read_table['comment'].str.lstrip('+')
    bkp_read_table['read_id'] = bkp_read_table['comment'].str.split(':CB_').str[0]
    bkp_read_table['cell_id'] = bkp_read_table['comment'].str.split(':CB_').str[1]
    bkp_read_table = bkp_read_table.drop(columns=['comment']).drop_duplicates()

    bkp_read_table = bkp_read_table[bkp_read_table['filtered'] == False]
    bkp_read_table = bkp_read_table[['prediction_id', 'read_id', 'cell_id']].drop_duplicates()
    bkp_read_table = bkp_read_table.groupby(['prediction_id', 'cell_id']).size().rename(
        'read_count').reset_index()

    csverve.write_dataframe_to_csv_and_yaml(
        bkp_read_table,
        output_csv,
        dtypes()['genotyping'],
        write_header=True
    )
