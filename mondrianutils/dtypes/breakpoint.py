def dtypes():
    consensus = {
        'chromosome_1': 'category',
        'position_1': 'int',
        'strand_1': 'str',
        'chromosome_2': 'category',
        'position_2': 'int',
        'strand_2': 'str',
        'type': 'category',
        'breakpoint_id': 'category',
        'caller': 'category',
        'grouped_breakpoint_id': 'category',
        'sample_id': 'category'
    }

    genotyping = {
        'prediction_id': int,
        'chromosome_1': 'category',
        'strand_1': str,
        'position_1': int,
        'chromosome_2': 'category',
        'strand_2': str,
        'position_2': int,
        'cell_id': 'category',
        'read_count': int

    }

    dtypes = locals()

    return dtypes
