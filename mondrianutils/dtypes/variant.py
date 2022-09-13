def dtypes():
    genotyping = {
        'cell_id': 'category',
        'chromosome': 'category',
        'position': 'int',
        'ref_count': 'int',
        'alt_count': 'int',
        'ref': 'category',
        'alt': 'category',
    }

    dtypes = locals()

    return dtypes
