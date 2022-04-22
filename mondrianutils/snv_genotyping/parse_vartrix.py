import gzip

import csverve.api as csverve
import pandas as pd


def load_barcodes(barcode_file):
    data = {}
    with open(barcode_file, 'rt') as reader:
        for i, line in enumerate(reader):
            line = line.strip().split()
            assert len(line) == 1
            line = line[0]
            assert line not in data

            data[i + 1] = line
    return data


def load_matrix(matrix_file, cells, variants):
    data = {}

    with open(matrix_file, 'rt') as reader:
        for line in reader:
            if not line.startswith("%"):
                break

        for line in reader:
            line = line.strip().split()

            variant = variants[int(line[0])]
            cell = cells[int(line[1])]
            count = int(line[2])

            if cell not in data:
                data[cell] = {}

            data[cell][variant] = count

    return data


def get_vcf_pos(vcf_file):
    data = set()

    with gzip.open(vcf_file, 'rt') as reader:
        for line in reader:
            if line.startswith('#'):
                continue
            line = line.strip().split()
            chrom = line[0]
            pos = int(line[1])
            data.add((chrom, pos))
    return data


def create_df(ref, alt, cells, variants, sparse=False):
    out = []

    for i, cell in cells.items():
        for i, variant in variants.items():
            chrom, pos = variant.split('_')
            pos = int(pos)

            pos += 1

            ref_count = ref.get(cell, {}).get(variant, 0)
            alt_count = alt.get(cell, {}).get(variant, 0)

            # only write zeros if sparse is true
            if ref_count == 0 and alt_count == 0 and not sparse:
                continue

            out.append([cell, chrom, pos, ref_count, alt_count])

    df = pd.DataFrame(out, columns=['cell_id', 'chromosome', 'position', 'ref_count', 'alt_count'])
    return df


def parse_vartrix(cells, variants, ref_counts, alt_counts, outfile, write_header=True, sparse=False):
    cells = load_barcodes(cells)
    variants = load_barcodes(variants)

    ref = load_matrix(ref_counts, cells, variants)
    alt = load_matrix(alt_counts, cells, variants)

    df = create_df(ref, alt, cells, variants, sparse=sparse)

    csverve.write_dataframe_to_csv_and_yaml(
        df, outfile,
        write_header=write_header,
        dtypes={
            'cell_id': 'str',
            'chromosome': 'str',
            'position': 'int',
            'ref_count': 'int',
            'alt_count': 'int'
        }
    )
