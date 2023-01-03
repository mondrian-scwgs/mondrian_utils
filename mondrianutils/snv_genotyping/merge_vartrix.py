import gzip
import os

import csverve.api as csverve
import mondrianutils.helpers as helpers
import pandas as pd
from mondrianutils.dtypes.variant import dtypes


def merge_idx_files(idx_files):
    idxs = set()
    for idx_file in idx_files:
        with open(idx_file, 'rt') as reader:
            for line in reader:
                idxs.add(line.strip())
    return sorted(idxs)


def write_idx_file(idx_data, idx_file):
    with open(idx_file, 'wt') as writer:
        for idx_val in idx_data:
            writer.write(idx_val + '\n')


def load_idx_file(idx_file):
    with open(idx_file, 'rt') as reader:
        idx_data = {i + 1: line.strip() for i, line in enumerate(reader)}
    return idx_data


def load_vcf(vcf_file):
    outdata = {}
    i = 1
    with gzip.open(vcf_file, 'rt') as reader:
        for line in reader:
            if line.startswith('#'):
                continue
            line = line.strip().split()
            chrom = line[0]
            pos = line[1]
            ref = line[3]
            alt = line[4]

            if (chrom, pos) in outdata:
                outdata[(chrom, pos)].append((i, ref, alt))
            else:
                outdata[(chrom, pos)] = [(i, ref, alt)]

            i += 1

    return outdata


def get_ref_alt(chrom, pos, idx, vcf_data):
    pos = str(int(pos) + 1)

    values = vcf_data[(chrom, pos)]

    values = [v for v in values if v[0] == idx]

    assert len(values) == 1

    return values[0][1], values[0][2]


def load_data(barcode_data, variant_data, vcf_data, count_file):
    counts = {}
    header = False

    with open(count_file, 'rt') as reader:

        for line in reader:
            if line.startswith('%'):
                continue
            # first line is the total counts
            if header is False:
                header = True
                continue
            var_idx, barcode_idx, count = line.strip().split()
            var_idx = int(var_idx)
            barcode_idx = int(barcode_idx)

            chrom, pos = variant_data[var_idx].split('_')

            ref, alt = get_ref_alt(chrom, pos, var_idx, vcf_data)

            key = (variant_data[var_idx], barcode_data[barcode_idx], ref, alt, var_idx)

            assert key not in counts, key

            counts[key] = count

    return counts


def write_counts_file(barcodes, variants, data, writer):
    for (variant, barcode, _, _, _), count in data.items():
        outstr = [str(variants[variant]), str(barcodes[barcode]), count]
        outstr = '\t'.join(outstr) + '\n'
        writer.write(outstr)


def write_parsed_format(ref_counts, alt_counts, parsed_output):
    for (variant, barcode, ref, alt, idx), ref_count in ref_counts.items():
        alt_count = alt_counts[(variant, barcode, ref, alt, idx)]

        del alt_counts[(variant, barcode, ref, alt, idx)]

        chrom, pos = variant.split('_')

        pos = str(int(pos) + 1)

        outstr = [barcode, chrom, pos, ref, alt, ref_count, alt_count]

        outstr = ','.join(outstr) + '\n'

        parsed_output.write(outstr)


def rewrite_counts_file(inputfile, outputfile, num_barcodes, num_variants):
    num_calls = 0
    with open(inputfile, 'rt') as reader:
        for num_calls, _ in enumerate(reader):
            pass

    with open(inputfile, 'rt') as reader, open(outputfile, 'wt') as writer:
        writer.write("%%MatrixMarket matrix coordinate real general\n")
        writer.write("% written by sprs\n")
        writer.write("{}\t{}\t{}\n".format(num_variants, num_barcodes, num_calls))
        for line in reader:
            writer.write(line)


def merge_vartrix(
        barcodes, variants, refs, alts, vcf_files,
        parsed_output, tempdir
):
    assert len(barcodes) == len(variants) == len(refs) == len(alts)
    helpers.makedirs(tempdir)

    temp_parsed = os.path.join(tempdir, 'parsed.csv')

    with open(temp_parsed, 'wt') as parsed_writer:
        parsed_writer.write('cell_id,chromosome,position,ref,alt,ref_count,alt_count\n')

        for barcode, variant, ref, alt, vcf in zip(barcodes, variants, refs, alts, vcf_files):
            assert os.path.dirname(barcode) == os.path.dirname(variant) == os.path.dirname(ref) == os.path.dirname(alt)

            vcf_data = load_vcf(vcf)
            barcode_data = load_idx_file(barcode)
            variant_data = load_idx_file(variant)

            ref_data = load_data(barcode_data, variant_data, vcf_data, ref)

            alt_data = load_data(barcode_data, variant_data, vcf_data, alt)

            write_parsed_format(ref_data, alt_data, parsed_writer)

    df = pd.read_csv(temp_parsed)
    csverve.write_dataframe_to_csv_and_yaml(
        df, parsed_output,
        skip_header=False,
        dtypes=dtypes()['genotyping']
    )


def parse_vartrix(barcode, variant, ref_counts, alt_counts, vcf_file, parsed_output, tempdir, skip_header=False):
    helpers.makedirs(tempdir)

    temp_parsed = os.path.join(tempdir, 'parsed.csv')

    with open(temp_parsed, 'wt') as parsed_writer:
        parsed_writer.write('cell_id,chromosome,position,ref,alt,ref_count,alt_count\n')

        vcf_data = load_vcf(vcf_file)
        barcode_data = load_idx_file(barcode)
        variant_data = load_idx_file(variant)

        ref_data = load_data(barcode_data, variant_data, vcf_data, ref_counts)

        alt_data = load_data(barcode_data, variant_data, vcf_data, alt_counts)

        write_parsed_format(ref_data, alt_data, parsed_writer)

    df = pd.read_csv(temp_parsed)
    csverve.write_dataframe_to_csv_and_yaml(
        df, parsed_output,
        skip_header=skip_header,
        dtypes=dtypes()['genotyping']
    )


def regenerate_vartrix_format(barcodes_file, variants_file, ref_matrix, alt_matrix, parsed_data, tempdir):
    helpers.makedirs(tempdir)

    barcodes = set()
    variants = set()

    with helpers.getFileHandle(parsed_data, 'rt') as reader:
        header = reader.readline().strip().split(',')
        header = {v: i for i, v in enumerate(header)}

        for line in reader:
            line = line.strip().split(',')
            cell_id = line[header['cell_id']]
            chrom = line[header['chromosome']]
            pos = int(line[header['position']])

            barcodes.add(cell_id)
            variants.add(f'{chrom}_{pos}')

    barcodes = sorted(barcodes)
    variants = sorted(variants)

    len_barcodes = len(barcodes)
    len_variants = len(variants)

    write_idx_file(barcodes, barcodes_file)
    write_idx_file(variants, variants_file)

    barcodes = {barcode: i + 1 for i, barcode in enumerate(barcodes)}
    variants = {variant: i + 1 for i, variant in enumerate(variants)}

    temp_merged_ref = os.path.join(tempdir, 'merged_ref.mtx')
    temp_merged_alt = os.path.join(tempdir, 'merged_alt.mtx')

    with open(temp_merged_ref, 'wt') as ref_writer, open(temp_merged_alt, 'wt') as alt_writer, helpers.getFileHandle(parsed_data, 'rt') as reader:
        header = reader.readline().strip().split(',')
        header = {v: i for i, v in enumerate(header)}
        for line in reader:
            line = line.strip().split(',')
            cell_id = line[header['cell_id']]
            chrom = line[header['chromosome']]
            pos = int(line[header['position']])
            ref_count = int(line[header['ref_count']])
            alt_count = int(line[header['alt_count']])

            outstr = [str(variants[f'{chrom}_{pos}']), str(barcodes[cell_id]), str(ref_count)]
            outstr = '\t'.join(outstr) + '\n'
            ref_writer.write(outstr)

            outstr = [str(variants[f'{chrom}_{pos}']), str(barcodes[cell_id]), str(alt_count)]
            outstr = '\t'.join(outstr) + '\n'
            alt_writer.write(outstr)

    rewrite_counts_file(temp_merged_ref, ref_matrix, len_barcodes, len_variants)
    rewrite_counts_file(temp_merged_alt, alt_matrix, len_barcodes, len_variants)
