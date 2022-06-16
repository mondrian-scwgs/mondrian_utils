import gzip
import os

import csverve.api as csverve
import pandas as pd


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

    with open(count_file, 'rt') as reader:

        for line in reader:
            if line.startswith('%'):
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
        outstr = [str(barcodes[barcode]), str(variants[variant]), count]
        outstr = '\t'.join(outstr) + '\n'
        writer.write(outstr)


def write_parsed_format(ref_counts, alt_counts, parsed_output):
    for (variant, barcode, ref, alt, idx), ref_count in ref_counts.items():
        alt_count = alt_counts[(variant, barcode, ref, alt, idx)]

        del alt_counts[(variant, barcode, ref, alt, idx)]

        chrom, pos = variant.split('_')

        outstr = [barcode, chrom, pos, ref, alt, ref_count, alt_count]

        outstr = ','.join(outstr) + '\n'

        parsed_output.write(outstr)


def rewrite_counts_file(inputfile, outputfile, num_barcodes, num_variants):
    with open(inputfile, 'rt') as reader:
        for num_calls, _ in enumerate(reader):
            pass

    with open(inputfile, 'rt') as reader, open(outputfile, 'wt') as writer:

        writer.write("%%MatrixMarket matrix coordinate real general\n")
        writer.write("% written by sprs\n")
        writer.write("{}\t{}\t{}\n".format(num_barcodes, num_variants, num_calls))

        for line in reader:
            writer.write(line)


def merge_vartrix(
        barcodes, variants, refs, alts, vcf_files,
        merged_barcodes, merged_variants, merged_ref,
        merged_alt, parsed_output, tempdir
):
    assert len(barcodes) == len(variants) == len(refs) == len(alts)

    try:
        os.makedirs(tempdir)
    except:
        pass

    all_barcodes = merge_idx_files(barcodes)
    all_variants = merge_idx_files(variants)

    write_idx_file(all_barcodes, merged_barcodes)
    write_idx_file(all_variants, merged_variants)

    barcodes_dict = {v: i + 1 for i, v in enumerate(all_barcodes)}
    variants_dict = {v: i + 1 for i, v in enumerate(all_variants)}

    temp_merged_ref = os.path.join(tempdir, 'merged_ref.mtx')
    temp_merged_alt = os.path.join(tempdir, 'merged_alt.mtx')
    temp_parsed = os.path.join(tempdir, 'parsed.csv')

    with open(temp_merged_ref, 'wt') as ref_writer, open(temp_merged_alt, 'wt') as alt_writer, open(temp_parsed,
                                                                                                    'wt') as parsed_writer:

        parsed_writer.write('cell_id,chromosome,position,ref,alt,ref_count,alt_count\n')

        for barcode, variant, ref, alt, vcf in zip(barcodes, variants, refs, alts, vcf_files):
            assert os.path.dirname(barcode) == os.path.dirname(variant) == os.path.dirname(ref) == os.path.dirname(alt)

            vcf_data = load_vcf(vcf)
            barcode_data = load_idx_file(barcode)
            variant_data = load_idx_file(variant)

            ref_data = load_data(barcode_data, variant_data, vcf_data, ref)
            write_counts_file(barcodes_dict, variants_dict, ref_data, ref_writer)

            alt_data = load_data(barcode_data, variant_data, vcf_data, alt)
            write_counts_file(barcodes_dict, variants_dict, alt_data, alt_writer)

            write_parsed_format(ref_data, alt_data, parsed_writer)

    rewrite_counts_file(temp_merged_ref, merged_ref, len(all_barcodes), len(all_variants))
    rewrite_counts_file(temp_merged_alt, merged_alt, len(all_barcodes), len(all_variants))

    df = pd.read_csv(temp_parsed)
    csverve.write_dataframe_to_csv_and_yaml(
        df, parsed_output,
        write_header=True,
        dtypes={
            'cell_id': 'str',
            'chromosome': 'str',
            'position': 'int',
            'ref_count': 'int',
            'alt_count': 'int'
        }
    )
