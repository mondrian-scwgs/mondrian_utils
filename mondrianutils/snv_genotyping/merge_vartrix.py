import os

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


def load_data(barcode_file, variant_file, count_file):
    with open(barcode_file, 'rt') as reader:
        barcode_data = {i+1: barcode.strip() for i, barcode in enumerate(reader)}

    with open(variant_file, 'rt') as reader:
        variant_data = {i+1: variant.strip() for i, variant in enumerate(reader)}

    counts = {}
    with open(count_file, 'rt') as reader:

        for line in reader:
            if line.startswith('%'):
                continue
            var_idx, barcode_idx, count = line.strip().split()
            var_idx = int(var_idx)
            barcode_idx = int(barcode_idx)
            key = (variant_data[var_idx], barcode_data[barcode_idx])
            counts[key] = count
    return counts


def write_counts_file(barcodes, variants, data, writer):
    for (variant, barcode), count in data.items():
        outstr = [barcodes[barcode], variants[variant], count]
        outstr = '\t'.join(outstr) + '\n'
        writer.write(outstr)


def merge_vartrix(barcodes, variants, refs, alts, vcf_files, merged_barcodes, merged_variants, merged_ref, merged_alt):
    assert len(barcodes) == len(variants) == len(refs) == len(alts)

    all_barcodes = merge_idx_files(barcodes)
    all_variants = merge_idx_files(variants)

    write_idx_file(all_barcodes, merged_barcodes)
    write_idx_file(all_variants, merged_variants)

    barcodes_dict = {v: i+1 for i, v in enumerate(barcodes)}
    variants_dict = {v: i+1 for i, v in enumerate(variants)}

    with open(merged_ref, 'wt') as ref_writer, open(merged_alt, 'wt') as alt_writer:
        for barcode, variant, ref, alt, vcf in zip(barcodes, variants, refs, alts, vcf_files):
            assert os.path.dirname(barcode) == os.path.dirname(variant) == os.path.dirname(ref) == os.path.dirname(alt)
            ref_data = load_data(barcode, variant, ref)
            write_counts_file(barcodes_dict, variants_dict, ref_data, ref_writer)

            alt_data = load_data(barcode, variant, alt)
            write_counts_file(barcodes_dict, variants_dict, alt_data, alt_writer)
