import csv
import gzip
import os

import argparse
import pandas as pd
import pysam


def build_repeats(rmsk, chr_map, repeats, satellites):
    with open(chr_map, 'rt') as chr_map_file:
        chr_map = {}
        for line in chr_map_file:
            line = line.strip().split()
            chr_map[line[0]] = line[1]

    with open(rmsk, 'rt') as repeat_reader, open(repeats, 'wt') as repeat_writer, open(satellites,
                                                                                       'wt') as satellite_writer:
        for row in csv.reader(repeat_reader, delimiter='\t'):
            if row[5] not in chr_map:
                continue
            chrom = chr_map[row[5]]
            start = str(int(row[6]) + 1)
            end = row[7]
            record_type = row[11]
            repeat_writer.write('\t'.join([chrom, start, end]) + '\n')
            if record_type == 'Satellite':
                satellite_writer.write('\t'.join([chrom, start, end]) + '\n')


def get_chrom_lengths(fai):
    chroms = [
        '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
        '11', '12', '13', '14', '15', '16', '17', '18', '19',
        '20', '21', '22', 'X'
    ]
    chroms = set(['chr' + v for v in chroms])

    data = {}
    with open(fai, 'rt') as reader:
        for line in reader:
            line = line.strip().split()
            if line[0] in chroms:
                data[line[0]] = int(line[1])
    return data


def get_legend_map(legend_dir):
    files = os.listdir(legend_dir)
    files = [v for v in files if v.endswith('impute.legend.gz')]

    legend_map = {}
    for filename in files:
        assert filename.startswith('ALL_1000G_phase1integrated_v3_')
        chrom = filename[len('ALL_1000G_phase1integrated_v3_'):-len("_impute.legend.gz")]
        legend_map[chrom] = os.path.join(legend_dir, filename)

    return legend_map


def create_snp_positions_grch37(fai, data_dir, snp_positions):
    chroms = get_chrom_lengths(fai)
    legend_map = get_legend_map(data_dir)

    with open(snp_positions, 'wt') as writer:
        for chromosome in chroms:
            phased_chromosome = 'chrX_nonPAR' if chromosome == 'chrX' else chromosome
            legend_filename = legend_map[phased_chromosome]

            with gzip.open(legend_filename, 'rt') as legend_file:
                for line in legend_file:
                    if line.startswith('id'):
                        continue
                    row = line.split()
                    position = row[1]
                    a0 = row[2]
                    a1 = row[3]
                    if len(a0) != 1 or len(a1) != 1:
                        continue
                    writer.write('\t'.join([chromosome, position, a0, a1]) + '\n')


def create_snp_positions_grch38(bcf_file_dir, output, chromosomes=None, chr_name_prefix='chr'):
    if chromosomes is None:
        chromosomes = [chr_name_prefix + str(v) for v in range(1, 23)] + ['X']

    snps = []
    for chromosome in chromosomes:
        bcf_filename = os.path.join(
            bcf_file_dir,
            f'1kGP_high_coverage_Illumina.{chromosome}.filtered.SNV_INDEL_SV_phased_panel.bcf'
        )
        for r in pysam.VariantFile(bcf_filename, 'r'):
            for alt in r.alts:
                # Translate chromosome names from grch38 1kg (chr prefix) to those in the remixt reference
                # as specified by the chr_name_prefix config param
                if chr_name_prefix == 'chr':
                    chrom = r.chrom
                elif chr_name_prefix == '':
                    assert r.chrom.startswith('chr')
                    chrom = r.chrom[3:]
                else:
                    raise ValueError(f'unrecognized chr_name_prefix {chr_name_prefix}')
                coord = r.pos
                ref = r.ref
                if ref not in ['A', 'C', 'T', 'G']:
                    continue
                if alt not in ['A', 'C', 'T', 'G']:
                    continue
                snps.append([chrom, coord, ref, alt])
    snps = pd.DataFrame(snps, columns=['chrom', 'coord', 'ref', 'alt'])
    snps.to_csv(output, index=False, header=False, sep='\t')


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    repeats = subparsers.add_parser('repeats')
    repeats.set_defaults(which='repeats')
    repeats.add_argument('--rmsk_file', required=True)
    repeats.add_argument('--chrom_map', required=True)
    repeats.add_argument('--repeats', required=True)
    repeats.add_argument('--satellites', required=True)

    snp_positions_grch37 = subparsers.add_parser('snp_positions_grch37')
    snp_positions_grch37.set_defaults(which='snp_positions_grch37')
    snp_positions_grch37.add_argument(
        '--reference_fai'
    )
    snp_positions_grch37.add_argument(
        '--data_dir'
    )
    snp_positions_grch37.add_argument(
        '--snp_positions'
    )

    snp_positions_grch38 = subparsers.add_parser('snp_positions_grch38')
    snp_positions_grch38.set_defaults(which='snp_positions_grch38')
    snp_positions_grch38.add_argument(
        '--data_dir'
    )
    snp_positions_grch38.add_argument(
        '--snp_positions'
    )
    snp_positions_grch38.add_argument(
        '--chromosomes', nargs='*', default=None
    )

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'repeats':
        build_repeats(args['rmsk_file'], args['chrom_map'], args['repeats'], args['satellites'])

    elif args['which'] == 'snp_positions_grch37':
        create_snp_positions_grch37(
            args['reference_fai'], args['data_dir'], args['snp_positions']
        )
    elif args['which'] == 'snp_positions_grch38':
        create_snp_positions_grch38(
            args['data_dir'], args['snp_positions'], chromosomes=args['chromosomes']
        )
    else:
        raise Exception()
