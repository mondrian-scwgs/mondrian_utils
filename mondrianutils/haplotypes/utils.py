import os

import argparse
import csverve.api as csverve
import remixt.analysis
import remixt.workflow
from mondrianutils import helpers

import remixt


def annotate_haps(haps_csv, refdir, tempdir, output_csv):
    helpers.makedirs(tempdir)
    temp_output = os.path.join(tempdir, 'output.csv')

    thousand_genomes = os.path.join(refdir, 'thousand_genomes_snps.tsv')

    annotation_data = {}

    with helpers.getFileHandle(thousand_genomes, 'rt') as db:
        for line in db:
            line = line.strip().split('\t')

            chrom, pos, ref, alt = line

            annotation_data[(chrom, pos)] = (ref, alt)

    with helpers.getFileHandle(haps_csv, 'rt') as reader, helpers.getFileHandle(temp_output, 'wt') as writer:

        header = reader.readline().strip()
        header += '\tref\talt\n'
        writer.write(header)

        for line in reader:
            line = line.strip()
            l_split = line.split('\t')

            chrom = l_split[0]
            pos = l_split[1]

            if (chrom, pos) in annotation_data:
                ref, alt = annotation_data[(chrom, pos)]
            else:
                ref = 'NA'
                alt = 'NA'

            line += '\t{}\t{}\n'.format(ref, alt)

            writer.write(line)

    csverve.rewrite_csv_file(
        temp_output, output_csv,
        dtypes={
            'chromosome': 'str', 'position': 'int', 'alt': 'str', 'ref': 'str',
            'allele': 'str', 'hap_label': 'str', 'allele_id': 'str', }
    )

def infer_snp_genotype(seqdata, output, chromosome, all_chromosomes, ref_data_dir, ref_genome, is_tumour=False):
    config = {
        'chromosomes': all_chromosomes,
        'extract_seqdata': {
            'genome_fasta_template': ref_genome,
            'genome_fai_template': ref_genome + '.fai',
        },
        'ref_data_dir': ref_data_dir,
    }

    if is_tumour:
        remixt.analysis.haplotype.infer_snp_genotype_from_tumour(output, seqdata, chromosome, config)
    else:
        remixt.analysis.haplotype.infer_snp_genotype_from_normal(output, seqdata, chromosome, config)

def infer_haps(output, snp_genotype, chromosome, tempdir, ref_data_dir, ref_genome):
    config = {
        'genome_fasta_template': ref_genome,
        'genome_fai_template': ref_genome + '.fai',
    }

    remixt.analysis.haplotype.infer_haps(
        output, snp_genotype, chromosome, tempdir, config, ref_data_dir
    )

def merge_haps(inputs, output):
    remixt.utils.merge_tables(inputs, output)

def parse_args():
    default_chroms = [str(v) for v in range(1, 23)] + ['X', 'Y']

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    extract_seqdata = subparsers.add_parser('extract_seqdata')
    extract_seqdata.set_defaults(which='extract_seqdata')
    extract_seqdata.add_argument(
        '--bam',
        required=True
    )
    extract_seqdata.add_argument(
        '--output',
        required=True
    )
    extract_seqdata.add_argument(
        '--ref_data_dir',
        required=True
    )
    extract_seqdata.add_argument(
        '--chromosomes',
        nargs='*',
        default=default_chroms
    )

    infer_snp_genotype = subparsers.add_parser('infer_snp_genotype')
    infer_snp_genotype.set_defaults(which='infer_snp_genotype')
    infer_snp_genotype.add_argument(
        '--seqdata',
        required=True
    )
    infer_snp_genotype.add_argument(
        '--output',
        required=True
    )
    infer_snp_genotype.add_argument(
        '--chromosome',
        required=True
    )
    infer_snp_genotype.add_argument(
        '--chromosomes',
        nargs='*',
        default=default_chroms
    )
    infer_snp_genotype.add_argument(
        '--ref_data_dir',
        required=True
    )
    infer_snp_genotype.add_argument(
        '--ref_genome',
        required=True
    )
    infer_snp_genotype.add_argument(
        '--is_tumour',
        action='store_true',
        default=False
    )

    infer_haps = subparsers.add_parser('infer_haps')
    infer_haps.set_defaults(which='infer_haps')
    infer_haps.add_argument(
        '--snp_genotype',
        required=True
    )
    infer_haps.add_argument(
        '--output',
        required=True
    )
    infer_haps.add_argument(
        '--chromosome',
        required=True
    )
    infer_haps.add_argument(
        '--ref_data_dir',
        required=True
    )
    infer_haps.add_argument(
        '--ref_genome',
        required=True
    )

    merge_haps = subparsers.add_parser('merge_haps')
    merge_haps.set_defaults(which='merge_haps')
    merge_haps.add_argument(
        '--inputs',
        nargs='*',
        required=True
    )
    merge_haps.add_argument(
        '--output',
        required=True
    )

    annotate_haps = subparsers.add_parser('annotate_haps')
    annotate_haps.set_defaults(which='annotate_haps')
    annotate_haps.add_argument(
        '--input',
        required=True
    )
    annotate_haps.add_argument(
        '--output',
        required=True
    )
    annotate_haps.add_argument(
        '--ref_data_dir',
        required=True
    )
    annotate_haps.add_argument(
        '--tempdir',
        required=True
    )

    args = vars(parser.parse_args())

    return args

def utils():
    args = parse_args()

    if args['which'] == 'extract_seqdata':
        remixt.workflow.create_extract_seqdata_workflow(
            args['bam'], args['output'], args['ref_data_dir'], args['chromosomes']
        )
    elif args['which'] == 'infer_snp_genotype':
        infer_snp_genotype(
            args['seqdata'], args['output'], args['chromosome'], args['chromosomes'],
            args['ref_data_dir'], args['ref_genome'], args['is_tumour']
        )
    elif args['which'] == 'infer_haps':
        infer_haps(
            args['output'], args['snp_genotype'], args['chromosome'], args['tempdir'],
            args['ref_data_dir'], args['ref_genome']
        )
    elif args['which'] == 'merge_haps':
        merge_haps(
            args['inputs'], args['output']
        )
    elif args['which'] == 'annotate_haps':
        annotate_haps(
            args['inputs'], args['ref_data_dir'], args['tempdir'], args['output']
        )
    else:
        raise Exception()
