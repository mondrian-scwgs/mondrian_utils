import json
import os

import argparse
import csverve.api as csverve
import pandas as pd
import remixt.analysis
import remixt.analysis.haplotype
import remixt.analysis.readcount
import remixt.analysis.segment
import remixt.seqdataio
import remixt.workflow
import yaml
from mondrianutils import helpers
from mondrianutils.dtypes.haplotypes import dtypes

import remixt


def add_cell_id_to_seqdata(seqdata, cellid):
    with pd.HDFStore(seqdata) as store:
        store.put('cell_id', pd.Series([cellid]))


def get_cell_id_from_seqdata(seqdata):
    with pd.HDFStore(seqdata) as store:
        if 'cell_id' not in store.keys():
            return
        cellid = store.get('cell_id')
        cellid = list(cellid)
        assert len(cellid) == 1
        return cellid[0]


def infer_type(files):
    with open(files, 'rt') as reader:
        files = json.load(reader)

    filetypes = sorted(set([v['left'] for v in files]))

    # more than one wf
    if 'haplotype_counts' in filetypes and 'infer_haplotype' in filetypes:
        return 'haplotype_calling'
    elif 'haplotype_counts' in filetypes:
        return 'haplotype_counting'
    elif 'infer_haplotype' in filetypes:
        return 'infer_haplotype'
    else:
        raise Exception()


def generate_metadata(
        files, metadata_yaml_files, samples, metadata_output
):
    wf_type = infer_type(files)
    data = helpers.metadata_helper(files, metadata_yaml_files, samples, wf_type)

    with open(metadata_output, 'wt') as writer:
        yaml.dump(data, writer, default_flow_style=False)


def convert_csv_to_tsv(csv_infile, tsv_outfile):
    df = csverve.read_csv(csv_infile)
    df.to_csv(tsv_outfile, sep='\t', index=False)


def finalize_tsv(infile, outfile, seqdata):
    df = pd.read_csv(infile, sep='\t')

    cellid = get_cell_id_from_seqdata(seqdata)
    if cellid:
        df['cell_id'] = cellid

    csverve.write_dataframe_to_csv_and_yaml(
        df, outfile, write_header=True, dtypes=dtypes()
    )


def annotate_haps(haps_csv, thousand_genomes, tempdir, output_csv):
    helpers.makedirs(tempdir)
    temp_output = os.path.join(tempdir, 'output.csv')

    annotation_data = {}

    with helpers.getFileHandle(thousand_genomes, 'rt') as db:
        for line in db:
            line = line.strip().split('\t')

            chrom, pos, ref, alt = line

            annotation_data[(chrom, pos)] = (ref, alt)

    with helpers.getFileHandle(haps_csv, 'rt') as reader, helpers.getFileHandle(temp_output, 'wt') as writer:

        header = reader.readline().strip().split('\t')
        header.extend(['ref', 'alt'])
        header = ','.join(header) + '\n'
        writer.write(header)

        for line in reader:
            line = line.strip().split('\t')

            chrom = line[0]
            pos = line[1]

            if (chrom, pos) in annotation_data:
                ref, alt = annotation_data[(chrom, pos)]
            else:
                ref = 'NA'
                alt = 'NA'

            line.extend([ref, alt])
            line = ','.join(line) + '\n'

            writer.write(line)

    csverve.rewrite_csv_file(
        temp_output, output_csv,
        dtypes={
            'chromosome': 'str', 'position': 'int', 'alt': 'str', 'ref': 'str',
            'allele': 'str', 'hap_label': 'str', 'allele_id': 'str', }
    )


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
        '--snp_positions',
        required=True
    )
    extract_seqdata.add_argument(
        '--tempdir',
        required=True
    )
    extract_seqdata.add_argument(
        '--chromosomes',
        nargs='*',
        default=default_chroms
    )
    extract_seqdata.add_argument(
        '--cell_id'
    )

    extract_chromosome_seqdata = subparsers.add_parser('extract_chromosome_seqdata')
    extract_chromosome_seqdata.set_defaults(which='extract_chromosome_seqdata')
    extract_chromosome_seqdata.add_argument(
        '--bam',
        required=True
    )
    extract_chromosome_seqdata.add_argument(
        '--output',
        required=True
    )
    extract_chromosome_seqdata.add_argument(
        '--snp_positions',
        required=True
    )
    extract_chromosome_seqdata.add_argument(
        '--chromosome',
        required=True
    )

    infer_snp_genotype = subparsers.add_parser('infer_snp_genotype_from_normal')
    infer_snp_genotype.set_defaults(which='infer_snp_genotype_from_normal')
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

    infer_haps = subparsers.add_parser('infer_haps')
    infer_haps.set_defaults(which='infer_haps')
    infer_haps.add_argument(
        '--snp_genotype',
        required=True
    )
    infer_haps.add_argument(
        '--thousand_genomes_impute_tar',
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
        '--tempdir',
        required=True
    )
    infer_haps.add_argument(
        '--genetic_map_filename_template',
        required=True
    )
    infer_haps.add_argument(
        '--haplotypes_filename_template',
        required=True
    )
    infer_haps.add_argument(
        '--legend_filename_template',
        required=True
    )
    infer_haps.add_argument(
        '--sample_filename',
        required=True
    )
    infer_haps.add_argument(
        '--phased_chromosome_x',
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

    merge_seqdata = subparsers.add_parser('merge_seqdata')
    merge_seqdata.set_defaults(which='merge_seqdata')
    merge_seqdata.add_argument(
        '--inputs',
        nargs='*',
        required=True
    )
    merge_seqdata.add_argument(
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
        '--thousand_genomes',
        required=True
    )
    annotate_haps.add_argument(
        '--output',
        required=True
    )
    annotate_haps.add_argument(
        '--tempdir',
        required=True
    )

    create_segments = subparsers.add_parser('create_segments')
    create_segments.set_defaults(which='create_segments')
    create_segments.add_argument(
        '--reference_fai',
        required=True
    )
    create_segments.add_argument(
        '--gap_table',
        required=True
    )
    create_segments.add_argument(
        '--output',
        required=True
    )

    convert_haplotypes_csv_to_tsv = subparsers.add_parser('convert_haplotypes_csv_to_tsv')
    convert_haplotypes_csv_to_tsv.set_defaults(which='convert_haplotypes_csv_to_tsv')
    convert_haplotypes_csv_to_tsv.add_argument(
        '--input',
        required=True
    )
    convert_haplotypes_csv_to_tsv.add_argument(
        '--output',
        required=True
    )

    haplotype_allele_readcount = subparsers.add_parser('haplotype_allele_readcount')
    haplotype_allele_readcount.set_defaults(which='haplotype_allele_readcount')
    haplotype_allele_readcount.add_argument(
        '--segments',
        required=True
    )
    haplotype_allele_readcount.add_argument(
        '--seqdata',
        required=True
    )
    haplotype_allele_readcount.add_argument(
        '--haplotypes',
        required=True
    )
    haplotype_allele_readcount.add_argument(
        '--output',
        required=True
    )
    haplotype_allele_readcount.add_argument(
        '--tempdir',
        required=True
    )

    generate_metadata = subparsers.add_parser('generate_metadata')
    generate_metadata.set_defaults(which='generate_metadata')
    generate_metadata.add_argument(
        '--files'
    )
    generate_metadata.add_argument(
        '--metadata_yaml_files', nargs='*'
    )
    generate_metadata.add_argument(
        '--samples', nargs='*'
    )
    generate_metadata.add_argument(
        '--metadata_output'
    )

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'extract_seqdata':
        bam_max_fragment_length = remixt.config.get_param({}, 'bam_max_fragment_length')
        bam_max_soft_clipped = remixt.config.get_param({}, 'bam_max_soft_clipped')
        bam_check_proper_pair = remixt.config.get_param({}, 'bam_check_proper_pair')
        remixt.seqdataio.create_seqdata(
            args['output'], args['bam'], args['snp_positions'],
            bam_max_fragment_length, bam_max_soft_clipped,
            bam_check_proper_pair, args['tempdir'], args['chromosomes']
        )
        if args['cell_id']:
            add_cell_id_to_seqdata(args['output'], args['cell_id'])

    elif args['which'] == 'extract_chromosome_seqdata':
        bam_max_fragment_length = remixt.config.get_param({}, 'bam_max_fragment_length')
        bam_max_soft_clipped = remixt.config.get_param({}, 'bam_max_soft_clipped')
        bam_check_proper_pair = remixt.config.get_param({}, 'bam_check_proper_pair')
        remixt.seqdataio.create_chromosome_seqdata(
            args['output'], args['bam'], args['snp_positions'],
            args['chromosome'], bam_max_fragment_length,
            bam_max_soft_clipped, bam_check_proper_pair
        )
    elif args['which'] == 'infer_snp_genotype_from_normal':
        remixt.analysis.haplotype.infer_snp_genotype_from_normal(
            args['output'], args['seqdata'], args['chromosome'], {}
        )
    elif args['which'] == 'infer_haps':
        helpers.untar(
            args['thousand_genomes_impute_tar'],
            os.path.join(args['tempdir'], 'thousand_genomes_impute_tar')
        )

        chrom = args['phased_chromosome_x'] if args['chromosome'].endswith('X') else args['chromosome']

        assert '{chromosome}' in args['genetic_map_filename_template']
        genetic_map = args['genetic_map_filename_template'].replace('{chromosome}', chrom)
        genetic_map = os.path.join(args['tempdir'], 'thousand_genomes_impute_tar', genetic_map)

        assert '{chromosome}' in args['haplotypes_filename_template']
        haplotypes = args['haplotypes_filename_template'].replace('{chromosome}', chrom)
        haplotypes = os.path.join(args['tempdir'], 'thousand_genomes_impute_tar', haplotypes)

        assert '{chromosome}' in args['legend_filename_template']
        legend = args['legend_filename_template'].replace('{chromosome}', chrom)
        legend = os.path.join(args['tempdir'], 'thousand_genomes_impute_tar', legend)

        sample = os.path.join(args['tempdir'], 'thousand_genomes_impute_tar', args['sample_filename'])

        config = {
            'genetic_map_template': genetic_map,
            'haplotypes_template': haplotypes,
            'legend_template': legend,
            'sample_template': sample
        }
        remixt.analysis.haplotype.infer_haps(
            args['output'], args['snp_genotype'], args['chromosome'], args['tempdir'],
            config, None
        )
    elif args['which'] == 'merge_haps':
        remixt.utils.merge_tables(
            args['output'], *args['inputs']
        )
    elif args['which'] == 'merge_seqdata':
        remixt.seqdataio.merge_seqdata(
            args['output'], args['inputs']
        )
    elif args['which'] == 'annotate_haps':
        annotate_haps(
            args['input'], args['thousand_genomes'], args['tempdir'], args['output']
        )
    elif args['which'] == 'create_segments':
        config = {
            'genome_fai_template': args['reference_fai'],
            'gap_table_template': args['gap_table']
        }
        ref_data_dir = os.path.dirname(args['reference_fai'])
        remixt.analysis.segment.create_segments(
            args['output'], config, ref_data_dir
        )
    elif args['which'] == 'haplotype_allele_readcount':
        helpers.makedirs(args['tempdir'])
        tempout = os.path.join(args['tempdir'], 'temp.tsv')
        remixt.analysis.readcount.haplotype_allele_readcount(
            tempout, args['segments'], args['seqdata'], args['haplotypes'], {}
        )
        finalize_tsv(tempout, args['output'], args['seqdata'])
    elif args['which'] == "convert_haplotypes_csv_to_tsv":
        convert_csv_to_tsv(
            args['input'], args['output']
        )
    elif args['which'] == 'generate_metadata':
        generate_metadata(
            args['files'], args['metadata_yaml_files'], args['samples'],
            args['metadata_output']
        )
    else:
        raise Exception()
