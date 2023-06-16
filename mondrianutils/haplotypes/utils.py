import argparse
import csverve.api as csverve
import json
import os
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
from mondrianutils.haplotypes import shapeit

import remixt


def add_cell_id_to_seqdata(seqdata, cellid):
    with pd.HDFStore(seqdata) as store:
        store.put('cell_id', pd.Series([cellid]))


def get_cell_id_from_seqdata(seqdata):
    with pd.HDFStore(seqdata) as store:
        if '/cell_id' not in store.keys():
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


def finalize_tsv(infile, outfile, seqdata, skip_header=False):
    df = pd.read_csv(infile, sep='\t')

    cellid = get_cell_id_from_seqdata(seqdata)
    if cellid:
        df['cell_id'] = cellid

    csverve.write_dataframe_to_csv_and_yaml(
        df, outfile, skip_header=skip_header, dtypes=dtypes()
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
        dtypes=dtypes()
    )


def convert_csv_to_tsv(csv_infile, tsv_outfile):
    df = csverve.read_csv(csv_infile)
    df.to_csv(tsv_outfile, sep='\t', index=False)


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
        '--chromosomes',
        nargs='*',
        default=default_chroms
    )
    create_segments.add_argument(
        '--output',
        required=True
    )
    create_segments.add_argument(
        '--tempdir',
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
    haplotype_allele_readcount.add_argument(
        '--skip_header',
        action='store_true',
        default=False
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

    run_shapeit = subparsers.add_parser('run_shapeit')
    run_shapeit.set_defaults(which='run_shapeit')
    run_shapeit.add_argument(
        '--input_bcf_file', required=True
    )
    run_shapeit.add_argument(
        '--genetic_map', required=True
    )
    run_shapeit.add_argument(
        '--regions_file', required=True
    )
    run_shapeit.add_argument(
        '--chromosome', required=True
    )
    run_shapeit.add_argument(
        '--tempdir', required=True
    )
    run_shapeit.add_argument(
        '--output', required=True
    )
    run_shapeit.add_argument(
        '--phased_chromosomes', nargs='*',
        default=tuple([f'chr{v}' for v in range(1, 23)] + ['chrX'])
    )
    run_shapeit.add_argument(
        '--is_female', default=False, action='store_true'
    )
    run_shapeit.add_argument(
        '--phased_chromosome_x', default='chrX'
    )
    run_shapeit.add_argument(
        '--shapeit_num_samples', default=100, type=int
    )
    run_shapeit.add_argument(
        '--shapeit_confidence_threshold', default=0.95, type=float
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
    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'run_shapeit':
        shapeit.run_shapeit4(
            args['input_bcf_file'], args['genetic_map'], args['regions_file'], args['chromosome'],
            args['tempdir'], args['output'],
            phased_chromosomes=args['phased_chromosomes'],
            is_female=args['is_female'],
            phased_chromx=args['phased_chromosome_x'],
            shapeit_num_samples=args['shapeit_num_samples'],
            shapeit_confidence_threshold=args['shapeit_confidence_threshold']
        )
    elif args['which'] == 'extract_seqdata':
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
    elif args['which'] == 'create_segments':
        chr_prefix = 'chr' if args['chromosomes'][0].startswith('chr') else ''
        helpers.makedirs(args['tempdir'])
        new_gap = os.path.join(args['tempdir'], 'gap.txt.gz')
        with helpers.getFileHandle(args['gap_table'], 'rt') as reader, helpers.getFileHandle(new_gap, 'wt') as writer:
            for line in reader:
                line = line.strip().split('\t')
                if chr_prefix == 'chr':
                    line[1] = line[1] if line[1].startswith('chr') else 'chr' + line[1]
                else:
                    line[1] = line[1][3:] if line[1].startswith('chr') else line[1]
                line = '\t'.join(line) + '\n'
                writer.write(line)

        config = {
            'genome_fai_template': args['reference_fai'],
            'gap_table_template': new_gap,
            'chromosomes': args['chromosomes'],
            'chr_name_prefix': chr_prefix
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
        finalize_tsv(tempout, args['output'], args['seqdata'], skip_header=args['skip_header'])
    elif args['which'] == 'generate_metadata':
        generate_metadata(
            args['files'], args['metadata_yaml_files'], args['samples'],
            args['metadata_output']
        )
    elif args['which'] == "convert_haplotypes_csv_to_tsv":
        convert_csv_to_tsv(
            args['input'], args['output']
        )
    else:
        raise Exception()
