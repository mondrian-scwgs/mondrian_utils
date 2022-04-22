import os

import argparse
import pysam
import yaml
from mondrianutils import __version__
from mondrianutils.snv_genotyping.parse_vartrix import parse_vartrix
from mondrianutils.snv_genotyping.snv_genotyper import SnvGenotyper


def generate_metadata(
        outputs, vartrix_outputs, metadata_input, metadata_output
):
    data = dict()
    data['files'] = {
        os.path.basename(outputs[0]): {'result_type': 'pysam_genotyping_counts', 'auxiliary': False},
        os.path.basename(outputs[1]): {'result_type': 'pysam_genotyping_counts', 'auxiliary': True},
        os.path.basename(vartrix_outputs[0]): {'result_type': 'vartrix_genotyping_counts', 'auxiliary': False},
        os.path.basename(vartrix_outputs[1]): {'result_type': 'vartrix_genotyping_counts', 'auxiliary': True},
    }

    with open(metadata_input, 'rt') as reader:
        meta = yaml.safe_load(reader)
        del meta['meta']['type']
        del meta['meta']['version']

    meta_dict = {
        'type': 'snv_genotyping',
        'version': __version__
    }

    data['meta'] = {**meta_dict, **meta['meta']}

    with open(metadata_output, 'wt') as writer:
        yaml.dump(data, writer, default_flow_style=False)


def generate_cell_barcodes_file(bamfile, output):
    bamfile = pysam.AlignmentFile(bamfile, 'rb')
    header = bamfile.header

    cells = []
    for line in str(header).split('\n'):
        if not line.startswith("@CO"):
            continue
        line = line.strip().split()
        cb = line[1]
        cell = cb.split(':')[1]
        cells.append(cell)

    with open(output, 'wt') as writer:
        for cell in cells:
            writer.write(cell + '\n')


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    snv_genotyper = subparsers.add_parser('snv_genotyper')
    snv_genotyper.set_defaults(which='snv_genotyper')
    snv_genotyper.add_argument('--bam', required=True)
    snv_genotyper.add_argument('--output', required=True)
    snv_genotyper.add_argument('--targets_vcf', required=True)
    snv_genotyper.add_argument('--cell_barcodes')
    snv_genotyper.add_argument('--interval')
    snv_genotyper.add_argument('--count_duplicates', default=False)
    snv_genotyper.add_argument('--sparse', default=False)
    snv_genotyper.add_argument('--ignore_untagged_reads', action='store_true', default=False)
    snv_genotyper.add_argument('--min_mqual', default=20)

    parse_vartrix = subparsers.add_parser('parse_vartrix')
    parse_vartrix.set_defaults(which='parse_vartrix')
    parse_vartrix.add_argument(
        '--barcodes'
    )
    parse_vartrix.add_argument(
        '--variants'
    )
    parse_vartrix.add_argument(
        '--ref_counts'
    )
    parse_vartrix.add_argument(
        '--alt_counts'
    )
    parse_vartrix.add_argument(
        '--outfile'
    )
    parse_vartrix.add_argument(
        '--skip_header',
        action='store_true',
        default=False
    )
    parse_vartrix.add_argument(
        '--sparse',
        action='store_true',
        default=False
    )

    generate_cell_barcodes = subparsers.add_parser('generate_cell_barcodes')
    generate_cell_barcodes.set_defaults(which='generate_cell_barcodes')
    generate_cell_barcodes.add_argument(
        '--bamfile'
    )
    generate_cell_barcodes.add_argument(
        '--output'
    )

    generate_metadata = subparsers.add_parser('generate_metadata')
    generate_metadata.set_defaults(which='generate_metadata')
    generate_metadata.add_argument(
        '--outputs', nargs=2
    )
    generate_metadata.add_argument(
        '--vartrix_outputs', nargs=2
    )
    generate_metadata.add_argument(
        '--metadata_input'
    )
    generate_metadata.add_argument(
        '--metadata_output'
    )

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'snv_genotyper':
        with SnvGenotyper(
                args['bam'], args['targets_vcf'], args['output'], cell_barcodes=args['cell_barcodes'],
                interval=args['interval'], count_duplicates=args['count_duplicates'],
                sparse=args['sparse'], min_mqual=args['min_mqual'], ignore_untagged_reads=args['ignore_untagged_reads']
        ) as genotyper:
            genotyper.genotyping()
    elif args['which'] == 'generate_metadata':
        generate_metadata(
            args['outputs'], args['vartrix_outputs'], args['metadata_input'], args['metadata_output']
        )
    elif args['which'] == "parse_vartrix":
        parse_vartrix(
            args['barcodes'], args['variants'], args['ref_counts'],
            args['alt_counts'], args['outfile'],
            write_header=(not args['skip_header']),
            sparse=args['sparse']
        )
    elif args['which'] == "generate_cell_barcodes":
        generate_cell_barcodes_file(args['bamfile'], args['output'])
    else:
        raise Exception()
