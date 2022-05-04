import os

import argparse
import yaml
from mondrianutils import __version__
from mondrianutils.sv_genotyping.sv_genotyper import SvGenotyper


def generate_metadata(
        outputs, metadata_input, metadata_output
):
    data = dict()
    data['files'] = {
        os.path.basename(outputs[0]): {'result_type': 'pysam_genotyping_counts', 'auxiliary': False},
        os.path.basename(outputs[1]): {'result_type': 'pysam_genotyping_counts', 'auxiliary': True},
    }

    with open(metadata_input, 'rt') as reader:
        meta = yaml.safe_load(reader)
        del meta['meta']['type']
        del meta['meta']['version']

    meta_dict = {
        'type': 'sv_genotyping',
        'version': __version__
    }

    data['meta'] = {**meta_dict, **meta['meta']}

    with open(metadata_output, 'wt') as writer:
        yaml.dump(data, writer, default_flow_style=False)


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    snv_genotyper = subparsers.add_parser('sv_genotyper')
    snv_genotyper.set_defaults(which='sv_genotyper')
    snv_genotyper.add_argument('--bam', required=True)
    snv_genotyper.add_argument('--destruct_reads', required=True)
    snv_genotyper.add_argument('--destruct_table', required=True)
    snv_genotyper.add_argument('--output')

    generate_metadata = subparsers.add_parser('generate_metadata')
    generate_metadata.set_defaults(which='generate_metadata')
    generate_metadata.add_argument(
        '--outputs', nargs=2
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

    if args['which'] == 'sv_genotyper':
        genotyper = SvGenotyper(
            args['bam'], args['destruct_reads'], args['destruct_table'], args['output']
        )
        genotyper.main()
    elif args['which'] == 'generate_metadata':
        generate_metadata(
            args['outputs'], args['metadata_input'], args['metadata_output']
        )
    else:
        raise Exception()
