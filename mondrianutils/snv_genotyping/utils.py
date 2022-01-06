import os

import argparse
import yaml
from mondrianutils.snv_genotyping.snv_genotyper import SnvGenotyper


def generate_metadata(
        outputs, metadata_input, metadata_output
):
    data = dict()
    data['files'] = {
        os.path.basename(outputs[0]): {'result_type': 'snv_genotyping_counts'},
        os.path.basename(outputs[1]): {'result_type': 'snv_genotyping_counts'},
    }

    with open(metadata_input, 'rt') as reader:
        meta = yaml.safe_load(reader)
        del meta['meta']['name']
        del meta['meta']['version']

    meta_dict = {
        'name': 'snv_genotyping',
        'version': 'v0.0.9'
    }

    data['meta'] = {**meta_dict, **meta['meta']}

    with open(metadata_output, 'wt') as writer:
        yaml.dump(data, writer, default_flow_style=False)


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
    snv_genotyper.add_argument('--interval')
    snv_genotyper.add_argument('--count_duplicates', default=False)
    snv_genotyper.add_argument('--sparse', default=False)
    snv_genotyper.add_argument('--min_mqual', default=20)


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

    if args['which'] == 'snv_genotyper':
        with SnvGenotyper(
                args['bam'], args['targets_vcf'], args['output'],
                interval=args['interval'], count_duplicates=args['count_duplicates'],
                sparse=args['sparse'], min_mqual=args['min_mqual']
        ) as genotyper:
            genotyper.genotyping()
    elif args['which'] == 'generate_metadata':
        generate_metadata(
            args['outputs'], args['metadata_input'], args['metadata_output']
        )
    else:
        raise Exception()
