import os

import argparse
import yaml
from mondrianutils import __version__

from . import consensus


def generate_metadata(
        consensus_files, destruct_files, destruct_reads_files,
        destruct_library_files, lumpy_vcf_files, svaba_vcf_files,
        gridss_vcf_files, metadata_yaml_files, samples, metadata_output
):
    assert len(samples) == len(metadata_yaml_files)
    data = dict()
    data['files'] = {
        os.path.basename(consensus_files[0]): {'result_type': 'breakpoint_consensus'},
        os.path.basename(consensus_files[1]): {'result_type': 'breakpoint_consensus'},
    }

    for destruct_file in destruct_files:
        destruct_file = os.path.basename(destruct_file)
        data['files'][destruct_file] = {'result_type': 'destruct_calls'}

    for destruct_reads_file in destruct_reads_files:
        destruct_reads_file = os.path.basename(destruct_reads_file)
        data['files'][destruct_reads_file] = {'result_type': 'destruct_reads'}

    for destruct_library_file in destruct_library_files:
        destruct_library_file = os.path.basename(destruct_library_file)
        data['files'][destruct_library_file] = {'result_type': 'destruct_library'}

    for lumpy_vcf in lumpy_vcf_files:
        lumpy_vcf = os.path.basename(lumpy_vcf)
        data['files'][lumpy_vcf] = {'result_type': 'lumpy_vcf'}

    for gridss_vcf in gridss_vcf_files:
        gridss_vcf = os.path.basename(gridss_vcf)
        data['files'][gridss_vcf] = {'result_type': 'gridss_vcf'}

    for svaba_vcf in svaba_vcf_files:
        svaba_vcf = os.path.basename(svaba_vcf)
        data['files'][svaba_vcf] = {'result_type': 'svaba_vcf'}

    data['meta'] = {
        'name': 'breakpoint_calling',
        'version': __version__,
    }
    for sample, metadata_yaml in zip(samples, metadata_yaml_files):
        with open(metadata_yaml, 'rt') as reader:
            meta = yaml.safe_load(reader)
            del meta['meta']['name']
            del meta['meta']['version']

        data['meta'][sample] = meta['meta']

    with open(metadata_output, 'wt') as writer:
        yaml.dump(data, writer, default_flow_style=False)


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    consensus = subparsers.add_parser('consensus')
    consensus.set_defaults(which='consensus')
    consensus.add_argument('--destruct', required=True)
    consensus.add_argument('--lumpy', required=True)
    consensus.add_argument('--gridss', required=True)
    consensus.add_argument('--svaba', required=True)
    consensus.add_argument('--sample_id', required=True)
    consensus.add_argument('--consensus', required=True)
    consensus.add_argument('--tempdir', required=True)

    generate_metadata = subparsers.add_parser('generate_metadata')
    generate_metadata.set_defaults(which='generate_metadata')
    generate_metadata.add_argument(
        '--consensus_files', nargs=2
    )
    generate_metadata.add_argument(
        '--destruct_files', nargs='*'
    )
    generate_metadata.add_argument(
        '--destruct_reads_files', nargs='*'
    )
    generate_metadata.add_argument(
        '--destruct_library_files', nargs='*'
    )
    generate_metadata.add_argument(
        '--lumpy_vcf_files', nargs='*'
    )
    generate_metadata.add_argument(
        '--svaba_vcf_files', nargs='*'
    )
    generate_metadata.add_argument(
        '--gridss_vcf_files', nargs='*'
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

    if args['which'] == 'consensus':
        consensus.consensus(
            args['destruct'], args['lumpy'], args['svaba'], args['gridss'],
            args['consensus'], args['sample_id'], args['tempdir']
        )
    elif args['which'] == 'generate_metadata':
        generate_metadata(
            args['consensus_files'], args['destruct_files'], args['destruct_reads_files'],
            args['destruct_library_files'], args['lumpy_vcf_files'], args['svaba_vcf_files'],
            args['gridss_vcf_files'], args['metadata_yaml_files'], args['samples'],
            args['metadata_output']
        )
    else:
        raise Exception()
