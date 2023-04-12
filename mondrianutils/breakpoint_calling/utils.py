import json

import argparse
import yaml
from mondrianutils import helpers

from . import consensus
from mondrianutils.breakpoint_calling import destruct_csv_to_vcf


def infer_type(files):
    with open(files, 'rt') as reader:
        files = json.load(reader)

    filetypes = sorted(set([v['left'] for v in files]))

    # more than one wf
    if 'svaba_vcf' in filetypes and 'breakpoint_consensus' in filetypes:
        return 'breakpoint_calling'
    elif 'breakpoint_consensus' in filetypes:
        return 'breakpoint_consensus'
    elif 'lumpy_vcf' in filetypes:
        return 'breakpoint_lumpy'
    elif 'svaba_vcf' in filetypes:
        return 'breakpoint_svaba'
    elif 'gridss_vcf' in filetypes:
        return 'breakpoint_gridss'
    elif 'destruct_calls' in filetypes:
        return 'breakpoint_destruct'
    else:
        raise Exception()


def generate_metadata(
        files, metadata_yaml_files, samples, metadata_output
):
    wf_type = infer_type(files)
    data = helpers.metadata_helper(files, metadata_yaml_files, samples, wf_type)

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
    consensus.add_argument('--region')
    consensus.add_argument('--blacklist_bed')

    csv_to_vcf = subparsers.add_parser('destruct_csv_to_vcf')
    csv_to_vcf.set_defaults(which='destruct_csv_to_vcf')
    csv_to_vcf.add_argument('--infile', required=True)
    csv_to_vcf.add_argument('--outfile', required=True)
    csv_to_vcf.add_argument('--reference', required=True)
    csv_to_vcf.add_argument('--sample_id', required=True)


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

    if args['which'] == 'consensus':
        consensus.consensus(
            args['destruct'], args['lumpy'], args['svaba'], args['gridss'],
            args['consensus'], args['sample_id'], args['tempdir'], region=args['region'],
            blacklist_bed=args['blacklist_bed']
        )
    elif args['which'] == 'generate_metadata':
        generate_metadata(
            args['files'], args['metadata_yaml_files'], args['samples'],
            args['metadata_output']
        )
    elif args['which'] == 'destruct_csv_to_vcf':
        destruct_csv_to_vcf.destruct_csv_to_vcf(
            args['infile'], args['outfile'], args['reference'],
            args['sample_id']
        )
    else:
        raise Exception()
