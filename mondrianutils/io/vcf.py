import os

from mondrianutils import helpers

import argparse

def get_num_calls(filepath):
    num_lines = 0
    with helpers.getFileHandle(filepath, 'rt') as reader:
        for line in reader:
            if not line.startswith('#'):
                num_lines += 1
    return num_lines


def get_header(filepath):
    header = []
    with helpers.getFileHandle(filepath, 'rt') as reader:
        for line in reader:
            if line.startswith('#'):
                header.append(line)
            else:
                break
    return header


def get_calls_grouped(filepath, num_calls):
    calls = []
    with helpers.getFileHandle(filepath, 'rt') as reader:
        for line in reader:
            if line.startswith('#'):
                continue

            if len(calls) >= num_calls:
                yield calls
                calls = []
            else:
                calls.append(line)

    yield calls


def split_vcf(infile, outdir, num_splits):
    num_lines = get_num_calls(infile)
    calls_per_file = num_lines // num_splits
    header = get_header(infile)

    helpers.makedirs(outdir)

    for i, calls in enumerate(get_calls_grouped(infile, calls_per_file)):
        outfile = os.path.join(outdir, '{}.vcf'.format(i))
        with open(outfile, 'wt') as writer:
            for line in header:
                writer.write(line)
            for line in calls:
                writer.write(line)


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    split_vcf = subparsers.add_parser('split_vcf')
    split_vcf.set_defaults(which='split_vcf')
    split_vcf.add_argument('--infile', required=True)
    split_vcf.add_argument('--outdir', required=True)
    split_vcf.add_argument('--num_splits', type=int, required=True)

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'split_vcf':
        split_vcf(
            args['infile'], args['outdir'], args['num_splits']
        )
    else:
        raise Exception()
