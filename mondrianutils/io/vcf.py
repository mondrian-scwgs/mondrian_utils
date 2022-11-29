import os

import argparse
import numpy as np
import pandas as pd
from mondrianutils import helpers
from mondrianutils.io import vcf_merge


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
                calls = [line]
            else:
                calls.append(line)

    yield calls


def split_vcf(infile, outdir, num_splits):
    num_lines = get_num_calls(infile)
    calls_per_file = max(1, num_lines // num_splits)

    header = get_header(infile)

    helpers.makedirs(outdir)

    for i, calls in enumerate(get_calls_grouped(infile, calls_per_file)):
        outfile = os.path.join(outdir, '{}.vcf'.format(i))
        with open(outfile, 'wt') as writer:
            for line in header:
                writer.write(line)
            for line in calls:
                writer.write(line)


def split_vcf_by_chrom(infile, outdir):
    header = get_header(infile)

    helpers.makedirs(outdir)

    filehandles = {}

    with helpers.getFileHandle(infile, 'rt') as reader:
        for line in reader:
            if line.startswith('#'):
                continue

            chrom = line.split()[0]
            outfile = os.path.join(outdir, f'{chrom}.vcf')

            if outfile in filehandles:
                filehandles[outfile].write(line)
            else:
                filehandles[outfile] = open(outfile, 'wt')
                [filehandles[outfile].write(v) for v in header]
                filehandles[outfile].write(line)

    for _, filehandle in filehandles.items():
        filehandle.close()


def remove_duplicates(input_vcf, output_vcf, include_ref_alt=False):
    seen = set()

    with helpers.getFileHandle(input_vcf, 'rt') as reader, helpers.getFileHandle(output_vcf, 'wt') as writer:

        for line in reader:
            if line.startswith('#'):
                writer.write(line)
                continue

            line_split = line.strip().split()
            chrom = line_split[0]
            pos = line_split[1]
            ref = line_split[3]
            alt = line_split[4]

            if include_ref_alt:
                key = (chrom, pos, ref, alt)
            else:
                key = (chrom, pos)

            if key in seen:
                continue

            seen.add(key)
            writer.write(line)


def _get_chrom_excluded(excluded, chrom):
    chrom_length = max(excluded[excluded['chrom'] == chrom]['end']) + 1
    chrom_excluded = np.zeros(chrom_length + 1, dtype=np.uint8)

    for start, end in excluded.loc[excluded['chrom'] == chrom, ['start', 'end']].values:
        start = min(start, chrom_length)
        end = min(end, chrom_length)
        chrom_excluded[start:end] = 1

    return chrom_excluded


def exclude_blacklist(input_vcf, output_vcf, exclusion_blacklist):
    excluded = pd.read_csv(exclusion_blacklist, sep="\t", )
    excluded.columns = ["chrom", "start", "end"]

    chrom_excluded = None

    with helpers.getFileHandle(input_vcf, 'rt') as reader, helpers.getFileHandle(output_vcf, 'wt') as writer:

        for line in reader:
            if line.startswith('#'):
                writer.write(line)
                continue

            line_split = line.strip().split()
            chrom = line_split[0]
            pos = int(line_split[1])

            if chrom_excluded is None or not chrom_excluded[0] == chrom:
                chrom_excluded = (chrom, _get_chrom_excluded(excluded, chrom))

            if chrom_excluded[1][pos]:
                continue

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

    split_vcf_by_chrom = subparsers.add_parser('split_vcf_by_chrom')
    split_vcf_by_chrom.set_defaults(which='split_vcf_by_chrom')
    split_vcf_by_chrom.add_argument('--infile', required=True)
    split_vcf_by_chrom.add_argument('--outdir', required=True)

    remove_duplicates = subparsers.add_parser('remove_duplicates')
    remove_duplicates.set_defaults(which='remove_duplicates')
    remove_duplicates.add_argument('--infile', required=True)
    remove_duplicates.add_argument('--outfile', required=True)
    remove_duplicates.add_argument('--include_ref_alt', action='store_true', default=False)

    exclude_blacklist = subparsers.add_parser('exclude_blacklist')
    exclude_blacklist.set_defaults(which='exclude_blacklist')
    exclude_blacklist.add_argument('--infile', required=True)
    exclude_blacklist.add_argument('--outfile', required=True)
    exclude_blacklist.add_argument('--exclusion_blacklist', default=False)

    merge_vcfs = subparsers.add_parser('merge_vcfs')
    merge_vcfs.set_defaults(which='merge_vcfs')
    merge_vcfs.add_argument('--infiles', nargs='*', required=True)
    merge_vcfs.add_argument('--outfile', required=True)

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'split_vcf':
        split_vcf(
            args['infile'], args['outdir'], args['num_splits']
        )
    elif args['which'] == 'split_vcf_by_chrom':
        split_vcf_by_chrom(
            args['infile'], args['outdir']
        )
    elif args['which'] == 'remove_duplicates':
        remove_duplicates(
            args['infile'], args['outfile'], include_ref_alt=args['include_ref_alt']
        )
    elif args['which'] == 'merge_vcfs':
        vcf_merge.merge_vcfs(
            args['infiles'], args['outfile']
        )
    elif args['which'] == 'exclude_blacklist':
        exclude_blacklist(
            args['infile'], args['outfile'], args['exclusion_blacklist']
        )
    else:
        raise Exception()
