import os

import click
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


def split_vcf_into_numsplits(infile, outdir, num_splits):
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


def split_vcf_by_lines(infile, outdir, num_lines):
    header = get_header(infile)

    helpers.makedirs(outdir)

    for i, calls in enumerate(get_calls_grouped(infile, num_lines)):
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

            if pos < len(chrom_excluded[1]) and chrom_excluded[1][pos]:
                continue

            writer.write(line)

