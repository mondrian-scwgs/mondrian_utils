'''
Created on Feb 19, 2018

@author: dgrewal
'''
import errno
import gzip
import json
import logging
import os
import subprocess
import tarfile
from subprocess import Popen, PIPE

import csverve.api as csverve
import numpy as np
import pandas as pd
import pysam
import yaml
from mondrianutils import __version__


def parse_region(region):
    if region is None:
        return None, None, None

    if ':' not in region:
        return region, None, None

    chrom, coords = region.split(':')

    assert '-' in coords

    beg, end = coords.split('-')

    return chrom, int(beg), int(end)


def get_chr_lengths(bam):
    names = bam.references
    lengths = bam.lengths
    return {name: length for name, length in zip(names, lengths)}


def get_bins_per_chromosome(chromosome_length, binsize):
    start = 1
    end = start - 1 + binsize
    bins = [(start, end)]
    while end < chromosome_length:
        start += binsize
        end = start - 1 + binsize
        bins.append((start, end))
    return bins


def get_cells(bam_reader):
    cells = []
    header = bam_reader.header
    for line in str(header).split('\n'):
        if not line.startswith("@CO"):
            continue
        line = line.strip().split()
        cb = line[1]
        cell = cb.split(':')[1]
        cells.append(cell)
    return cells


def _sum_two_dataframes(df):
    def sum_nan(y):
        if np.isnan(y).all():
            return float("nan")
        return np.nansum(y)

    if len(df) == 1:
        return df.iloc[0]
    else:
        return df.apply(sum_nan)


def _sort_bins(bins):
    bins = bins.drop_duplicates()

    # just the sort order in here, dont need to change for mouse.
    # might need to extend if reference has more genomes than human
    chromosomes = [str(v) for v in range(1, 23)] + ['X', 'Y']
    if bins.iloc[0][0].startswith('chr'):
        chromosomes = ['chr' + v for v in chromosomes]

    bins["chr"] = pd.Categorical(bins["chr"], chromosomes)

    bins = bins.sort_values(['chr', 'start', ])

    bins = [tuple(v) for v in bins.values.tolist()]

    return bins


def _fix_all_chunks(chunks):
    colnames = list(chunks[0].columns.values)

    new_chunks = [chunks[0]]

    for chunk in chunks[1:]:
        current_colnames = list(chunk.columns.values)

        if not current_colnames == colnames:
            missing_cols = [v for v in colnames if v not in current_colnames]
            temp_df = pd.DataFrame(columns=missing_cols)
            chunk = chunk.join(temp_df)

            chunk.columns = colnames

        new_chunks.append(chunk)

    return new_chunks


def load_and_pivot_reads_data(filepath, column_name):
    df = []
    for chunk in csverve.read_csv(filepath, chunksize=1e6):

        for col in chunk.columns:
            if chunk[col].dtype.name == 'category':
                chunk[col] = chunk[col].astype('str')

        chunk["bin"] = list(zip(chunk.chr, chunk.start, chunk.end))

        chunk = chunk.pivot(
            index='cell_id',
            columns='bin',
            values=column_name)

        df.append(chunk)

    df = _fix_all_chunks(df)
    df = pd.concat(df)
    df = df.groupby(df.index).apply(_sum_two_dataframes)

    bins = pd.DataFrame(
        df.columns.values.tolist(),
        columns=[
            'chr',
            'start',
            'end'])

    bins = _sort_bins(bins)

    df = df[bins]

    return df


def get_sample_from_bam(infile):
    bam = pysam.AlignmentFile(infile, "rb")
    header = bam.header

    samples = set()
    for rg in header['RG']:
        samples.add(rg['SM'])

    assert len(samples) == 1

    return sorted(samples)[0]


def chunks(bamfiles, numcores):
    output = []
    for i in range(0, len(bamfiles), numcores):
        output.append(bamfiles[i:i + numcores])
    return output


def get_merge_command(bams, output, ncores=1):
    if len(bams) == 1:
        command = ['cp', bams[0], output]
    else:
        command = ['sambamba', 'merge', '-t', str(ncores), output]
        command.extend(bams)

    return command


def merge_bams(infiles, outfile, tempdir, ncores):
    assert len(infiles) > 0

    if len(infiles) < ncores * 2:
        run_cmd(get_merge_command(infiles, outfile, ncores=ncores))
        return

    chunked_infiles = chunks(list(infiles), ncores)

    commands = []
    outputs = []
    for i, chunk in enumerate(chunked_infiles):
        chunk_tempdir = os.path.join(tempdir, str(i))
        makedirs(chunk_tempdir)
        output = os.path.join(chunk_tempdir, 'merged.bam')
        outputs.append(output)
        commands.append(get_merge_command(chunk, output))

    parallel_temp_dir = os.path.join(tempdir, 'gnu_parallel_temp')
    run_in_gnu_parallel(commands, parallel_temp_dir, ncores)

    command = get_merge_command(outputs, outfile, ncores=ncores)
    run_cmd(command)


def untar(input_tar, outdir):
    makedirs(outdir)
    with tarfile.open(input_tar) as tar:
        tar.extractall(path=outdir)


def metadata_helper(files_json, metadata_yamls, samples, wf_type):
    with open(files_json, 'rt') as files_json:
        jsondata = json.load(files_json)

    files_dict = {}

    for item in jsondata:
        filetype = str(item['left'])
        filepaths = item['right']
        for filepath in filepaths:
            filepath = os.path.basename(str(filepath))
            files_dict[filepath] = {
                'result_type': filetype, 'auxiliary': get_auxiliary_files(filepath)
            }

    metadata = {'type': wf_type, 'version': __version__}

    assert len(samples) == len(metadata_yamls)
    for sample, metadata_yaml in zip(samples, metadata_yamls):
        with open(metadata_yaml, 'rt') as reader:
            meta = yaml.safe_load(reader)
            del meta['meta']['type']
            del meta['meta']['version']

        metadata[sample] = meta['meta']

    return {'files': files_dict, 'meta': metadata}


def get_auxiliary_files(filepath):
    if filepath.endswith('.yaml'):
        return True
    elif filepath.endswith('.csi'):
        return True
    elif filepath.endswith('.tbi'):
        return True
    elif filepath.endswith('.bai'):
        return True
    else:
        return False


def run_cmd(cmd, output=None):
    stdout = PIPE
    if output:
        stdout = open(output, "w")

    p = Popen(cmd, stdout=stdout, stderr=PIPE)

    cmdout, cmderr = p.communicate()
    retc = p.returncode

    cmdout = "" if cmdout is None else cmdout.decode('utf-8')
    cmderr = "" if cmderr is None else cmderr.decode('utf-8')

    if output:
        stdout.close()

    console_out = '-' * 30 + ' cmdout ' + '-' * 30 + '\n'
    console_out += '\n'.join(cmdout.split('\n')) + '\n'
    console_out += '-' * 70

    console_err = '-' * 30 + ' cmderr ' + '-' * 30 + '\n'
    console_err += '\n'.join(cmderr.split('\n')) + '\n'
    console_err += '-' * 70

    if retc:
        raise Exception("command failed.\n {}\n {}".format(console_out, console_err))

    print(console_out)
    print(console_err)


class getFileHandle(object):
    def __init__(self, filename, mode='rt'):
        self.filename = filename
        self.mode = mode

    def __enter__(self):
        if self.get_file_format(self.filename) in ["csv", 'plain-text']:
            self.handle = open(self.filename, self.mode)
        elif self.get_file_format(self.filename) == "gzip":
            self.handle = gzip.open(self.filename, self.mode)
        elif self.get_file_format(self.filename) == "h5":
            self.handle = pd.HDFStore(self.filename, self.mode)
        return self.handle

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.handle.close()

    def get_file_format(self, filepath):
        if filepath.endswith('.tmp'):
            filepath = filepath[:-4]

        _, ext = os.path.splitext(filepath)

        if ext == ".csv":
            return "csv"
        elif ext == ".gz":
            return "gzip"
        elif ext == ".h5" or ext == ".hdf5":
            return "h5"
        elif ext == '.yaml':
            return 'plain-text'
        else:
            logging.getLogger("single_cell.helpers").warning(
                "Couldn't detect output format. extension {}".format(ext)
            )
            return "plain-text"


def makedirs(directory, isfile=False):
    if isfile:
        directory = os.path.dirname(directory)
        if not directory:
            return

    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def build_shell_script(command, tag, tempdir):
    outfile = os.path.join(tempdir, "{}.sh".format(tag))
    with open(outfile, 'w') as scriptfile:
        scriptfile.write("#!/bin/bash\n")
        if isinstance(command, list) or isinstance(command, tuple):
            command = ' '.join(map(str, command)) + '\n'
        scriptfile.write(command)
    return outfile


def run_in_gnu_parallel(commands, tempdir, ncores):
    makedirs(tempdir)

    scriptfiles = []

    for tag, command in enumerate(commands):
        scriptfiles.append(build_shell_script(command, tag, tempdir))

    parallel_outfile = os.path.join(tempdir, "commands.txt")
    with open(parallel_outfile, 'w') as outfile:
        for scriptfile in scriptfiles:
            outfile.write("sh {}\n".format(scriptfile))

    subprocess.run(['parallel', '--jobs', str(ncores)], stdin=open(parallel_outfile))


def make_tarfile(output_filename, source_dir):
    assert output_filename.endswith('.tar.gz')

    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))
