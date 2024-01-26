import os
import re
import shutil
import pandas as pd
from itertools import islice
import csverve.api as csverve
from subprocess import Popen, PIPE
from collections import defaultdict
from mondrianutils import helpers
from mondrianutils.dtypes.alignment import dtypes


class FastqReader(object):

    def __init__(self, filepath):
        self.file_path = filepath

    def get_read_iterator(self):
        with helpers.getFileHandle(self.file_path) as fq_reader:
            while True:
                fastq_read = list(islice(fq_reader, 4))

                fastq_read = [line for line in fastq_read]

                if not fastq_read:
                    break

                assert len(fastq_read) == 4, 'fastq file format error'

                if not fastq_read[0].startswith('@'):
                    raise ValueError('Expected @ as first character of read name')

                if not fastq_read[2].startswith('+'):
                    raise ValueError('Expected = as first character of read comment')

                yield fastq_read


def _get_read_name(fastq_line1):
    return re.split('/| |\t|#FQST:', fastq_line1.split()[0])[0].rstrip()


class PairedFastqReader(object):
    def __init__(self, r1_path, r2_path):
        self.reader_r1 = FastqReader(r1_path)
        self.reader_r2 = FastqReader(r2_path)

    def get_read_pair_iterator(self):
        read_r1_iter = self.reader_r1.get_read_iterator()
        read_r2_iter = self.reader_r2.get_read_iterator()

        for read_r1, read_r2 in zip(read_r1_iter, read_r2_iter):

            if not read_r1 and not read_r2:
                break

            if not read_r1 or not read_r2:
                raise Exception('mismatching number of reads in R1 and R2')

            assert _get_read_name(read_r1[0]) == _get_read_name(read_r2[0])

            yield read_r1, read_r2


class TaggedFastqReader(FastqReader):
    def __init__(self, fastq_path):
        super(TaggedFastqReader, self).__init__(fastq_path)
        self.indices = None

    def get_read_tag(self, fastq_read):

        read_id = fastq_read[0]
        fq_tag = read_id[read_id.index('FQST:'):]
        fq_tag = fq_tag.strip().split(':')

        if not self.indices:
            if len(fq_tag) > 2:
                self.indices = {i: v for i, v in enumerate(fq_tag[1:-1])}
            else:
                raise Exception('First line in fastq file should have filter explanation')

        flag = map(int, list(fq_tag[-1]))

        flag_map = {self.indices[i]: v for i, v in enumerate(flag)}

        return flag_map

    def add_tag_to_read_comment(self, read, tag=None):
        read_name = _get_read_name(read[0])

        if not tag:
            tag = self.get_read_tag(read)

        tag = ['{}_{}'.format(k, v) for k, v in tag.items()]
        tag = ','.join(tag)
        tag = 'FS:Z:' + tag

        comment = tag

        read[0] = read_name + '\t' + comment + '\n'

        return read

    def filter_read_iterator(self, reference):
        for read in self.get_read_iterator():
            read_tags = self.get_read_tag(read)

            # skip if read maps to multiple genomes
            if len([v for v in read_tags.values() if v]) > 1:
                continue

            if read_tags[reference]:
                yield read
                continue

    def gather_counts(self):
        key_order = None
        counts = defaultdict(int)

        for read in self.get_read_iterator():

            read_tags = self.get_read_tag(read)
            if not key_order:
                key_order = sorted(read_tags.keys())

            flags = zip(key_order, [read_tags[key] for key in key_order])

            counts[flags] += 1

        return counts


class PairedTaggedFastqReader(PairedFastqReader, TaggedFastqReader):
    def __init__(self, fastq_r1, fastq_r2):
        super(PairedTaggedFastqReader, self).__init__(fastq_r1, fastq_r2)
        self.indices = None

    def filter_read_iterator(self, reference):
        for read_1, read_2 in self.get_read_pair_iterator():

            tags_r1 = self.get_read_tag(read_1)
            tags_r2 = self.get_read_tag(read_2)

            # skip if doesnt match
            if not tags_r1[reference] and not tags_r2[reference]:
                continue

            # skip if read maps to multiple genomes
            if len([v for v in tags_r1.values() if v]) > 1:
                continue
            if len([v for v in tags_r2.values() if v]) > 1:
                continue

            r1_nomatch = set(tags_r1.values()) == {0}
            r2_nomatch = set(tags_r2.values()) == {0}

            if tags_r1[reference] and tags_r2[reference]:
                yield read_1, read_2
            elif tags_r1[reference] and r2_nomatch:
                yield read_1, read_2
            elif r1_nomatch and tags_r2[reference]:
                yield read_1, read_2

    def gather_counts(self):
        key_order = None
        counts = {'R1': defaultdict(int), 'R2': defaultdict(int)}

        for read_1, read_2 in self.get_read_pair_iterator():
            tags_r1 = self.get_read_tag(read_1)
            tags_r2 = self.get_read_tag(read_2)

            if not key_order:
                key_order = sorted(tags_r1.keys())

            r1_flags = tuple(zip(key_order, [tags_r1[key] for key in key_order]))
            r2_flags = tuple(zip(key_order, [tags_r2[key] for key in key_order]))

            counts["R1"][r1_flags] += 1
            counts["R2"][r2_flags] += 1

        return counts


def merge_fastq_screen_counts(
        all_detailed_counts, all_summary_counts, merged_detailed_counts, merged_summary_counts
):
    helpers.makedirs(os.path.dirname(merged_summary_counts))
    helpers.makedirs(os.path.dirname(merged_detailed_counts))

    if isinstance(all_detailed_counts, dict):
        all_detailed_counts = all_detailed_counts.values()

    detailed_data = []
    for countsfile in all_detailed_counts:
        if os.stat(countsfile).st_size == 0:
            continue
        detailed_data.append(pd.read_csv(countsfile))

    df = pd.concat(detailed_data)

    index_cols = [v for v in df.columns.values if v != "count"]

    df['count'] = df.groupby(index_cols)['count'].transform('sum')

    df = df.drop_duplicates(subset=index_cols)

    fastqscreen_genomes = [v for v in df.columns.values if v not in ['cell_id', 'readend', 'count']]

    csverve.write_dataframe_to_csv_and_yaml(
        df, merged_detailed_counts,
        dtypes(fastqscreen_genomes=fastqscreen_genomes)['fastqscreen_detailed'],
        skip_header=False
    )

    if isinstance(all_summary_counts, dict):
        all_summary_counts = all_summary_counts.values()

    summary_counts = [pd.read_csv(countsfile) for countsfile in all_summary_counts]

    df = pd.concat(summary_counts)

    update_cols = [v for v in df.columns.values if v != 'cell_id']

    for colname in update_cols:
        df[colname] = df.groupby('cell_id')[colname].transform('sum')

    df = df.drop_duplicates(subset=['cell_id'])

    fastqscreen_genomes = [v for v in df.columns.values if v.startswith('fastqscreen_')]
    fastqscreen_genomes = sorted(set([v.split('_')[1] for v in fastqscreen_genomes]))
    fastqscreen_genomes = [v for v in fastqscreen_genomes if v not in ['nohit', 'total']]

    csverve.write_dataframe_to_csv_and_yaml(
        df, merged_summary_counts,
        dtypes(fastqscreen_genomes=fastqscreen_genomes)['metrics'], skip_header=False
    )


def run_cmd(cmd, output=None):
    stdout = PIPE
    if output:
        stdout = open(output, "w")

    p = Popen(cmd, stdout=stdout, stderr=PIPE)

    cmdout, cmderr = p.communicate()
    retc = p.returncode

    if retc:
        raise Exception(
            "command failed. stderr:{}, stdout:{}".format(
                cmdout,
                cmderr))

    if output:
        stdout.close()


def run_fastq_screen_paired_end(fastq_r1, fastq_r2, tempdir, params, num_threads):
    def get_basename(filepath):
        filepath_base = os.path.basename(filepath)

        if filepath_base.endswith('.fastq.gz'):
            filepath_base = filepath_base[:-len('.fastq.gz')]
        elif filepath_base.endswith('.fq.gz'):
            filepath_base = filepath_base[:-len('.fq.gz')]
        elif filepath_base.endswith('.fastq'):
            filepath_base = filepath_base[:-len('.fastq')]
        elif filepath_base.endswith('.fq'):
            filepath_base = filepath_base[:-len('.fq')]
        else:
            raise Exception('unknown file format. {}'.format(filepath))
        return filepath_base

    basename = get_basename(fastq_r1)
    tagged_fastq_r1 = os.path.join(tempdir, '{}.tagged.fastq.gz'.format(basename))

    basename = get_basename(fastq_r2)
    tagged_fastq_r2 = os.path.join(tempdir, '{}.tagged.fastq.gz'.format(basename))

    # fastq screen fails if run on empty files
    with helpers.getFileHandle(fastq_r1) as reader:
        if not reader.readline():
            shutil.copy(fastq_r1, tagged_fastq_r1)
            shutil.copy(fastq_r2, tagged_fastq_r2)
            return tagged_fastq_r1, tagged_fastq_r2

    config = os.path.join(tempdir, 'fastq_screen.config')

    with open(config, 'w') as config_writer:
        for genome in params['genomes']:
            genome_name = genome['name']
            genome_path = genome['path']
            outstr = '\t'.join(['DATABASE', genome_name, genome_path]) + '\n'
            config_writer.write(outstr)

        if not num_threads == 1:
            outstr = '\t'.join(['THREADS', str(num_threads)]) + '\n'
            config_writer.write(outstr)

    cmd = [
        'fastq_screen',
        '--aligner', params['aligner'],
        '--conf', config,
        '--outdir', tempdir,
        '--tag',
        fastq_r1,
        fastq_r2,
    ]

    run_cmd(cmd)

    return tagged_fastq_r1, tagged_fastq_r2


def write_detailed_counts(counts, outfile, cell_id, fastqscreen_params):
    header = None

    genomes = [genome['name'] for genome in fastqscreen_params['genomes']]

    with helpers.getFileHandle(outfile, 'wt') as writer:

        for read_end, read_end_counts in counts.items():

            if not read_end_counts and not header:
                outstr = ['cell_id', 'readend'] + genomes + ['count']
                writer.write(','.join(outstr) + '\n')
                header = 1
                continue

            if not header:
                outstr = ['cell_id', 'readend']
                outstr += [v[0] for v in list(read_end_counts.keys())[0]]
                outstr += ['count']
                writer.write(','.join(outstr) + '\n')
                header = 1

            for flags, count in read_end_counts.items():
                outstr = [cell_id, read_end]
                outstr += [v[1] for v in flags]
                outstr += [count]
                writer.write(','.join(map(str, outstr)) + '\n')


def write_summary_counts(counts, outfile, cell_id, fastqscreen_params):
    genomes = [genome['name'] for genome in fastqscreen_params['genomes']]

    summary_counts = {'nohit': 0, 'total_reads': 0}
    for genome in genomes:
        summary_counts[genome] = 0
        summary_counts['{}_multihit'.format(genome)] = 0

    for read_end, read_end_counts in counts.items():
        for flags, count in read_end_counts.items():
            summary_counts['total_reads'] += count
            hit_orgs = [v[0] for v in flags if v[1] > 0]

            for org in hit_orgs:
                summary_counts[org] += count

            if len(hit_orgs) > 1:
                for org in hit_orgs:
                    summary_counts['{}_multihit'.format(org)] += count
            elif len(hit_orgs) == 0:
                summary_counts['nohit'] += count

    with helpers.getFileHandle(outfile, 'wt') as writer:
        if not summary_counts:
            columns = ['cell_id']
            columns += ['fastqscreen_' + genome for genome in genomes]
            columns += ['fastqscreen_nohit', 'fastqscreen_total_reads']
            header = ','.join(columns) + '\n'
            writer.write(header)
            data = [0] * len(columns)
            data[0] = cell_id
            data = [str(v) for v in data]
            data = ','.join(data) + '\n'
            writer.write(data)
            return

        keys = sorted(summary_counts.keys())
        header = ['cell_id'] + ['fastqscreen_{}'.format(key) for key in keys]
        header = ','.join(header) + '\n'
        writer.write(header)

        values = [cell_id] + [summary_counts[v] for v in keys]
        values = ','.join(map(str, values)) + '\n'
        writer.write(values)


def re_tag_reads(infile, outfile):
    reader = TaggedFastqReader(infile)

    with helpers.getFileHandle(outfile, 'wt') as writer:

        for read in reader.get_read_iterator():
            read = reader.add_tag_to_read_comment(read)

            for line in read:
                writer.write(line)


def organism_filter(
        fastq_r1, fastq_r2, filtered_fastq_r1, filtered_fastq_r2,
        detailed_metrics, summary_metrics, tempdir, cell_id,
        reference, reference_name, supplementary_references, supplementary_reference_names,
        num_threads
):
    genomes = [{'name': reference_name, 'path': reference}]

    for ref, name in zip(supplementary_references, supplementary_reference_names):
        genomes.append(
            {'name': name, 'path': ref}
        )

    params = {
        'strict_validation': True,
        'filter_contaminated_reads': False,
        'aligner': 'bwa',
        'genomes': genomes
    }

    # fastq screen tries to skip if files from old runs are available
    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)

    helpers.makedirs(tempdir)

    tagged_fastq_r1, tagged_fastq_r2 = run_fastq_screen_paired_end(
        fastq_r1, fastq_r2, tempdir, params, num_threads
    )

    reader = PairedTaggedFastqReader(tagged_fastq_r1, tagged_fastq_r2)
    counts = reader.gather_counts()

    write_detailed_counts(counts, detailed_metrics, cell_id, params)
    write_summary_counts(counts, summary_metrics, cell_id, params)

    # use the full tagged fastq downstream
    # with organism type information in readname
    re_tag_reads(tagged_fastq_r1, filtered_fastq_r1)
    re_tag_reads(tagged_fastq_r2, filtered_fastq_r2)
