import os
import shutil

import csverve.api as csverve
import pandas as pd
import pysam
from mondrianutils import helpers
from mondrianutils.dtypes.hmmcopy import dtypes as hmmcopy_dtypes


def _filter_read(read, mapping_quality_threshold=20):
    if read.mapping_quality < mapping_quality_threshold:
        return True
    elif read.is_secondary:
        return True
    elif read.is_supplementary:
        return True
    elif read.is_duplicate:
        return True
    elif read.is_unmapped:
        return True
    elif read.is_qcfail:
        return True
    else:
        return False


def _iterate_reads_in_pairs(bam, chromosome, start, end, mapping_quality_threshold=20):
    read_pairs = {}
    reads_processed = set()

    for read in bam.fetch(chromosome, start, end):

        if _filter_read(read, mapping_quality_threshold=mapping_quality_threshold):
            continue

        cell_id = read.get_tag('CB')
        qname = read.query_name

        assert qname not in reads_processed

        if (cell_id, qname) not in read_pairs:
            read_pairs[(cell_id, qname)] = read.get_reference_positions()
        else:
            reads_processed.add(qname)
            yield cell_id, set(read_pairs.pop((cell_id, qname)) + read.get_reference_positions())

    for (cell_id, qname), positions in read_pairs.items():
        yield cell_id, positions


def get_per_cell_fraction(bam, chromosome, start, end, mapping_quality_threshold=20):
    data = {}
    iterator = _iterate_reads_in_pairs(
        bam, chromosome, start, end,
        mapping_quality_threshold=mapping_quality_threshold
    )

    for cell_id, positions in iterator:

        if cell_id not in data:
            data[cell_id] = [0] * (end - start + 1)

        for pos in positions:
            if start <= pos <= end:
                data[cell_id][pos - start] += 1

    fractions = {}
    for cell, celldata in data.items():
        overlapping = len([v for v in celldata if v > 1])
        fractions[cell] = (overlapping / (end - start))

    return fractions


def _overlapping_fraction_per_bin_serial(
        bam,
        output,
        chromosomes=[str(v) for v in range(1, 23)] + ['X', 'Y'],
        binsize=500000,
        mapping_quality=20
):
    bam = pysam.AlignmentFile(bam, 'rb')

    cells = helpers.get_cells(bam)

    chrlengths = helpers.get_chr_lengths(bam)

    total_data = []
    for chromosome in chromosomes:
        bins = helpers.get_bins_per_chromosome(chrlengths[chromosome], binsize)

        for (start, end) in bins:
            fractions = get_per_cell_fraction(bam, chromosome, start, end, mapping_quality_threshold=mapping_quality)

            for cell in cells:
                total_data.append([chromosome, start, end, cell, fractions.get(cell, 0)])

    df = pd.DataFrame(total_data, columns=['chr', 'start', 'end', 'cell_id', 'fraction_overlapping_reads'])

    csverve.write_dataframe_to_csv_and_yaml(
        df, output, hmmcopy_dtypes()['reads']
    )


def overlapping_fraction_per_bin(
        bam, output, tempdir,
        chromosomes=[str(v) for v in range(1, 23)] + ['X', 'Y'],
        binsize=500000, mapping_quality=20, ncores=8
):
    if ncores == 1:
        _overlapping_fraction_per_bin_serial(
            bam, output,
            chromosomes=chromosomes,
            binsize=binsize,
            mapping_quality=mapping_quality
        )
    else:
        if tempdir is None:
            raise Exception('please specify --tempdir')
        if os.path.exists(tempdir):
            shutil.rmtree(tempdir)
        os.makedirs(tempdir)

        split_commands = []
        for chromosome in chromosomes:
            cmd = [
                'io_utils', 'overlapping-fraction-per-bin',
                '--bam', bam, '--chromosomes', chromosome,
                '--output', os.path.join(tempdir, f'output_{chromosome}.csv.gz'),
                '--binsize', binsize, '--mapping_quality', mapping_quality,
                '--ncores', 1
            ]
            split_commands.append(cmd)

        helpers.run_in_gnu_parallel(split_commands, tempdir, ncores)

        csverve.concatenate_csv(
            [os.path.join(tempdir, f'output_{v}.csv.gz') for v in chromosomes],
            output
        )


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


def init_out_bams(outdir, reader, cells):
    helpers.makedirs(outdir)

    outdata = {}

    for cell in cells:
        outpath = os.path.join(outdir, "{}.bam".format(cell))

        outfile = pysam.AlignmentFile(outpath, "wb", template=reader)

        outdata[cell] = outfile

    return outdata


def close_all_bams(out_bams):
    for cellid, writer in out_bams.items():
        writer.close()


def _split_bam_by_barcode_serial(infile, outdir, chromosomes):
    reader = pysam.AlignmentFile(infile, "rb")

    cells = get_cells(reader)

    outputs = init_out_bams(outdir, reader, cells)

    if chromosomes is None or chromosomes == []:
        chromosomes = helpers.get_chr_lengths(reader).keys()

    for chromosome in chromosomes:
        for pileupobj in reader.fetch(chromosome, until_eof=True):
            cell_id = pileupobj.get_tag('CB')

            outputs[cell_id].write(pileupobj)

    close_all_bams(outputs)


def merge_cmd(bams, output):
    if len(bams) == 1:
        command = ['cp', bams[0], output]
    else:
        command = [
            'sambamba', 'merge', output
        ]
        command.extend(bams)

    return command


def _split_in_parallel(infile, outdir, chromosomes, tempdir, ncores=8):
    split_commands = []
    for chromosome in chromosomes:
        cmd = [
            'io_utils', 'split-bam-by-barcode',
            '--infile', infile, '--outdir',
            os.path.join(outdir, str(chromosome)), '--chromosomes ', chromosome, '--ncores', 1
        ]
        split_commands.append(cmd)

    split_scripts_tempdir = os.path.join(tempdir, 'scripts_split')

    helpers.run_in_gnu_parallel(split_commands, split_scripts_tempdir, ncores)


def _merge_in_parallel(datadir, outdir, tempdir, cells, chromosomes, ncores=8):
    merge_commands = []
    for cell in cells:

        bam_files = [os.path.join(datadir, chromosome, f'{cell}.bam') for chromosome in chromosomes]

        for bamfile in bam_files:
            assert os.path.exists(bamfile)

        cmd = merge_cmd(bam_files, os.path.join(outdir, f'{cell}.bam'))
        merge_commands.append(cmd)

    helpers.run_in_gnu_parallel(merge_commands, tempdir, ncores)


def split_bam_by_barcode(infile, outdir, tempdir, chromosomes, ncores=8):
    if ncores == 1:
        _split_bam_by_barcode_serial(infile, outdir, chromosomes)
        return
    else:
        if tempdir is None:
            raise Exception('please specify --tempdir')
        if chromosomes is None:
            raise Exception('please specify --chromosomes')

    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)
    if os.path.exists(outdir):
        shutil.rmtree(outdir)

    os.makedirs(tempdir)
    os.makedirs(outdir)

    cells = get_cells(pysam.AlignmentFile(infile, 'rb'))

    _split_in_parallel(
        infile,
        os.path.join(tempdir, 'files_split'),
        chromosomes,
        os.path.join(tempdir, 'scripts_split'),
        ncores=ncores
    )

    _merge_in_parallel(
        os.path.join(tempdir, 'files_split'),
        outdir,
        os.path.join(tempdir, 'scripts_merge'),
        cells,
        chromosomes,
        ncores=ncores
    )
