import os
import shutil

import argparse
import pysam
from mondrianutils import helpers


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


def _split_bam_by_barcode_serial(infile, outdir):
    reader = pysam.AlignmentFile(infile, "rb")

    cells = get_cells(reader)

    outputs = init_out_bams(outdir, reader, cells)

    for pileupobj in reader.fetch(until_eof=True):
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
            'bam_utils', 'split_bam_by_barcode',
            '--infile', infile, '--outdir',
            os.path.join(outdir, str(chromosome)), '--chromosome', chromosome, '--ncores', 1
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
        _split_bam_by_barcode_serial(infile, outdir)
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


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    split_bam = subparsers.add_parser('split_bam_by_barcode')
    split_bam.set_defaults(which='split_bam_by_barcode')
    split_bam.add_argument('--infile', required=True)
    split_bam.add_argument('--outdir', required=True)
    split_bam.add_argument('--tempdir')
    split_bam.add_argument('--chromosomes', nargs='*')
    split_bam.add_argument('--ncores', default=8, type=int)

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'split_bam_by_barcode':
        split_bam_by_barcode(
            args['infile'], args['outdir'], args['tempdir'],
            args['chromosomes'], ncores=args['ncores']
        )
    else:
        raise Exception()
