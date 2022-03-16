import os

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


def split_bam_by_barcode(infile, outdir):
    reader = pysam.AlignmentFile(infile, "rb")

    cells = get_cells(reader)

    outputs = init_out_bams(outdir, reader, cells)

    for pileupobj in reader.fetch(until_eof=True):
        cell_id = pileupobj.get_tag('CB')

        outputs[cell_id].write(pileupobj)

    close_all_bams(outputs)


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    split_bam = subparsers.add_parser('split_bam_by_barcode')
    split_bam.set_defaults(which='split_bam_by_barcode')
    split_bam.add_argument('--infile', required=True)
    split_bam.add_argument('--outdir', required=True)

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'split_bam_by_barcode':
        split_bam_by_barcode(
            args['infile'], args['outdir']
        )
    else:
        raise Exception()
