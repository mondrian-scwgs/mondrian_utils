import os

import argparse
import pysam
import yaml
from mondrianutils import helpers
from mondrianutils.io import identify_normal_cells
import copy


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


def update_header_comments(header, cells):
    if not isinstance(header, dict):
        newheader = copy.deepcopy(header.to_dict())
    else:
        newheader = copy.deepcopy(header)

    newheader['CO'] = [f'CB:{v}' for v in cells]

    return newheader


def separate_normal_and_tumour_cells(infile, normal_output, tumour_output, normal_cells_yaml):
    with open(normal_cells_yaml, 'rt') as reader:
        normal_cells = set(yaml.safe_load(reader))

    reader = pysam.AlignmentFile(infile, "rb")

    tumour_cells = [v for v in get_cells(reader) if v not in normal_cells]

    normal_header = update_header_comments(reader.header, normal_cells)
    tumour_header = update_header_comments(reader.header, tumour_cells)

    normal_writer = pysam.AlignmentFile(normal_output, "wb", header=normal_header)
    tumour_writer = pysam.AlignmentFile(tumour_output, "wb", header=tumour_header)

    for pileupobj in reader.fetch(until_eof=True):
        cell_id = pileupobj.get_tag('CB')

        if cell_id in normal_cells:
            normal_writer.write(pileupobj)
        else:
            tumour_writer.write(pileupobj)


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    split_bam = subparsers.add_parser('split_bam_by_barcode')
    split_bam.set_defaults(which='split_bam_by_barcode')
    split_bam.add_argument('--infile', required=True)
    split_bam.add_argument('--outdir', required=True)

    identify_normal_cells = subparsers.add_parser('identify_normal_cells')
    identify_normal_cells.set_defaults(which='identify_normal_cells')
    identify_normal_cells.add_argument('--reads_data', required=True)
    identify_normal_cells.add_argument('--metrics_data', required=True)
    identify_normal_cells.add_argument('--output_yaml', required=True)
    identify_normal_cells.add_argument('--reference_name', required=True)
    identify_normal_cells.add_argument('--min_reads', default=500000)
    identify_normal_cells.add_argument('--min_quality', default=0.85)
    identify_normal_cells.add_argument('--allowed_aneuploidy_score', default=0)
    identify_normal_cells.add_argument('--relative_aneuploidy_threshold', default=0.05)
    identify_normal_cells.add_argument('--ploidy_threshold', default=2.5)

    separate_normal_and_tumour_cells = subparsers.add_parser('separate_normal_and_tumour_cells')
    separate_normal_and_tumour_cells.set_defaults(which='separate_normal_and_tumour_cells')
    separate_normal_and_tumour_cells.add_argument('--infile', required=True)
    separate_normal_and_tumour_cells.add_argument('--normal_output', required=True)
    separate_normal_and_tumour_cells.add_argument('--tumour_output', required=True)
    separate_normal_and_tumour_cells.add_argument('--normal_cells_yaml', required=True)

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'split_bam_by_barcode':
        split_bam_by_barcode(
            args['infile'], args['outdir']
        )
    elif args['which'] == 'identify_normal_cells':
        identify_normal_cells.identify_normal_cells(
            args['reads_data'],
            args['metrics_data'],
            args['output_yaml'],
            args['reference_name'],
            min_reads=args['min_reads'],
            min_quality=args['min_quality'],
            allowed_aneuploidy_score=args['allowed_aneuploidy_score'],
            relative_aneuploidy_threshold=args['relative_aneuploidy_threshold'],
            ploidy_threshold=args['ploidy_threshold']
        )
    elif args['which'] == 'separate_normal_and_tumour_cells':
        separate_normal_and_tumour_cells(
            args['infile'],
            args['normal_output'],
            args['tumour_output'],
            args['normal_cells_yaml'],
        )
    else:
        raise Exception()
