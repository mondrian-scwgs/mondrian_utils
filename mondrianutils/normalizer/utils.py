import copy
import os

import argparse
import pysam
import yaml
from mondrianutils import __version__
from mondrianutils import helpers
from mondrianutils.normalizer import heatmap
from mondrianutils.normalizer import identify_normal_cells


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


def update_header_comments(header, cells):
    if not isinstance(header, dict):
        newheader = copy.deepcopy(header.to_dict())
    else:
        newheader = copy.deepcopy(header)

    newheader['CO'] = [f'CB:{v}' for v in cells]

    return newheader


def update_normal_readgroup(header):
    rg_map = {}

    for readgroup in header['RG']:
        old_sample_id = readgroup['SM']
        new_sample_id = old_sample_id + '-N'

        readgroup['SM'] = new_sample_id

        old_readgroup = readgroup['ID']
        new_readgroup = readgroup['ID'].replace(old_sample_id, new_sample_id)

        readgroup['ID'] = new_readgroup

        rg_map[old_readgroup] = new_readgroup

    return header, rg_map


def update_readgroup_in_read(read, readgroup_mapping):
    assert read.get_tag('RG') in readgroup_mapping

    read.set_tag('RG', readgroup_mapping[read.get_tag('RG')])
    return read


def separate_normal_and_tumour_cells(infile, normal_output, tumour_output, normal_cells_yaml):
    with open(normal_cells_yaml, 'rt') as reader:
        normal_cells = set(yaml.safe_load(reader)['cells'])

    reader = pysam.AlignmentFile(infile, "rb")

    tumour_cells = [v for v in get_cells(reader) if v not in normal_cells]

    normal_header = update_header_comments(reader.header, normal_cells)
    normal_header, readgroup_mapping = update_normal_readgroup(normal_header)
    tumour_header = update_header_comments(reader.header, tumour_cells)

    normal_writer = pysam.AlignmentFile(normal_output, "wb", header=normal_header)
    tumour_writer = pysam.AlignmentFile(tumour_output, "wb", header=tumour_header)

    for pileupobj in reader.fetch(until_eof=True):
        cell_id = pileupobj.get_tag('CB')

        if cell_id in normal_cells:
            updated_read = update_readgroup_in_read(pileupobj, readgroup_mapping)
            normal_writer.write(updated_read)
        else:
            tumour_writer.write(pileupobj)


def separate_tumour_and_normal_metadata(
        tumour_bam, normal_bam, heatmap, normal_cells_yaml,
        metadata_input, metadata_output
):
    with open(metadata_input, 'rt') as reader:
        data = yaml.safe_load(reader)

    out_data = dict()
    out_data['meta'] = dict(
        type='hmmcopy',
        version=__version__,
        lane_ids=data['meta']['lane_ids'],
        sample_ids=data['meta']['sample_ids'],
        library_ids=data['meta']['library_ids'],
        cell_ids=data['meta']['cell_ids'],
    )

    files = {
        os.path.basename(normal_cells_yaml): {
            'result_type': 'yaml',
            'auxiliary': helpers.get_auxiliary_files(normal_cells_yaml)
        }
    }

    if normal_bam:
        files[os.path.basename(normal_bam[0])] = {
            'result_type': 'merged_cells_bam', 'filtering': 'passed',
            'auxiliary': helpers.get_auxiliary_files(normal_bam[0])
        }
        files[os.path.basename(normal_bam[1])] = {
            'result_type': 'merged_cells_bam', 'filtering': 'passed',
            'auxiliary': helpers.get_auxiliary_files(normal_bam[1])
        }

    if tumour_bam:
        files[os.path.basename(tumour_bam[0])] = {
            'result_type': 'merged_cells_bam', 'filtering': 'passed',
            'auxiliary': helpers.get_auxiliary_files(tumour_bam[0])
        }
        files[os.path.basename(tumour_bam[1])] = {
            'result_type': 'merged_cells_bam', 'filtering': 'passed',
            'auxiliary': helpers.get_auxiliary_files(tumour_bam[1])
        }

    for heatmap_file in heatmap:
        files[os.path.basename(heatmap_file)] = {
            'result_type': 'hmmcopy_heatmap_plots',
            'auxiliary': helpers.get_auxiliary_files(heatmap_file)
        }

    out_data['files'] = files

    with open(metadata_output, 'wt') as writer:
        yaml.dump(out_data, writer, default_flow_style=False)


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    identify_normal_cells = subparsers.add_parser('identify_normal_cells')
    identify_normal_cells.set_defaults(which='identify_normal_cells')
    identify_normal_cells.add_argument('--reads_data', required=True)
    identify_normal_cells.add_argument('--metrics_data', required=True)
    identify_normal_cells.add_argument('--output_yaml', required=True)
    identify_normal_cells.add_argument('--output_csv', required=True)
    identify_normal_cells.add_argument('--reference_name', required=True)
    identify_normal_cells.add_argument('--blacklist_file', required=True)
    identify_normal_cells.add_argument('--min_reads', default=500000)
    identify_normal_cells.add_argument('--min_quality', default=0.85)
    identify_normal_cells.add_argument('--allowed_aneuploidy_score', default=0)
    identify_normal_cells.add_argument('--relative_aneuploidy_threshold', default=0.005)
    identify_normal_cells.add_argument('--ploidy_threshold', default=2.5)

    separate_normal_and_tumour_cells = subparsers.add_parser('separate_normal_and_tumour_cells')
    separate_normal_and_tumour_cells.set_defaults(which='separate_normal_and_tumour_cells')
    separate_normal_and_tumour_cells.add_argument('--infile', required=True)
    separate_normal_and_tumour_cells.add_argument('--normal_output', required=True)
    separate_normal_and_tumour_cells.add_argument('--tumour_output', required=True)
    separate_normal_and_tumour_cells.add_argument('--normal_cells_yaml', required=True)

    aneuploidy_heatmap = subparsers.add_parser('aneuploidy_heatmap')
    aneuploidy_heatmap.set_defaults(which='aneuploidy_heatmap')
    aneuploidy_heatmap.add_argument('--metrics', required=True)
    aneuploidy_heatmap.add_argument('--reads', required=True)
    aneuploidy_heatmap.add_argument('--output', required=True)
    aneuploidy_heatmap.add_argument('--aneuploidy_score', default=0.005)

    separate_tumour_and_normal_metadata = subparsers.add_parser('separate_tumour_and_normal_metadata')
    separate_tumour_and_normal_metadata.set_defaults(which='separate_tumour_and_normal_metadata')
    separate_tumour_and_normal_metadata.add_argument(
        '--normal_bam', nargs=2
    )
    separate_tumour_and_normal_metadata.add_argument(
        '--tumour_bam', nargs=2
    )
    separate_tumour_and_normal_metadata.add_argument(
        '--heatmap', nargs='*'
    )
    separate_tumour_and_normal_metadata.add_argument(
        '--normal_cells_yaml'
    )
    separate_tumour_and_normal_metadata.add_argument(
        '--metadata_input'
    )
    separate_tumour_and_normal_metadata.add_argument(
        '--metadata_output'
    )

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'identify_normal_cells':
        identify_normal_cells.identify_normal_cells(
            args['reads_data'],
            args['metrics_data'],
            args['output_yaml'],
            args['output_csv'],
            args['reference_name'],
            args['blacklist_file'],
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
    elif args['which'] == 'aneuploidy_heatmap':
        heatmap.aneuploidy_heatmap(
            args['reads'],
            args['metrics'],
            args['output'],
            aneuploidy_score=args['aneuploidy_score']
        )
    elif args['which'] == 'separate_tumour_and_normal_metadata':
        separate_tumour_and_normal_metadata(
            args['normal_bam'], args['tumour_bam'], args['heatmap'],
            args['normal_cells_yaml'], args['metadata_input'],
            args['metadata_output']
        )
    else:
        raise Exception()
