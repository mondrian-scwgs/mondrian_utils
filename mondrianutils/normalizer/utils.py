import copy
import os

import click
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


