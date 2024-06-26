import os

import pysam
import yaml
from mondrianutils import helpers
from mondrianutils import __version__


def generate_metadata(
        outputs, vartrix_outputs, metadata_input, metadata_output
):
    data = dict()
    data['files'] = {}

    if outputs is not None:
        data['files'][os.path.basename(outputs[0])] = {
            'result_type': 'pysam_genotyping_counts', 'auxiliary': False
        }
        data['files'][os.path.basename(outputs[1])] = {
            'result_type': 'pysam_genotyping_counts', 'auxiliary': True
        }
    if vartrix_outputs is not None:
        for filepath in vartrix_outputs:
            data['files'][os.path.basename(filepath)] = {
                'result_type': 'vartrix_genotyping_counts', 'auxiliary': helpers.get_auxiliary_files(filepath)
            }

    with open(metadata_input, 'rt') as reader:
        meta = yaml.safe_load(reader)
        del meta['meta']['type']
        del meta['meta']['version']

    meta_dict = {
        'type': 'snv_genotyping',
        'version': __version__
    }

    data['meta'] = {**meta_dict, **meta['meta']}

    with open(metadata_output, 'wt') as writer:
        yaml.dump(data, writer, default_flow_style=False)


def generate_cell_barcodes_file(bamfile, output):
    bamfile = pysam.AlignmentFile(bamfile, 'rb')
    header = bamfile.header

    cells = []
    for line in str(header).split('\n'):
        if not line.startswith("@CO\tCB:"):
            continue
        line = line.strip().split()
        cb = line[1]
        cell = cb.split(':')[1]
        cells.append(cell)

    with open(output, 'wt') as writer:
        for cell in cells:
            writer.write(cell + '\n')
