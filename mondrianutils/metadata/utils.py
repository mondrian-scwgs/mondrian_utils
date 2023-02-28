import os

import argparse
import mondrianutils.helpers as helpers
import yaml
from mondrianutils import __version__


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
        os.path.basename(heatmap): {
            'result_type': 'hmmcopy_heatmap_plots',
            'auxiliary': helpers.get_auxiliary_files(heatmap)
        },
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

    out_data['files'] = files

    with open(metadata_output, 'wt') as writer:
        yaml.dump(out_data, writer, default_flow_style=False)


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    separate_tumour_and_normal = subparsers.add_parser('separate_tumour_and_normal')
    separate_tumour_and_normal.set_defaults(which='separate_tumour_and_normal')
    separate_tumour_and_normal.add_argument(
        '--normal_bam', nargs=2
    )
    separate_tumour_and_normal.add_argument(
        '--tumour_bam', nargs=2
    )
    separate_tumour_and_normal.add_argument(
        '--heatmap'
    )
    separate_tumour_and_normal.add_argument(
        '--normal_cells_yaml'
    )
    separate_tumour_and_normal.add_argument(
        '--metadata_input'
    )
    separate_tumour_and_normal.add_argument(
        '--metadata_output'
    )

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'separate_tumour_and_normal':
        separate_tumour_and_normal_metadata(
            args['normal_bam'], args['tumour_bam'], args['heatmap'],
            args['normal_cells_yaml'], args['metadata_input'],
            args['metadata_output']
        )
    else:
        raise Exception()
