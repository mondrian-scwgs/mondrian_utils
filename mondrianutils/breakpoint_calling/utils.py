import os
import json
import yaml

from mondrianutils import helpers
from mondrianutils import __version__


def infer_type(files):
    with open(files, 'rt') as reader:
        files = json.load(reader)

    filetypes = sorted(set([v['left'] for v in files]))

    # more than one wf
    if 'svaba_vcf' in filetypes and 'breakpoint_consensus' in filetypes:
        return 'breakpoint_calling'
    elif 'breakpoint_consensus' in filetypes:
        return 'breakpoint_consensus'
    elif 'lumpy_vcf' in filetypes:
        return 'breakpoint_lumpy'
    elif 'svaba_vcf' in filetypes:
        return 'breakpoint_svaba'
    elif 'gridss_vcf' in filetypes:
        return 'breakpoint_gridss'
    elif 'destruct_calls' in filetypes:
        return 'breakpoint_destruct'
    else:
        raise Exception()


def generate_metadata(
        files, metadata_yaml_files, samples, metadata_output
):
    wf_type = infer_type(files)
    data = helpers.metadata_helper(files, metadata_yaml_files, samples, wf_type)

    with open(metadata_output, 'wt') as writer:
        yaml.dump(data, writer, default_flow_style=False)


def generate_per_sample_metadata(
        destruct, consensus, lumpy, svaba, gridss, metadata_input, metadata_output
):
    with open(metadata_input, 'rt') as reader:
        data = yaml.safe_load(reader)

    out_data = dict()
    out_data['meta'] = dict(
        type='haplotype_count',
        version=__version__,
        lane_ids=data['meta']['lane_ids'],
        sample_ids=data['meta']['sample_ids'],
        library_ids=data['meta']['library_ids'],
        cell_ids=data['meta']['cell_ids'],
    )

    files = {
        os.path.basename(lumpy): {
            'result_type': 'breakpoint_lumpy',
            'auxiliary': helpers.get_auxiliary_files(lumpy)
        },
        os.path.basename(svaba): {
            'result_type': 'breakpoint_svaba',
            'auxiliary': helpers.get_auxiliary_files(svaba)
        },
        os.path.basename(gridss): {
            'result_type': 'breakpoint_gridss',
            'auxiliary': helpers.get_auxiliary_files(gridss)
        }
    }

    for filepath in consensus:
        files[os.path.basename(filepath)] = {
            'result_type': 'breakpoint_consensus', 'auxiliary': helpers.get_auxiliary_files(filepath)
        }
    for filepath in destruct:
        files[os.path.basename(filepath)] = {
            'result_type': 'breakpoint_destruct', 'auxiliary': helpers.get_auxiliary_files(filepath)
        }
    out_data['files'] = files
    with open(metadata_output, 'wt') as writer:
        yaml.dump(out_data, writer, default_flow_style=False)
