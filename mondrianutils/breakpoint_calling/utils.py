import json
import yaml

from mondrianutils import helpers

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


