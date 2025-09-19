import os
import yaml
import mondrianutils.helpers as helpers
from mondrianutils import __version__


def generate_metadata(
        bam, control, contaminated, metrics, gc_metrics,
        reads, params, segments, heatmap, qc_html,
        alignment_tarfile, hmmcopy_tarfile,
        metadata_input, metadata_output
):
    with open(metadata_input, 'rt') as reader:
        data = yaml.safe_load(reader)

    lane_data = data['meta']['lanes']

    samples = set()
    libraries = set()
    cells = []
    for cell in data['meta']['cells']:
        cells.append(cell)
        samples.add(data['meta']['cells'][cell]['sample_id'])
        libraries.add(data['meta']['cells'][cell]['library_id'])

    data = dict()
    data['files'] = {
        os.path.basename(metrics[0]): {
            'result_type': 'qc_metrics',
            'auxiliary': helpers.get_auxiliary_files(metrics[0])
        },
        os.path.basename(metrics[1]): {
            'result_type': 'qc_metrics',
            'auxiliary': helpers.get_auxiliary_files(metrics[1])
        },
        os.path.basename(gc_metrics[0]): {
            'result_type': 'alignment_gc_metrics',
            'auxiliary': helpers.get_auxiliary_files(gc_metrics[0])
        },
        os.path.basename(gc_metrics[1]): {
            'result_type': 'alignment_gc_metrics',
            'auxiliary': helpers.get_auxiliary_files(gc_metrics[1])
        },
        os.path.basename(bam[0]): {
            'result_type': 'merged_cells_bam', 'filtering': 'passed',
            'auxiliary': helpers.get_auxiliary_files(bam[0])
        },
        os.path.basename(bam[1]): {
            'result_type': 'merged_cells_bam', 'filtering': 'passed',
            'auxiliary': helpers.get_auxiliary_files(bam[1])
        },
        os.path.basename(control[0]): {
            'result_type': 'merged_cells_bam', 'filtering': 'control',
            'auxiliary': helpers.get_auxiliary_files(control[0])
        },
        os.path.basename(control[1]): {
            'result_type': 'merged_cells_bam', 'filtering': 'control',
            'auxiliary': helpers.get_auxiliary_files(control[1])
        },
        os.path.basename(contaminated[0]): {
            'result_type': 'merged_cells_bam', 'filtering': 'contaminated',
            'auxiliary': helpers.get_auxiliary_files(contaminated[0])
        },
        os.path.basename(contaminated[1]): {
            'result_type': 'merged_cells_bam', 'filtering': 'contaminated',
            'auxiliary': helpers.get_auxiliary_files(contaminated[1])
        },
        os.path.basename(heatmap): {
            'result_type': 'hmmcopy_heatmap_plots',
            'auxiliary': helpers.get_auxiliary_files(heatmap)
        },
        os.path.basename(qc_html): {
            'result_type': 'qc_report_html',
            'auxiliary': helpers.get_auxiliary_files(qc_html)
        },
        os.path.basename(reads[0]): {
            'result_type': 'hmmcopy_reads',
            'auxiliary': helpers.get_auxiliary_files(reads[0])
        },
        os.path.basename(reads[1]): {
            'result_type': 'hmmcopy_reads',
            'auxiliary': helpers.get_auxiliary_files(reads[1])
        },
        os.path.basename(params[0]): {
            'result_type': 'hmmcopy_params',
            'auxiliary': helpers.get_auxiliary_files(params[0])
        },
        os.path.basename(params[1]): {
            'result_type': 'hmmcopy_params',
            'auxiliary': helpers.get_auxiliary_files(params[1])
        },
        os.path.basename(segments[0]): {
            'result_type': 'hmmcopy_segments',
            'auxiliary': helpers.get_auxiliary_files(segments[0])
        },
        os.path.basename(segments[1]): {
            'result_type': 'hmmcopy_segments',
            'auxiliary': helpers.get_auxiliary_files(segments[1])
        },
        os.path.basename(hmmcopy_tarfile): {
            'result_type': 'hmmcopy_metrics_tar',
            'auxiliary': helpers.get_auxiliary_files(hmmcopy_tarfile)
        },
        os.path.basename(alignment_tarfile): {
            'result_type': 'alignment_metrics_tar',
            'auxiliary': helpers.get_auxiliary_files(alignment_tarfile)
        },

    }

    data['meta'] = {
        'type': 'alignment',
        'version': __version__,
        'sample_ids': sorted(samples),
        'library_ids': sorted(libraries),
        'cell_ids': sorted(cells),
        'lane_ids': lane_data
    }

    with open(metadata_output, 'wt') as writer:
        yaml.dump(data, writer, default_flow_style=False)

