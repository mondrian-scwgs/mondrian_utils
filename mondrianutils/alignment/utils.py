import json
import os
import subprocess
from collections import defaultdict

import csverve.api as csverve
import mondrianutils.helpers as helpers
import pysam
import yaml
from mondrianutils import __version__
from mondrianutils.dtypes.alignment import dtypes


class MultipleSamplesPerRun(Exception):
    pass


class MissingField(Exception):
    pass


def _check_sample_id_uniqueness(meta_data):
    non_control_samples = [v['sample_id'] for k, v in meta_data['cells'].items() if v['is_control'] == False]
    non_control_samples = sorted(set(non_control_samples))

    # allow 0 sample ids in case all cells are control
    if not len(non_control_samples) <= 1:
        raise MultipleSamplesPerRun(
            f'only one sample id expected in non control cells, found {non_control_samples}'
        )


def _check_metadata_required_field(meta_data, field_name):
    try:
        [v[field_name] for k, v in meta_data['cells'].items()]
    except KeyError:
        raise MissingField(f'{field_name} is required for each cell in meta section of metadata input yaml')


def _check_lanes_and_flowcells(meta_data, input_data):
    lane_data = defaultdict(set)

    for val in input_data:
        for lane in val['lanes']:
            lane_data[lane['flowcell_id']].add(lane['lane_id'])

    for flowcell, lanes in lane_data.items():
        if flowcell not in meta_data['lanes']:
            raise MissingField(
                f'missing flowcell {flowcell} in metadata yaml'
            )
        for lane_id in lanes:
            if lane_id not in meta_data['lanes'][flowcell]:
                raise MissingField(
                    f'missing lane {lane_id} for flowcell {flowcell} in metadata yaml'
                )


def input_validation(meta_yaml, input_json):
    with open(meta_yaml, 'rt') as reader:
        meta_data = yaml.safe_load(reader)
        meta_data = meta_data['meta']

    with open(input_json, 'rt') as reader:
        input_data = json.load(reader)

    _check_metadata_required_field(meta_data, 'is_control')
    _check_metadata_required_field(meta_data, 'library_id')
    _check_metadata_required_field(meta_data, 'sample_id')
    # these are required in hmmcopy
    _check_metadata_required_field(meta_data, 'pick_met')
    _check_metadata_required_field(meta_data, 'condition')

    _check_sample_id_uniqueness(meta_data)

    _check_lanes_and_flowcells(meta_data, input_data)

    for flowcell_id, all_lane_data in meta_data['lanes'].items():
        for lane_id, lane_data in all_lane_data.items():
            if 'sequencing_centre' not in lane_data:
                raise MissingField(
                    f'sequencing centre missing for flowcell {flowcell_id} lane {lane_id}'
                )


def get_new_header(cells, bamfile, new_header):
    subprocess.run(['samtools', 'view', '-H', bamfile, '-o', new_header])
    with open(new_header, 'at') as header:
        for cell in cells:
            header.write('@CO\tCB:{}\n'.format(cell))


def reheader(infile, new_header, outfile):
    subprocess.run(
        ['picard', 'ReplaceSamHeader', 'I={}'.format(infile),
         'HEADER={}'.format(new_header), 'O={}'.format(outfile)
         ]
    )


def get_pass_files(infiles, cell_ids, metrics):
    metrics = csverve.read_csv(metrics)
    assert set(cell_ids) == set(list(metrics['cell_id']))

    cells_to_skip = set(list(metrics[metrics['is_contaminated']]['cell_id']))
    infiles = {cell: infile for cell, infile in zip(cell_ids, infiles) if cell not in cells_to_skip}

    cells_to_skip = set(list(metrics[metrics['is_control']]['cell_id']))
    infiles = {cell: infile for cell, infile in infiles.items() if cell not in cells_to_skip}

    return infiles


def get_control_files(infiles, cell_ids, metrics):
    metrics = csverve.read_csv(metrics)
    assert set(cell_ids) == set(list(metrics['cell_id']))
    control_cells = set(list(metrics[metrics['is_control'] == True]['cell_id']))
    infiles = {cell: infile for cell, infile in zip(cell_ids, infiles) if cell in control_cells}
    return infiles


def get_contaminated_files(infiles, cell_ids, metrics):
    metrics = csverve.read_csv(metrics)
    assert set(cell_ids) == set(list(metrics['cell_id']))

    cells_to_skip = set(list(metrics[metrics['is_control']]['cell_id']))
    infiles = {cell: infile for cell, infile in zip(cell_ids, infiles) if cell not in cells_to_skip}

    contaminated_cells = set(list(metrics[metrics['is_contaminated']]['cell_id']))
    infiles = {cell: infile for cell, infile in infiles.items() if cell in contaminated_cells}

    return infiles


def samtools_index(infile):
    cmd = ['samtools', 'index', infile]
    helpers.run_cmd(cmd)


def igvtools_count(infile, reference):
    cmd = ['igvtools', 'count', infile, infile + '.tdf', reference]
    helpers.run_cmd(cmd)


def merge_bams(infiles, tempdir, ncores, outfile, reference, empty_bam_content):
    if len(infiles.values()) == 0:
        pysam.AlignmentFile(outfile, "wb", header=empty_bam_content).close()
    else:
        helpers.merge_bams(list(infiles.values()), outfile, tempdir, ncores)

    samtools_index(outfile)
    igvtools_count(outfile, reference)


def get_bam_header(bam):
    infile = pysam.AlignmentFile(bam, "rb")

    header = infile.header


    return header


def merge_cells_by_type(
        infiles, reference, cell_ids, metrics,
        control_outfile, contaminated_outfile, pass_outfile,
        tempdir, ncores
):
    header = get_bam_header(infiles[0])
    # controls
    control_bams = get_control_files(infiles, cell_ids, metrics)
    control_tempdir = os.path.join(tempdir, 'control')
    helpers.makedirs(control_tempdir)
    merge_bams(control_bams, control_tempdir, ncores, control_outfile, reference, header)

    # contaminated
    contaminated_bams = get_contaminated_files(infiles, cell_ids, metrics)
    contaminated_tempdir = os.path.join(tempdir, 'contaminated')
    helpers.makedirs(contaminated_tempdir)
    merge_bams(contaminated_bams, contaminated_tempdir, ncores, contaminated_outfile, reference, header)

    # pass
    pass_bams = get_pass_files(infiles, cell_ids, metrics)
    pass_tempdir = os.path.join(tempdir, 'pass')
    helpers.makedirs(pass_tempdir)
    merge_bams(pass_bams, pass_tempdir, ncores, pass_outfile, reference, header)


def tag_bam_with_cellid(infile, outfile, cell_id):
    infile = pysam.AlignmentFile(infile, "rb")
    outfile = pysam.AlignmentFile(outfile, "wb", template=infile)

    iter = infile.fetch(until_eof=True)
    for read in iter:
        read.set_tag('CB', cell_id, replace=False)
        outfile.write(read)
    infile.close()
    outfile.close()


def add_contamination_status(
        infile, outfile,
        reference, threshold=0.05
):
    def get_col_data(df, organism):
        return df['fastqscreen_{}'.format(organism)] - df['fastqscreen_{}_multihit'.format(organism)]

    data = csverve.read_csv(infile)

    data = data.set_index('cell_id', drop=False)

    organisms = [v for v in data.columns.values if v.startswith('fastqscreen_')]
    organisms = sorted(set([v.split('_')[1] for v in organisms]))
    organisms = [v for v in organisms if v not in ['nohit', 'total']]

    if reference not in organisms:
        raise Exception("Could not find the fastq screen counts")

    alts = [col for col in organisms if not col == reference]

    data['is_contaminated'] = False

    for altcol in alts:
        perc_alt = get_col_data(data, altcol) / data['fastqscreen_total_reads']
        data.loc[perc_alt > threshold, 'is_contaminated'] = True

    col_type = dtypes()['metrics']['is_contaminated']

    data['is_contaminated'] = data['is_contaminated'].astype(col_type)
    csverve.write_dataframe_to_csv_and_yaml(
        data, outfile, dtypes(fastqscreen_genomes=organisms)['metrics']
    )


def add_metadata(metrics, metadata_yaml, output):
    df = csverve.read_csv(metrics)

    metadata = yaml.safe_load(open(metadata_yaml, 'rt'))

    cells = set(df['cell_id'])

    metadata['meta']['cells'] = {k: v for k, v in metadata['meta']['cells'].items() if k in cells}

    for cellid, cell_info in metadata['meta']['cells'].items():
        for colname, val in cell_info.items():
            df.loc[df['cell_id'] == cellid, colname] = val

    organisms = [v for v in df.columns.values if v.startswith('fastqscreen_')]
    organisms = sorted(set([v.split('_')[1] for v in organisms]))
    organisms = [v for v in organisms if v not in ['nohit', 'total']]

    csverve.write_dataframe_to_csv_and_yaml(
        df, output,
        dtypes=dtypes(fastqscreen_genomes=organisms)['metrics']
    )


def generate_metadata(
        bam, control, contaminated, metrics, gc_metrics,
        tarfile, metadata_input, metadata_output
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
            'result_type': 'alignment_metrics',
            'auxiliary': helpers.get_auxiliary_files(metrics[0])
        },
        os.path.basename(metrics[1]): {
            'result_type': 'alignment_metrics',
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
        os.path.basename(tarfile): {
            'result_type': 'alignment_metrics_plots',
            'auxiliary': helpers.get_auxiliary_files(tarfile)
        }
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


def supplementary_reference_cmdline(jsonfile):
    refs = []
    with open(jsonfile, 'rt') as reader:
        data = json.load(reader)
        for entry in data:
            refs.append(f'{entry["genome_name"]},{entry["genome_version"]},{entry["reference"]}')

        refs = [f'--supplementary_references {v}' for v in refs]
        print(" ".join(refs))


def fastqs_cmdline(jsonfile):
    pairs = []
    with open(jsonfile, 'rt') as reader:
        data = json.load(reader)
        for entry in data:
            pairs.append(f"{entry['lane_id']},{entry['flowcell_id']},{entry['fastq1']},{entry['fastq2']}")
    pairs = [f'--fastq_pairs {v}' for v in pairs]
    print(" ".join(pairs))
