import json
import os
import subprocess
from collections import defaultdict

import click
import csverve.api as csverve
import mondrianutils.helpers as helpers
import pysam
import yaml
from mondrianutils import __version__
from mondrianutils.alignment.collect_gc_metrics import collect_gc_metrics
from mondrianutils.alignment.collect_metrics import collect_metrics
from mondrianutils.alignment.complete_alignment import alignment
from mondrianutils.alignment.coverage_metrics import get_coverage_metrics
from mondrianutils.dtypes.alignment import dtypes
from mondrianutils.alignment.fastqscreen import merge_fastq_screen_counts
from mondrianutils.alignment.fastqscreen import organism_filter
from mondrianutils.alignment.trim_galore import trim_galore


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


def get_cell_id_from_bam(infile):
    infile = pysam.AlignmentFile(infile, "rb")

    iter = infile.fetch(until_eof=True)
    for read in iter:
        return read.get_tag('CB')


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


def merge_cells(infiles, tempdir, ncores, outfile, reference, empty_bam_content):
    if len(infiles.values()) == 0:
        pysam.AlignmentFile(outfile, "wb", header=empty_bam_content).close()
    else:
        final_merge_output = os.path.join(tempdir, 'merged_all.bam')
        helpers.merge_bams(list(infiles.values()), final_merge_output, tempdir, ncores)

        new_header = os.path.join(tempdir, 'header.sam')
        get_new_header(infiles.keys(), final_merge_output, new_header)

        reheader(final_merge_output, new_header, outfile)

    samtools_index(outfile)
    igvtools_count(outfile, reference)


def get_bam_header(bam):
    infile = pysam.AlignmentFile(bam, "rb")

    header = infile.header

    if 'CO' in header:
        del header['CO']

    return header


def generate_bams(
        infiles, reference, cell_ids, metrics,
        control_outfile, contaminated_outfile, pass_outfile,
        tempdir, ncores
):
    header = get_bam_header(infiles[0])
    # controls
    control_bams = get_control_files(infiles, cell_ids, metrics)
    control_tempdir = os.path.join(tempdir, 'control')
    helpers.makedirs(control_tempdir)
    merge_cells(control_bams, control_tempdir, ncores, control_outfile, reference, header)

    # contaminated
    contaminated_bams = get_contaminated_files(infiles, cell_ids, metrics)
    contaminated_tempdir = os.path.join(tempdir, 'contaminated')
    helpers.makedirs(contaminated_tempdir)
    merge_cells(contaminated_bams, contaminated_tempdir, ncores, contaminated_outfile, reference, header)

    # pass
    pass_bams = get_pass_files(infiles, cell_ids, metrics)
    pass_tempdir = os.path.join(tempdir, 'pass')
    helpers.makedirs(pass_tempdir)
    merge_cells(pass_bams, pass_tempdir, ncores, pass_outfile, reference, header)


def tag_bam_with_cellid(infile, outfile, cell_id):
    infile = pysam.AlignmentFile(infile, "rb")
    outfile = pysam.AlignmentFile(outfile, "wb", template=infile)

    iter = infile.fetch(until_eof=True)
    for read in iter:
        read.set_tag('CB', cell_id, replace=False)
        outfile.write(read)
    infile.close()
    outfile.close()


def _get_col_data(df, organism):
    return df['fastqscreen_{}'.format(organism)] - df['fastqscreen_{}_multihit'.format(organism)]


def add_contamination_status(
        infile, outfile,
        reference, threshold=0.05
):
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
        perc_alt = _get_col_data(data, altcol) / data['fastqscreen_total_reads']
        data.loc[perc_alt > threshold, 'is_contaminated'] = True

    col_type = dtypes()['metrics']['is_contaminated']

    data['is_contaminated'] = data['is_contaminated'].astype(col_type)
    csverve.write_dataframe_to_csv_and_yaml(
        data, outfile, dtypes(fastqscreen_genomes=organisms)['metrics']
    )


def add_metadata(metrics, metadata_yaml, output):
    df = csverve.read_csv(metrics)

    metadata = yaml.safe_load(open(metadata_yaml, 'rt'))

    cells = metadata['meta']['cells'].keys()

    assert set(cells) == set(df['cell_id'])

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
        fastqscreen, tarfile, metadata_input, metadata_output
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
        os.path.basename(fastqscreen[0]): {
            'result_type': 'detailed_fastqscreen_breakdown',
            'auxiliary': helpers.get_auxiliary_files(fastqscreen[0])
        },
        os.path.basename(fastqscreen[1]): {
            'result_type': 'detailed_fastqscreen_breakdown',
            'auxiliary': helpers.get_auxiliary_files(fastqscreen[1])
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


def _json_file_parser(filepath):
    return json.load(open(filepath, 'rt'))


@click.group()
def cli():
    pass


@cli.command()
@click.option('--r1', help='specify R1 fastq')
@click.option('--r2', help='specify R2 fastq')
@click.option('--output_r1', help='specify output R1 fastq')
@click.option('--output_r2', help='specify output R2 fastq')
@click.option('--detailed_metrics', help='specify detailed metrics file')
@click.option('--summary_metrics', help='specify summary metrics file')
@click.option('--tempdir', help='specify temporary directory')
@click.option('--cell_id', help='specify cell ID')
@click.option('--human_reference', help='specify human reference fasta')
@click.option('--mouse_reference', help='specify mouse reference fasta')
@click.option('--salmon_reference', help='specify salmon reference fasta')
def fastqscreen_cmd(r1, r2, output_r1, output_r2, detailed_metrics, summary_metrics, tempdir, cell_id, human_reference,
                    mouse_reference, salmon_reference):
    organism_filter(r1, r2, output_r1, output_r2, detailed_metrics, summary_metrics, tempdir, cell_id, human_reference,
                    mouse_reference, salmon_reference)


@cli.command()
@click.option('--detailed_counts', multiple=True, help='detailed counts files')
@click.option('--summary_counts', multiple=True, help='summary counts files')
@click.option('--merged_detailed', help='merged detailed counts file')
@click.option('--merged_summary', help='merged summary counts file')
def merge_fastqscreen_counts_cmd(detailed_counts, summary_counts, merged_detailed, merged_summary):
    merge_fastq_screen_counts(detailed_counts, summary_counts, merged_detailed, merged_summary)


@cli.command()
@click.option('--wgs_metrics', help='path to WGS metrics file')
@click.option('--insert_metrics', help='path to insert metrics file')
@click.option('--flagstat', help='path to flagstat file')
@click.option('--markdups_metrics', help='path to markdups metrics file')
@click.option('--coverage_metrics', help='path to coverage metrics file')
@click.option('--output', help='output file')
@click.option('--cell_id', help='cell ID')
def collect_metrics_cmd(wgs_metrics, insert_metrics, flagstat, markdups_metrics, coverage_metrics, output, cell_id):
    collect_metrics(wgs_metrics, insert_metrics, flagstat, markdups_metrics, coverage_metrics, output, cell_id)


@cli.command()
@click.option('--infile', help='input file')
@click.option('--outfile', help='output file')
@click.option('--cell_id', help='cell ID')
def collect_gc_metrics_cmd(infile, outfile, cell_id):
    collect_gc_metrics(infile, outfile, cell_id)


@cli.command()
@click.option('--infile', help='input file')
@click.option('--outfile', help='output file')
@click.option('--cell_id', help='cell ID')
def tag_bam_with_cellid_cmd(infile, outfile, cell_id):
    tag_bam_with_cellid(infile, outfile, cell_id)


@cli.command()
@click.option('--infile', help='input file')
@click.option('--outfile', help='output file')
@click.option('--reference', help='reference file')
def add_contamination_status_cmd(infile, outfile, reference):
    add_contamination_status(infile, outfile, reference)


@cli.command()
@click.option('--infiles', multiple=True, help='input files')
@click.option('--cell_ids', multiple=True, help='cell IDs')
@click.option('--reference', help='reference file')
@click.option('--control_outfile', help='output file for control cells')
@click.option('--contaminated_outfile', help='output file for contaminated cells')
@click.option('--pass_outfile', help='output file for cells that pass')
@click.option('--metrics', help='metrics file')
@click.option('--tempdir', help='temporary directory')
@click.option('--ncores', type=int, help='number of cores')
def merge_cells_cmd(infiles, cell_ids, reference, control_outfile, contaminated_outfile, pass_outfile, metrics, tempdir,
                    ncores):
    generate_bams(infiles, reference, cell_ids, metrics, control_outfile, contaminated_outfile, pass_outfile, tempdir,
                  ncores)


# @cli.command()
# @click.option('--training_data', help='training data file')
# @click.option('--metrics', help='metrics file')
# @click.option('--output', help='output file')
# def classify_fastqscreen_cmd(training_data, metrics, output):
#     classify_fastqscreen(training_data, metrics, output)


@cli.command()
@click.option('--bamfile', help='BAM file')
@click.option('--output', help='output file')
def coverage_metrics_cmd(bamfile, output):
    get_coverage_metrics(bamfile, output)


@cli.command()
@click.option('--metrics', nargs=2, help='metrics files')
@click.option('--gc_metrics', nargs=2, help='GC metrics files')
@click.option('--bam', nargs=2, help='BAM files')
@click.option('--control', nargs=2, help='control files')
@click.option('--contaminated', nargs=2, help='contaminated files')
@click.option('--fastqscreen_detailed', nargs=2, help='fastqscreen detailed files')
@click.option('--tarfile', help='tarfile')
@click.option('--metadata_input', help='metadata input file')
@click.option('--metadata_output', help='metadata output file')
def generate_metadata_cmd(
        metrics,
        gc_metrics,
        bam,
        control,
        contaminated,
        fastqscreen_detailed,
        tarfile,
        metadata_input,
        metadata_output,
):
    generate_metadata(
        bam, control, contaminated, metrics, gc_metrics, fastqscreen_detailed, tarfile, metadata_input, metadata_output
    )


@cli.command()
@click.option('--metrics', help='metrics file')
@click.option('--metadata', help='metadata YAML file')
@click.option('--output', help='output file')
def add_metadata_cmd(metrics, metadata, output):
    add_metadata(metrics, metadata, output)


@cli.command()
@click.option('--r1', help='Input read 1 file')
@click.option('--r2', help='Input read 2 file')
@click.option('--output_r1', help='Output trimmed read 1 file')
@click.option('--output_r2', help='Output trimmed read 2 file')
@click.option('--adapter1', help='Adapter sequence for read 1')
@click.option('--adapter2', help='Adapter sequence for read 2')
@click.option('--tempdir', help='Temporary directory')
def trim_galore_cmd(r1, r2, output_r1, output_r2, adapter1, adapter2, tempdir):
    trim_galore(r1, r2, output_r1, output_r2, adapter1, adapter2, tempdir)


@cli.command()
@click.option('--meta_yaml', required=True, help='Path to the metadata YAML file')
@click.option('--input_data_json', required=True, help='Path to the input data JSON file')
def input_validation_cmd(meta_yaml, input_data_json):
    input_validation(meta_yaml, input_data_json)


@cli.command()
@click.option('--fastq_files', help='Comma-separated list of FASTQ files')
@click.option('--metadata_yaml', help='Path to the metadata YAML file')
@click.option('--reference', help='Path to the reference file')
@click.option('--reference_name', help='Reference name')
@click.option('--reference_version', help='Reference version')
@click.option('--supplementary_references_json', help='Path to supplementary references JSON file')
@click.option('--tempdir', help='Path to the temporary directory')
@click.option('--adapter1', help='Adapter sequence for read 1 trimming')
@click.option('--adapter2', help='Adapter sequence for read 2 trimming')
@click.option('--cell_id', help='Cell ID')
@click.option('--wgs_metrics_mqual', help='Path to the WGS metrics mqual file')
@click.option('--wgs_metrics_bqual', help='Path to the WGS metrics bqual file')
@click.option('--wgs_metrics_count_unpaired', help='Path to the WGS metrics count unpaired file')
@click.option('--bam_output', help='Path to the BAM output file')
@click.option('--metrics_output', help='Path to the metrics output file')
@click.option('--metrics_gc_output', help='Path to the GC metrics output file')
@click.option('--fastqscreen_detailed_output', help='Path to the FastQScreen detailed output file')
@click.option('--fastqscreen_summary_output', help='Path to the FastQScreen summary output file')
@click.option('--tar_output', help='Path to the TAR output file')
@click.option('--num_threads', default=1, help='Number of threads')
@click.option('--run_fastqc', is_flag=True, help='Run FastQC')
def alignment_cmd(
        fastq_files, metadata_yaml, reference, reference_name, reference_version,
        supplementary_references_json, tempdir, adapter1, adapter2, cell_id, wgs_metrics_mqual,
        wgs_metrics_bqual, wgs_metrics_count_unpaired, bam_output, metrics_output,
        metrics_gc_output, fastqscreen_detailed_output, fastqscreen_summary_output,
        tar_output, num_threads, run_fastqc
):
    alignment(
        fastq_files, metadata_yaml, reference,
        reference_name, reference_version, supplementary_references_json, tempdir,
        adapter1, adapter2, cell_id, wgs_metrics_mqual,
        wgs_metrics_bqual, wgs_metrics_count_unpaired,
        bam_output, metrics_output, metrics_gc_output,
        fastqscreen_detailed_output, fastqscreen_summary_output,
        tar_output, num_threads, run_fastqc=run_fastqc
    )


if __name__ == '__main__':
    cli()
