import click
import mondrianutils.alignment


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
def fastqscreen(r1, r2, output_r1, output_r2, detailed_metrics, summary_metrics, tempdir, cell_id, human_reference,
                mouse_reference, salmon_reference):
    mondrianutils.alignment.organism_filter(
        r1, r2, output_r1, output_r2, detailed_metrics, summary_metrics, tempdir, cell_id, human_reference,
        mouse_reference, salmon_reference
    )


@cli.command()
@click.option('--detailed_counts', multiple=True, help='detailed counts files')
@click.option('--summary_counts', multiple=True, help='summary counts files')
@click.option('--merged_detailed', help='merged detailed counts file')
@click.option('--merged_summary', help='merged summary counts file')
def merge_fastqscreen_counts(detailed_counts, summary_counts, merged_detailed, merged_summary):
    mondrianutils.alignment.merge_fastq_screen_counts(
        detailed_counts, summary_counts, merged_detailed, merged_summary
    )


@cli.command()
@click.option('--wgs_metrics', help='path to WGS metrics file')
@click.option('--insert_metrics', help='path to insert metrics file')
@click.option('--flagstat', help='path to flagstat file')
@click.option('--markdups_metrics', help='path to markdups metrics file')
@click.option('--coverage_metrics', help='path to coverage metrics file')
@click.option('--output', help='output file')
@click.option('--cell_id', help='cell ID')
def collect_metrics(wgs_metrics, insert_metrics, flagstat, markdups_metrics, coverage_metrics, output, cell_id):
    mondrianutils.alignment.collect_metrics(
        wgs_metrics, insert_metrics, flagstat, markdups_metrics, coverage_metrics, output, cell_id
    )


@cli.command()
@click.option('--infile', help='input file')
@click.option('--outfile', help='output file')
@click.option('--cell_id', help='cell ID')
def collect_gc_metrics(infile, outfile, cell_id):
    mondrianutils.alignment.collect_gc_metrics(
        infile, outfile, cell_id
    )


@cli.command()
@click.option('--infile', help='input file')
@click.option('--outfile', help='output file')
@click.option('--cell_id', help='cell ID')
def tag_bam_with_cellid(infile, outfile, cell_id):
    mondrianutils.alignment.tag_bam_with_cellid(
        infile, outfile, cell_id
    )


@cli.command()
@click.option('--infile', help='input file')
@click.option('--outfile', help='output file')
@click.option('--reference', help='reference file')
def add_contamination_status(infile, outfile, reference):
    mondrianutils.alignment.add_contamination_status(
        infile, outfile, reference
    )


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
def merge_cells(infiles, cell_ids, reference, control_outfile, contaminated_outfile, pass_outfile, metrics, tempdir,
                    ncores):
    mondrianutils.alignment.generate_bams(
        infiles, reference, cell_ids, metrics, control_outfile, contaminated_outfile, pass_outfile, tempdir,
        ncores
    )


@cli.command()
@click.option('--bamfile', help='BAM file')
@click.option('--output', help='output file')
def coverage_metrics(bamfile, output):
    mondrianutils.alignment.get_coverage_metrics(
        bamfile, output
    )


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
def generate_metadata(
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
    mondrianutils.alignment.generate_metadata(
        bam, control, contaminated, metrics, gc_metrics, fastqscreen_detailed, tarfile, metadata_input, metadata_output
    )


@cli.command()
@click.option('--metrics', help='metrics file')
@click.option('--metadata', help='metadata YAML file')
@click.option('--output', help='output file')
def add_metadata(metrics, metadata, output):
    mondrianutils.alignment.add_metadata(metrics, metadata, output)


@cli.command()
@click.option('--r1', help='Input read 1 file')
@click.option('--r2', help='Input read 2 file')
@click.option('--output_r1', help='Output trimmed read 1 file')
@click.option('--output_r2', help='Output trimmed read 2 file')
@click.option('--adapter1', help='Adapter sequence for read 1')
@click.option('--adapter2', help='Adapter sequence for read 2')
@click.option('--tempdir', help='Temporary directory')
def trim_galore(r1, r2, output_r1, output_r2, adapter1, adapter2, tempdir):
    mondrianutils.alignment.trim_galore(r1, r2, output_r1, output_r2, adapter1, adapter2, tempdir)


@cli.command()
@click.option('--meta_yaml', required=True, help='Path to the metadata YAML file')
@click.option('--input_data_json', required=True, help='Path to the input data JSON file')
def input_validation(meta_yaml, input_data_json):
    mondrianutils.alignment.input_validation(meta_yaml, input_data_json)


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
def alignment(
        fastq_files, metadata_yaml, reference, reference_name, reference_version,
        supplementary_references_json, tempdir, adapter1, adapter2, cell_id, wgs_metrics_mqual,
        wgs_metrics_bqual, wgs_metrics_count_unpaired, bam_output, metrics_output,
        metrics_gc_output, fastqscreen_detailed_output, fastqscreen_summary_output,
        tar_output, num_threads, run_fastqc
):
    mondrianutils.alignment.alignment(
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
