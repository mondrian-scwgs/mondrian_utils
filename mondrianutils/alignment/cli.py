import click
import mondrianutils.alignment


@click.group()
def cli():
    pass


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
    mondrianutils.alignment.merge_cells_by_type(
        infiles, reference, cell_ids, metrics, control_outfile, contaminated_outfile, pass_outfile, tempdir,
        ncores
    )


@cli.command()
@click.option('--metrics', nargs=2, help='metrics files')
@click.option('--gc_metrics', nargs=2, help='GC metrics files')
@click.option('--bam', nargs=2, help='BAM files')
@click.option('--control', nargs=2, help='control files')
@click.option('--contaminated', nargs=2, help='contaminated files')
@click.option('--tarfile', help='tarfile')
@click.option('--metadata_input', help='metadata input file')
@click.option('--metadata_output', help='metadata output file')
def generate_metadata(
        metrics,
        gc_metrics,
        bam,
        control,
        contaminated,
        tarfile,
        metadata_input,
        metadata_output,
):
    mondrianutils.alignment.generate_metadata(
        bam, control, contaminated, metrics, gc_metrics, tarfile, metadata_input, metadata_output
    )


@cli.command()
@click.option('--meta_yaml', required=True, help='Path to the metadata YAML file')
@click.option('--input_data_json', required=True, help='Path to the input data JSON file')
def input_validation(meta_yaml, input_data_json):
    mondrianutils.alignment.input_validation(meta_yaml, input_data_json)


@cli.command()
@click.option('--fastq_pairs', multiple=True, help='Comma-separated list of FASTQ files')
@click.option('--metadata_yaml', help='Path to the metadata YAML file')
@click.option('--reference', help='Path to the reference file')
@click.option('--supplementary_references', multiple=True, help='Path to supplementary references')
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
@click.option('--tar_output', help='Path to the TAR output file')
@click.option('--num_threads', default=1, help='Number of threads')
@click.option('--run_fastqc', is_flag=True, help='Run FastQC')
@click.option('--run_tss_enrichment/--no_tss_enrichment', default=True, help='Run TSS enrichment analysis (default: enabled)')
def alignment(
        fastq_pairs, metadata_yaml, reference, supplementary_references,
        tempdir, adapter1, adapter2, cell_id, wgs_metrics_mqual,
        wgs_metrics_bqual, wgs_metrics_count_unpaired, bam_output, metrics_output,
        metrics_gc_output, tar_output, num_threads, run_fastqc, run_tss_enrichment
):
    mondrianutils.alignment.alignment(
        fastq_pairs, metadata_yaml,
        reference, supplementary_references,
        tempdir, adapter1, adapter2, cell_id, wgs_metrics_mqual,
        wgs_metrics_bqual, wgs_metrics_count_unpaired,
        bam_output, metrics_output, metrics_gc_output,
        tar_output, num_threads, run_fastqc=run_fastqc,
        run_tss_enrichment=run_tss_enrichment
    )


if __name__ == '__main__':
    cli()
