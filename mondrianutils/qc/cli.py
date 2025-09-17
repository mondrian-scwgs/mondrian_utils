import click
import mondrianutils.qc


@click.group()
def cli():
    pass


@cli.command()
@click.option('--bam', nargs=2, help='BAM files')
@click.option('--control', nargs=2, help='control files')
@click.option('--contaminated', nargs=2, help='contaminated files')
@click.option('--gc_metrics', nargs=2, help='GC metrics files')
@click.option('--alignment_tar')
@click.option('--metrics', nargs=2)
@click.option('--params', nargs=2)
@click.option('--reads', nargs=2)
@click.option('--segments', nargs=2)
@click.option('--hmmcopy_tar')
@click.option('--heatmap')
@click.option('--qc_report_html')
@click.option('--metadata_input')
@click.option('--metadata_output')
def generate_metadata(
        bam, control, contaminated, gc_metrics, alignment_tar,
        metrics, params, reads, segments, hmmcopy_tar,
        heatmap, qc_report_html,
        metadata_input, metadata_output
):
    mondrianutils.qc.generate_metadata(
        bam, control, contaminated, metrics, gc_metrics,
        reads, params, segments, heatmap, qc_report_html,
        alignment_tar, hmmcopy_tar,
        metadata_input, metadata_output
    )


@cli.command()
@click.option('--stats_file', required=True, help='Input BAM stats file')
@click.option('--output_file', required=True, help='Output pickle file')
def parse_bamstats(stats_file, output_file):
    """Parse BAM statistics from samtools stats output."""
    import pickle
    result = mondrianutils.qc.parse_bamstat(stats_file)
    with open(output_file, 'wb') as f:
        pickle.dump(result, f)


@cli.command()
@click.option('--results_dir', required=True, help='Kraken2 results directory')
@click.option('--cell_id', required=True, help='Cell ID')
@click.option('--output_table', required=True, help='Output table file')
@click.option('--output_human', required=True, help='Output human reads file')
@click.option('--output_nonhuman', required=True, help='Output non-human reads file')
def parse_kraken_output(results_dir, cell_id, output_table, output_human, output_nonhuman):
    """Parse Kraken2 output and separate human vs non-human reads."""
    mondrianutils.qc.parse_kraken_output(
        results_dir, cell_id, output_table, output_human, output_nonhuman
    )


@cli.command()
@click.option('--pipeline_outputs_dir', required=True, help='Pipeline output directory')
@click.option('--hmmcopy_metrics_filename', required=True, help='HMMcopy metrics file')
@click.option('--library_id', required=True, help='Library ID')
@click.option('--output_dir', required=True, help='Output directory')
@click.option('--ncbi_taxonomy_database', default='/data1/shahs3/users/myersm2/repos/contamination/kraken_db/ncbi_taxonomy/taxa.sqlite', help='Path to ete3.NCBITaxa sqlite database file')
@click.option('--min_percent_aggregate', default=0.0, help='Minimum percent of non-human reads required to include a taxon in data table')
@click.option('--min_percent_show', default=2.0, help='Minimum percent of non-human reads required to show a taxon in bar graph')
@click.option('--min_num_taxa_condense', default=25, help='Minimum number of taxa required to perform tree-cutting procedure')
def generate_contamination_table_figures(
        pipeline_outputs_dir,
        hmmcopy_metrics_filename,
        library_id,
        output_dir,
        ncbi_taxonomy_database,
        min_percent_aggregate,
        min_percent_show,
        min_num_taxa_condense
):
    """Generate contamination analysis tables and figures from pipeline outputs."""
    mondrianutils.qc.generate_contamination_table_figures(
        pipeline_outputs_dir,
        hmmcopy_metrics_filename,
        library_id,
        output_dir,
        ncbi_taxonomy_database,
        min_percent_aggregate,
        min_percent_show,
        min_num_taxa_condense
    )
