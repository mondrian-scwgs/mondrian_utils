import click

import mondrianutils.normalizer


@click.group()
def cli():
    pass


@cli.command()
@click.option('--reads_data', required=True)
@click.option('--metrics_data', required=True)
@click.option('--normal_copy', required=True)
@click.option('--output_yaml', required=True)
@click.option('--output_csv', required=True)
@click.option('--output_plot', required=True)
@click.option('--min_reads', default=500000)
@click.option('--min_quality', default=0.85)
@click.option('--aneuploidy_score_threshold', default=0.005)
@click.option('--ploidy_threshold', default=2.5)
def identify_normal_cells(
        reads_data, metrics_data, normal_copy, output_yaml, output_csv, output_plot,
        min_reads, min_quality, aneuploidy_score_threshold, ploidy_threshold
):
    mondrianutils.normalizer.identify_normal_cells.identify_normal_cells(
        reads_data,
        metrics_data,
        normal_copy,
        output_yaml,
        output_csv,
        output_plot,
        min_reads=min_reads,
        min_quality=min_quality,
        aneuploidy_score_threshold=aneuploidy_score_threshold,
        ploidy_threshold=ploidy_threshold
    )


@cli.command()
@click.option('--infile', required=True)
@click.option('--normal_output', required=True)
@click.option('--tumour_output', required=True)
@click.option('--normal_cells_yaml', required=True)
def separate_normal_and_tumour_cells(
        infile, normal_output, tumour_output, normal_cells_yaml
):
    mondrianutils.normalizer.separate_normal_and_tumour_cells(
        infile,
        normal_output,
        tumour_output,
        normal_cells_yaml,
    )


@cli.command()
@click.option('--normal_bam', nargs=2)
@click.option('--tumour_bam', nargs=2)
@click.option('--heatmap', multiple=True)
@click.option('--normal_cells_yaml')
@click.option('--metadata_input')
@click.option('--metadata_output')
def separate_tumour_and_normal_metadata(
        normal_bam, tumour_bam, heatmap, normal_cells_yaml, metadata_input,
        metadata_output
):
    mondrianutils.normalizer.separate_tumour_and_normal_metadata(
        normal_bam, tumour_bam, heatmap,
        normal_cells_yaml, metadata_input,
        metadata_output
    )
