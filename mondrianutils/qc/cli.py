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

