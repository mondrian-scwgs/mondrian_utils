import os
import click
import mondrianutils.hmmcopy
import csverve
import mondrianutils.helpers as helpers
from mondrianutils.dtypes.hmmcopy import dtypes as hmmcopy_dtypes


@click.group()
def cli():
    pass


@cli.command()
@click.option('--infile', required=True)
@click.option('--outdir', required=True)
@click.option('--tempdir')
@click.option('--chromosomes', multiple=True, default=[str(v) for v in range(1, 23)] + ['X', 'Y'],
              help='specify target chromosomes')
@click.option('-w', '--window_size', type=int, default=500000, help='specify bin size')
@click.option('-m', '--mapping_quality_threshold', type=int, default=20,
              help='threshold for the mapping quality, reads with quality lower than threshold will be ignored')
@click.option('--exclude_list', default=None, help='regions to skip')
@click.option('--tag_name', default='CB', type=str, help='name of the cell barcode tag')
@click.option('--discover_cells', is_flag=True, default=False, type=bool, help='discover cell barcodes instead of reading from header')
@click.option('--ignore_missing_tags', is_flag=True, default=False, type=bool, help='ignore missing tags')
@click.option('--ncores', default=12, type=int, help='number of cores to use')
def readcounter(infile, outdir, tempdir, chromosomes, window_size, mapping_quality_threshold, exclude_list, tag_name, discover_cells, ignore_missing_tags, ncores):
    mondrianutils.hmmcopy.readcounter(
        infile,
        outdir,
        chromosomes,
        exclude_list=exclude_list,
        ncores=ncores,
        mapping_quality_threshold=mapping_quality_threshold,
        window_size=window_size,
        tag_name=tag_name,
        discover_cells=discover_cells,
        ignore_missing_tags=ignore_missing_tags,
    )


@cli.command()
@click.option('--readcount_wig', required=True)
@click.option('--gc_wig_file', required=True)
@click.option('--map_wig_file', required=True)
@click.option('--alignment_metrics', required=True)
@click.option('--metrics', required=True)
@click.option('--params', required=True)
@click.option('--reads', required=True)
@click.option('--segments', required=True)
@click.option('--output_tarball', required=True)
@click.option('--reference', required=True)
@click.option('--segments_output', required=True)
@click.option('--bias_output', required=True)
@click.option('--cell_id', required=True)
@click.option('--tempdir', required=True)
@click.option('--quality_classifier_training_data', required=True)
@click.option('--quality_classifier_model')
@click.option('--map_cutoff', default=0.9, type=float)
def run_hmmcopy(
        readcount_wig, gc_wig_file, map_wig_file, alignment_metrics,
        metrics, params, reads, segments, output_tarball,
        reference, segments_output, bias_output, cell_id, tempdir,
        quality_classifier_training_data, quality_classifier_model,
        map_cutoff
):
    mondrianutils.hmmcopy.complete_hmmcopy(
        readcount_wig, gc_wig_file, map_wig_file, alignment_metrics,
        metrics, params, reads, segments, output_tarball,
        reference, segments_output, bias_output,
        cell_id, tempdir,
        quality_classifier_training_data,
        quality_classifier_model=quality_classifier_model,
        mappability_cutoff=map_cutoff
    )


@cli.command()
@click.option('--bam_file', required=True)
@click.option('--gc_wig_file', required=True)
@click.option('--map_wig_file', required=True)
@click.option('--alignment_metrics', required=True)
@click.option('--exclude_list', required=True)
@click.option('--chromosomes', multiple=True, default=[str(v) for v in range(1, 23)] + ['X', 'Y'],
              help='specify target chromosomes')
@click.option('--metrics', required=True)
@click.option('--params', required=True)
@click.option('--reads', required=True)
@click.option('--segments', required=True)
@click.option('--output_tarball', required=True)
@click.option('--reference', required=True)
@click.option('--segments_output', required=True)
@click.option('--bias_output', required=True)
@click.option('--tempdir', required=True)
@click.option('--quality_classifier_training_data', required=True)
@click.option('--quality_classifier_model')
@click.option('--map_cutoff', default=0.9, type=float)
@click.option('--binsize', default=500000, type=int)
@click.option('--mapping_quality_threshold', default=20, type=int)
def run_cell_hmmcopy(
        bam_file, gc_wig_file, map_wig_file, alignment_metrics, exclude_list, chromosomes,
        metrics, params, reads, segments, output_tarball,
        reference, segments_output, bias_output, tempdir,
        quality_classifier_training_data, quality_classifier_model,
        map_cutoff, binsize, mapping_quality_threshold
):
    mondrianutils.hmmcopy.cell_hmmcopy(
        bam_file, gc_wig_file, map_wig_file, alignment_metrics, exclude_list,
        chromosomes, metrics, params, reads, segments, output_tarball,
        reference, segments_output, bias_output,
        tempdir, quality_classifier_training_data,
        quality_classifier_model=quality_classifier_model,
        mappability_cutoff=map_cutoff,
        binsize=binsize,
        mapping_quality_threshold=mapping_quality_threshold
    )


@cli.command()
@click.option('--reads')
@click.option('--metrics')
@click.option('--chromosomes', default=[str(v) for v in range(1, 23)] + ['X', 'Y'], multiple=True)
@click.option('--output')
@click.option('--sidebar_column', default='pick_met')
@click.option('--disable_clustering', default=False, is_flag=True, help='Disable clustering')
def heatmap(reads, metrics, chromosomes, output, sidebar_column, disable_clustering):
    mondrianutils.hmmcopy.plot_heatmap(
        reads, metrics, chromosomes, output,
        sidebar_column=sidebar_column, disable_clustering=disable_clustering
    )


@cli.command()
@click.option('--segs_pdf', multiple=True)
@click.option('--segs_samples', multiple=True)
@click.option('--metrics')
@click.option('--pass_output')
@click.option('--fail_output')
@click.option('--tempdir')
def create_segs_tar(segs_pdf, segs_samples, metrics, pass_output, fail_output, tempdir):
    mondrianutils.hmmcopy.create_segs_tar(
        segs_pdf, metrics,
        pass_output, fail_output,
        tempdir, segs_samples=segs_samples
    )


@cli.command()
@click.option('--tempdir')
@click.option('--html')
@click.option('--metrics')
@click.option('--gc_metrics')
def generate_html_report(tempdir, html, metrics, gc_metrics):
    mondrianutils.hmmcopy.generate_html_report(
        tempdir, html, metrics, gc_metrics
    )


@cli.command()
@click.option('--reads')
@click.option('--metrics')
@click.option('--output')
@click.option('--chromosomes', default=[str(v) for v in range(1, 23)] + ['X', 'Y'], multiple=True)
def add_clustering_order(reads, metrics, output, chromosomes):
    mondrianutils.hmmcopy.add_clustering_order(reads, metrics, output, chromosomes=chromosomes)


@cli.command()
@click.option('--reads')
@click.option('--alignment_metrics')
@click.option('--hmmcopy_metrics')
@click.option('--tempdir')
@click.option('--output')
def cell_cycle_classifier(reads, alignment_metrics, hmmcopy_metrics, tempdir, output):
    mondrianutils.hmmcopy.cell_cycle_classifier(
        reads, alignment_metrics, hmmcopy_metrics, tempdir, output
    )


@cli.command()
@click.option('--metrics', nargs=2)
@click.option('--params', nargs=2)
@click.option('--reads', nargs=2)
@click.option('--segments', nargs=2)
@click.option('--segments_tar_pass')
@click.option('--segments_tar_fail')
@click.option('--heatmap')
@click.option('--metadata_input')
@click.option('--metadata_output')
def generate_metadata(metrics, params, reads, segments, segments_tar_pass, segments_tar_fail, heatmap,
                      metadata_input, metadata_output):
    mondrianutils.hmmcopy.generate_metadata(
        metrics, params, reads, segments,
        segments_tar_pass, segments_tar_fail,
        heatmap, metadata_input, metadata_output
    )


if __name__ == '__main__':
    cli()
