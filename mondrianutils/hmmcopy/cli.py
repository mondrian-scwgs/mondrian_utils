import click
import mondrianutils.hmmcopy

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
@click.option('--ncores', default=12, type=int, help='number of cores to use')
def readcounter(infile, outdir, tempdir, chromosomes, window_size, mapping_quality_threshold, exclude_list, ncores):
    mondrianutils.hmmcopy.readcounter(
        infile,
        outdir,
        tempdir,
        chromosomes,
        exclude_list=exclude_list,
        ncores=ncores,
        mapping_quality_threshold=mapping_quality_threshold,
        window_size=window_size
    )


@cli.command()
@click.option('--readcount_wig', required=True)
@click.option('--gc_wig_file', required=True)
@click.option('--map_wig_file', required=True)
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
@click.option('--map_cutoff', default=0.9, type=float)
def run_hmmcopy(readcount_wig, gc_wig_file, map_wig_file, metrics, params, reads, segments, output_tarball,
                    reference, segments_output, bias_output, cell_id, tempdir, map_cutoff):
    mondrianutils.hmmcopy.complete_hmmcopy(
        readcount_wig, gc_wig_file, map_wig_file,
        metrics, params,
        reads, segments, output_tarball,
        reference, segments_output, bias_output,
        cell_id, tempdir,
        mappability_cutoff=map_cutoff
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
@click.option('--infile', required=True)
@click.option('--outfile', required=True)
def add_mappability(infile, outfile):
    mondrianutils.hmmcopy.add_mappability(infile, outfile)


@cli.command()
@click.option('--hmmcopy_metrics', required=True)
@click.option('--alignment_metrics', required=True)
@click.option('--training_data', required=True)
@click.option('--output', required=True)
@click.option('--tempdir', required=True)
def add_quality(hmmcopy_metrics, alignment_metrics, training_data, output, tempdir):
    mondrianutils.hmmcopy.add_quality(hmmcopy_metrics, alignment_metrics, tempdir, output,
                training_data)


@cli.command()
@click.option('--segs_pdf', multiple=True)
@click.option('--metrics')
@click.option('--pass_output')
@click.option('--fail_output')
@click.option('--tempdir')
def create_segs_tar(segs_pdf, metrics, pass_output, fail_output, tempdir):
    mondrianutils.hmmcopy.create_segs_tar(
        segs_pdf, metrics,
        pass_output, fail_output,
        tempdir
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
