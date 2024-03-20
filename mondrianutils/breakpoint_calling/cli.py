import click

import mondrianutils.breakpoint_calling

@click.group()
def cli():
    pass


@cli.command()
@click.option('--destruct', required=True)
@click.option('--lumpy', required=True)
@click.option('--gridss', required=True)
@click.option('--svaba', required=True)
@click.option('--sample_id', required=True)
@click.option('--consensus_output', required=True)
@click.option('--tempdir', required=True)
@click.option('--region')
@click.option('--blacklist_bed')
def breakpoint_consensus(destruct, lumpy, gridss, svaba, sample_id, consensus_output, tempdir, region, blacklist_bed):
    mondrianutils.breakpoint_calling.consensus(
        destruct, lumpy, svaba, gridss, consensus_output,
        sample_id, tempdir, region=region,
        blacklist_bed=blacklist_bed
    )


@cli.command()
@click.option('--infile', required=True)
@click.option('--outfile', required=True)
@click.option('--reference', required=True)
@click.option('--sample_id', required=True)
def breakpoint_destruct_csv_to_vcf(infile, outfile, reference, sample_id):
    mondrianutils.breakpoint_calling.destruct_csv_to_vcf(
        infile, outfile, reference, sample_id
    )


@cli.command()
@click.option('--reads', required=True)
@click.option('--output', required=True)
def breakpoint_destruct_extract_cell_counts(reads, output):
    mondrianutils.breakpoint_calling.get_counts(
        reads, output
    )


@cli.command()
@click.option('--files')
@click.option('--metadata_yaml_files', multiple=True)
@click.option('--samples', multiple=True)
@click.option('--metadata_output')
def breakpoint_generate_metadata(files, metadata_yaml_files, samples, metadata_output):
    mondrianutils.breakpoint_calling.generate_metadata(
        files, metadata_yaml_files, samples, metadata_output
    )

@cli.command()
@click.option('--destruct', multiple=True)
@click.option('--consensus', multiple=True)
@click.option('--lumpy')
@click.option('--svaba')
@click.option('--gridss')
@click.option('--metadata_input')
@click.option('--metadata_output')
def breakpoint_generate_per_sample_metadata(
        destruct, consensus, lumpy, svaba, gridss, metadata_input, metadata_output
):
    mondrianutils.breakpoint_calling.generate_per_sample_metadata(
        destruct, consensus, lumpy, svaba, gridss, metadata_input, metadata_output
    )


if __name__ == '__main__':
    cli()
