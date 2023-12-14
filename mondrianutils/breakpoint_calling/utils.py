import json
import yaml
import click

from mondrianutils import helpers
from .consensus import consensus
from mondrianutils.breakpoint_calling import destruct_csv_to_vcf
from mondrianutils.breakpoint_calling import destruct_extract_cell_counts


def infer_type(files):
    with open(files, 'rt') as reader:
        files = json.load(reader)

    filetypes = sorted(set([v['left'] for v in files]))

    # more than one wf
    if 'svaba_vcf' in filetypes and 'breakpoint_consensus' in filetypes:
        return 'breakpoint_calling'
    elif 'breakpoint_consensus' in filetypes:
        return 'breakpoint_consensus'
    elif 'lumpy_vcf' in filetypes:
        return 'breakpoint_lumpy'
    elif 'svaba_vcf' in filetypes:
        return 'breakpoint_svaba'
    elif 'gridss_vcf' in filetypes:
        return 'breakpoint_gridss'
    elif 'destruct_calls' in filetypes:
        return 'breakpoint_destruct'
    else:
        raise Exception()


def generate_metadata(
        files, metadata_yaml_files, samples, metadata_output
):
    wf_type = infer_type(files)
    data = helpers.metadata_helper(files, metadata_yaml_files, samples, wf_type)

    with open(metadata_output, 'wt') as writer:
        yaml.dump(data, writer, default_flow_style=False)


@click.group()
def cli():
    pass


@cli.command()
@click.option('--destruct', required=True)
@click.option('--lumpy', required=True)
@click.option('--gridss', required=True)
@click.option('--svaba', required=True)
@click.option('--sample_id', required=True)
@click.option('--consensus', required=True)
@click.option('--tempdir', required=True)
@click.option('--region')
@click.option('--blacklist_bed')
def breakpoint_consensus(destruct, lumpy, gridss, svaba, sample_id, consensus_output, tempdir, region, blacklist_bed):
    consensus(
        destruct, lumpy, svaba, gridss, consensus_output, sample_id, tempdir, region=region,
        blacklist_bed=blacklist_bed
    )


@cli.command()
@click.option('--infile', required=True)
@click.option('--outfile', required=True)
@click.option('--reference', required=True)
@click.option('--sample_id', required=True)
def destruct_csv_to_vcf_cmd(infile, outfile, reference, sample_id):
    destruct_csv_to_vcf.destruct_csv_to_vcf(infile, outfile, reference, sample_id)


@cli.command()
@click.option('--reads', required=True)
@click.option('--output', required=True)
def destruct_extract_cell_counts_cmd(reads, output):
    destruct_extract_cell_counts.get_counts(reads, output)


@cli.command()
@click.option('--files')
@click.option('--metadata_yaml_files', multiple=True)
@click.option('--samples', multiple=True)
@click.option('--metadata_output')
def generate_metadata_cmd(files, metadata_yaml_files, samples, metadata_output):
    generate_metadata(files, metadata_yaml_files, samples, metadata_output)


if __name__ == '__main__':
    cli()
