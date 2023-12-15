import os

import click
import yaml
from mondrianutils import __version__
from mondrianutils.sv_genotyping.sv_genotyper import SvGenotyper


def generate_metadata(
        outputs, metadata_input, metadata_output
):
    data = dict()
    data['files'] = {
        os.path.basename(outputs[0]): {'result_type': 'pysam_genotyping_counts', 'auxiliary': False},
        os.path.basename(outputs[1]): {'result_type': 'pysam_genotyping_counts', 'auxiliary': True},
    }

    with open(metadata_input, 'rt') as reader:
        meta = yaml.safe_load(reader)
        del meta['meta']['type']
        del meta['meta']['version']

    meta_dict = {
        'type': 'sv_genotyping',
        'version': __version__
    }

    data['meta'] = {**meta_dict, **meta['meta']}

    with open(metadata_output, 'wt') as writer:
        yaml.dump(data, writer, default_flow_style=False)


@click.group()
def cli():
    pass


@cli.command()
@click.option('--bam', required=True)
@click.option('--destruct_reads', required=True)
@click.option('--destruct_table', required=True)
@click.option('--output')
def sv_genotyper_cmd(bam, destruct_reads, destruct_table, output):
    genotyper = SvGenotyper(
        bam, destruct_reads, destruct_table, output
    )
    genotyper.main()


@cli.command()
@click.option('--outputs', nargs=2)
@click.option('--metadata_input')
@click.option('--metadata_output')
def generate_metadata_cmd(outputs, metadata_input, metadata_output):
    generate_metadata(outputs, metadata_input, metadata_output)


if __name__ == '__main__':
    cli()
