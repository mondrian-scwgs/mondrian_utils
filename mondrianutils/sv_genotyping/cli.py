import click

import mondrianutils.sv_genotyping


@click.group()
def cli():
    pass


@cli.command()
@click.option('--bam', required=True)
@click.option('--destruct_reads', required=True)
@click.option('--destruct_table', required=True)
@click.option('--output')
def sv_genotyper(bam, destruct_reads, destruct_table, output):
    genotyper = mondrianutils.sv_genotyping.SvGenotyper(
        bam, destruct_reads, destruct_table, output
    )
    genotyper.main()


@cli.command()
@click.option('--outputs', nargs=2)
@click.option('--metadata_input')
@click.option('--metadata_output')
def generate_metadata(outputs, metadata_input, metadata_output):
    mondrianutils.sv_genotyping.generate_metadata(outputs, metadata_input, metadata_output)


if __name__ == '__main__':
    cli()
