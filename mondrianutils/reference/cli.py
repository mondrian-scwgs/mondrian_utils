import click

import mondrianutils.reference

@click.group()
def cli():
    pass


@cli.command()
@click.option('--rmsk_file', required=True)
@click.option('--chrom_map', required=True)
@click.option('--repeats', required=True)
@click.option('--satellites', required=True)
def repeats(rmsk_file, chrom_map, repeats, satellites):
    mondrianutils.reference.build_repeats(rmsk_file, chrom_map, repeats, satellites)


@cli.command()
@click.option('--reference_fai', required=True)
@click.option('--data_dir', required=True)
@click.option('--snp_positions', required=True)
def snp_positions_grch37(reference_fai, data_dir, snp_positions):
    mondrianutils.reference.create_snp_positions_grch37(reference_fai, data_dir, snp_positions)


@cli.command()
@click.option('--data_dir', required=True)
@click.option('--snp_positions', required=True)
@click.option('--chromosomes', multiple=True, default=None)
def snp_positions_grch38(data_dir, snp_positions, chromosomes):
    mondrianutils.reference.create_snp_positions_grch38(data_dir, snp_positions, chromosomes=chromosomes)


@cli.command()
@click.option('--reference', required=True, type=click.Path(exists=True))
@click.option('--output', required=True, type=click.Path())
@click.option('--chromosomes', multiple=True, required=True)
@click.option('--interval_size', type=int, default=1e6)
def get_intervals(reference, output, chromosomes, interval_size):
    mondrianutils.reference.get_intervals(reference, output, chromosomes, interval_size)


if __name__ == '__main__':
    cli()
