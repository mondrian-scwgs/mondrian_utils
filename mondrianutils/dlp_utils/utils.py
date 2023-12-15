import click

from mondrianutils.dlp_utils.dlp_bams_to_mondrian_bam import dlp_bams_to_mondrian_bam


@click.group()
def cli():
    pass

@click.command()
@click.option('--dlp_bam_dir', required=True)
@click.option('--output', required=True)
@click.option('--tempdir', required=True)
@click.option('--cores', default=8)
def dlp_bams_to_mondrian_bam_cmd(dlp_bam_dir, output, tempdir, cores):
    dlp_bams_to_mondrian_bam(dlp_bam_dir, output, tempdir, ncores=cores)

if __name__ == '__main__':
    cli()

