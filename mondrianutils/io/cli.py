import click
import mondrianutils.io


@click.group()
def cli():
    pass


@cli.command()
@click.option('--infile', required=True, help='Path to the input BAM file')
@click.option('--outdir', required=True, help='Path to the output directory')
@click.option('--tempdir', help='Path to the temporary directory')
@click.option('--chromosomes', multiple=True, help='List of chromosomes')
@click.option('--ncores', default=8, type=int, help='Number of cores')
def split_bam_by_barcode(infile, outdir, tempdir, chromosomes, ncores):
    mondrianutils.io.split_bam_by_barcode(infile, outdir, tempdir, chromosomes, ncores)


@cli.command()
@click.option('--bam', required=True, help='Path to the input BAM file')
@click.option('--output', required=True, help='Path to the output file')
@click.option('--tempdir', help='Path to the temporary directory')
@click.option('--chromosomes', default=[str(v) for v in range(1, 23)] + ['X', 'Y'], multiple=True,
              help='List of chromosomes')
@click.option('--binsize', default=500000, type=int, help='Size of bins')
@click.option('--mapping_quality', default=20, type=int, help='Mapping quality threshold')
@click.option('--ncores', default=8, type=int, help='Number of cores')
def overlapping_fraction_per_bin(bam, output, tempdir, chromosomes, binsize, mapping_quality, ncores):
    mondrianutils.io.overlapping_fraction_per_bin(
        bam, output, tempdir,
        chromosomes=chromosomes,
        binsize=binsize,
        mapping_quality=mapping_quality,
        ncores=ncores
    )


@cli.command()
@click.option('--infile', required=True, help='Path to the input CSV file')
@click.option('--dtypes', required=True, help='Path to the file containing data types')
@click.option('--outfile', required=True, help='Path to the output CSV file')
def rewrite_csv(infile, dtypes, outfile):
    mondrianutils.io.rewrite_csv(infile, outfile, dtypes)


@cli.command()
@click.option('--infiles', multiple=True, required=True, help='List of input PDF files')
@click.option('--outfile', required=True, help='Path to the output PDF file')
def merge_pdfs(infiles, outfile):
    mondrianutils.io.merge_pdfs(infiles, outfile)


@cli.command()
@click.option('--infiles', multiple=True, required=True, help='List of input PDF files')
@click.option('--outfile', required=True, help='Path to the output PDF file')
def merge_pdfs_scaled(infiles, outfile):
    mondrianutils.io.merge_pdfs_with_scaling(infiles, outfile)


@cli.command()
@click.option('--infile', required=True, help='Path to the input VCF file')
@click.option('--outdir', required=True, help='Directory for output VCF files')
@click.option('--num_splits', type=int, help='Number of splits for split_vcf command')
@click.option('--num_lines', type=int, help='Number of splits for split_vcf command')
def split_vcf(infile, outdir, num_splits=None, num_lines=None):
    if num_splits is not None and num_lines is not None:
        raise Exception("specify num_splits or num_lines, but not both")

    if num_lines is not None:
        mondrianutils.io.split_vcf_by_lines(infile, outdir, num_lines)
    else:
        mondrianutils.io.split_vcf_into_numsplits(infile, outdir, num_splits)


@cli.command()
@click.option('--infile', required=True, help='Path to the input VCF file')
@click.option('--outdir', required=True, help='Directory for output VCF files')
def split_vcf_by_chrom(infile, outdir):
    mondrianutils.io.split_vcf_by_chrom(infile, outdir)


@cli.command()
@click.option('--infile', required=True, help='Path to the input VCF file')
@click.option('--outfile', required=True, help='Path to the output VCF file')
@click.option('--include_ref_alt', is_flag=True, default=False,
              help='Include ref/alt information in duplicates removal')
def remove_duplicates(infile, outfile, include_ref_alt):
    mondrianutils.io.remove_duplicates(
        infile, outfile, include_ref_alt=include_ref_alt
    )


@cli.command()
@click.option('--infile', required=True, help='Path to the input VCF file')
@click.option('--outfile', required=True, help='Path to the output VCF file')
@click.option('--exclusion_blacklist', help='Path to the exclusion blacklist file')
def exclude_blacklist(infile, outfile, exclusion_blacklist):
    mondrianutils.io.exclude_blacklist(infile, outfile, exclusion_blacklist)


@cli.command()
@click.option('--infiles', multiple=True, required=True, help='List of input VCF files')
@click.option('--outfile', required=True, help='Path to the output VCF file')
def merge_vcfs(infiles, outfile):
    mondrianutils.io.merge_vcfs(infiles, outfile)


if __name__ == '__main__':
    cli()
