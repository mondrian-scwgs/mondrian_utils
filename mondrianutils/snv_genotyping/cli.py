import click

import mondrianutils.snv_genotyping


@click.group()
def cli():
    pass


@cli.command()
@click.option('--bam', required=True)
@click.option('--output', required=True)
@click.option('--targets_vcf', required=True)
@click.option('--cell_barcodes')
@click.option('--interval')
@click.option('--count_duplicates', default=False)
@click.option('--sparse', is_flag=True, default=False)
@click.option('--ignore_untagged_reads', is_flag=True, default=False)
@click.option('--min_mqual', default=20)
@click.option('--skip_header', is_flag=True, default=False)
def snv_genotyper(
        bam, output, targets_vcf, cell_barcodes, interval, count_duplicates, sparse,
        ignore_untagged_reads, min_mqual, skip_header
):
    with mondrianutils.snv_genotyping.SnvGenotyper(
            bam, targets_vcf, output, cell_barcodes=cell_barcodes,
            interval=interval, count_duplicates=count_duplicates,
            sparse=sparse, min_mqual=min_mqual, ignore_untagged_reads=ignore_untagged_reads,
            skip_header=skip_header
    ) as genotyper:
        genotyper.genotyping()


@cli.command()
@click.option('--bamfile')
@click.option('--output')
def generate_cell_barcodes(
        bamfile, output
):
    mondrianutils.snv_genotyping.generate_cell_barcodes_file(bamfile, output)


@cli.command()
@click.option('--outputs', nargs=2)
@click.option('--vartrix_outputs', nargs=6)
@click.option('--pysam_genotyper', nargs=2)
@click.option('--metadata_input')
@click.option('--metadata_output')
def generate_metadata(
        outputs, vartrix_outputs, pysam_genotyper_outputs, metadata_input, metadata_output
):
    mondrianutils.snv_genotyping.generate_metadata(
        outputs, vartrix_outputs, pysam_genotyper_outputs, metadata_input, metadata_output
    )


@cli.command()
@click.option('--barcodes', multiple=True, required=True)
@click.option('--variants', multiple=True, required=True)
@click.option('--ref_matrices', multiple=True, required=True)
@click.option('--alt_matrices', multiple=True, required=True)
@click.option('--vcf_files', multiple=True, required=True)
@click.option('--parsed_output', required=True)
@click.option('--tempdir', required=True)
def merge_vartrix(
        barcodes, variants, ref_matrices, alt_matrices, vcf_files, parsed_output, tempdir
):
    mondrianutils.snv_genotyping.merge_vartrix(
        barcodes, variants, ref_matrices, alt_matrices, vcf_files,
        parsed_output, tempdir
    )


@cli.command()
@click.option('--barcode', required=True)
@click.option('--variant', required=True)
@click.option('--ref_matrix', required=True)
@click.option('--alt_matrix', required=True)
@click.option('--vcf_file', required=True)
@click.option('--parsed_output', required=True)
@click.option('--tempdir', required=True)
@click.option('--skip_header', is_flag=True, default=False)
def parse_vartrix(
        barcode, variant, ref_matrix, alt_matrix, vcf_file, parsed_output, tempdir, skip_header
):
    mondrianutils.snv_genotyping.parse_vartrix(
        barcode, variant, ref_matrix, alt_matrix, vcf_file,
        parsed_output, tempdir, skip_header=skip_header
    )


@cli.command()
@click.option('--barcodes', required=True)
@click.option('--variants', required=True)
@click.option('--ref_matrix', required=True)
@click.option('--alt_matrix', required=True)
@click.option('--parsed_data', required=True)
@click.option('--tempdir', required=True)
def regenerate_vartrix_format(
        barcodes, variants, ref_matrix, alt_matrix, parsed_data, tempdir
):
    mondrianutils.snv_genotyping.regenerate_vartrix_format(
        barcodes, variants, ref_matrix, alt_matrix,
        parsed_data, tempdir
    )


if __name__ == '__main__':
    cli()
