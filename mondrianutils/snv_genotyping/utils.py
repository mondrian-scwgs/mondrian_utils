import os

import click
import pysam
import yaml
from mondrianutils import __version__
from mondrianutils.snv_genotyping.merge_vartrix import merge_vartrix
from mondrianutils.snv_genotyping.merge_vartrix import regenerate_vartrix_format
from mondrianutils.snv_genotyping.merge_vartrix import parse_vartrix
from mondrianutils.snv_genotyping.snv_genotyper import SnvGenotyper


def generate_metadata(
        outputs, vartrix_outputs, metadata_input, metadata_output
):
    data = dict()
    data['files'] = {
        os.path.basename(outputs[0]): {'result_type': 'pysam_genotyping_counts', 'auxiliary': False},
        os.path.basename(outputs[1]): {'result_type': 'pysam_genotyping_counts', 'auxiliary': True},
        os.path.basename(vartrix_outputs[0]): {'result_type': 'vartrix_genotyping_counts', 'auxiliary': False},
        os.path.basename(vartrix_outputs[1]): {'result_type': 'vartrix_genotyping_counts', 'auxiliary': True},
    }

    with open(metadata_input, 'rt') as reader:
        meta = yaml.safe_load(reader)
        del meta['meta']['type']
        del meta['meta']['version']

    meta_dict = {
        'type': 'snv_genotyping',
        'version': __version__
    }

    data['meta'] = {**meta_dict, **meta['meta']}

    with open(metadata_output, 'wt') as writer:
        yaml.dump(data, writer, default_flow_style=False)


def generate_cell_barcodes_file(bamfile, output):
    bamfile = pysam.AlignmentFile(bamfile, 'rb')
    header = bamfile.header

    cells = []
    for line in str(header).split('\n'):
        if not line.startswith("@CO\tCB:"):
            continue
        line = line.strip().split()
        cb = line[1]
        cell = cb.split(':')[1]
        cells.append(cell)

    with open(output, 'wt') as writer:
        for cell in cells:
            writer.write(cell + '\n')


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
def snv_genotyper_cmd(
        bam, output, targets_vcf, cell_barcodes, interval, count_duplicates, sparse,
        ignore_untagged_reads, min_mqual, skip_header
):
    with SnvGenotyper(
            bam, targets_vcf, output, cell_barcodes=cell_barcodes,
            interval=interval, count_duplicates=count_duplicates,
            sparse=sparse, min_mqual=min_mqual, ignore_untagged_reads=ignore_untagged_reads,
            skip_header=skip_header
    ) as genotyper:
        genotyper.genotyping()


@cli.command()
@click.option('--bamfile')
@click.option('--output')
def generate_cell_barcodes_cmd(
        bamfile, output
):
    generate_cell_barcodes_file(bamfile, output)


@cli.command()
@click.option('--outputs', nargs=2)
@click.option('--vartrix_outputs', nargs=6)
@click.option('--metadata_input')
@click.option('--metadata_output')
def generate_metadata_cmd(
        outputs, vartrix_outputs, metadata_input, metadata_output
):
    generate_metadata(
        outputs, vartrix_outputs, metadata_input, metadata_output
    )


@cli.command()
@click.option('--barcodes', multiple=True, required=True)
@click.option('--variants', multiple=True, required=True)
@click.option('--ref_matrices', multiple=True, required=True)
@click.option('--alt_matrices', multiple=True, required=True)
@click.option('--vcf_files', multiple=True, required=True)
@click.option('--parsed_output', required=True)
@click.option('--tempdir', required=True)
def merge_vartrix_cmd(
        barcodes, variants, ref_matrices, alt_matrices, vcf_files, parsed_output, tempdir
):
    merge_vartrix(
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
def parse_vartrix_cmd(
        barcode, variant, ref_matrix, alt_matrix, vcf_file, parsed_output, tempdir, skip_header
):
    parse_vartrix(
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
def regenerate_vartrix_format_cmd(
        barcodes, variants, ref_matrix, alt_matrix, parsed_data, tempdir
):
    regenerate_vartrix_format(
        barcodes, variants, ref_matrix, alt_matrix,
        parsed_data, tempdir
    )


if __name__ == '__main__':
    cli()
