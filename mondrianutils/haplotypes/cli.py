import os
import click
import mondrianutils.haplotypes
import remixt.analysis
import remixt.analysis.haplotype
import remixt.analysis.readcount
import remixt.analysis.segment
import remixt.seqdataio
import remixt.workflow
from mondrianutils import helpers


@click.group()
def cli():
    pass


@cli.command()
@click.option('--bam', required=True)
@click.option('--output', required=True)
@click.option('--snp_positions', required=True)
@click.option('--tempdir', required=True)
@click.option('--chromosomes', multiple=True, default=[str(v) for v in range(1, 23)] + ['X', 'Y'])
@click.option('--cell_id')
def extract_seqdata(bam, output, snp_positions, tempdir, chromosomes, cell_id):
    bam_max_fragment_length = remixt.config.get_param({}, 'bam_max_fragment_length')
    bam_max_soft_clipped = remixt.config.get_param({}, 'bam_max_soft_clipped')
    bam_check_proper_pair = remixt.config.get_param({}, 'bam_check_proper_pair')
    remixt.seqdataio.create_seqdata(
        output, bam, snp_positions,
        bam_max_fragment_length, bam_max_soft_clipped,
        bam_check_proper_pair, tempdir, chromosomes
    )
    if cell_id:
        mondrianutils.haplotypes.add_cell_id_to_seqdata(output, cell_id)


@cli.command()
@click.option('--reference_fai', required=True)
@click.option('--gap_table', required=True)
@click.option('--chromosomes', multiple=True, default=[str(v) for v in range(1, 23)] + ['X', 'Y'])
@click.option('--output', required=True)
@click.option('--tempdir', required=True)
def create_segments(reference_fai, gap_table, chromosomes, output, tempdir):
    chr_prefix = 'chr' if chromosomes[0].startswith('chr') else ''
    helpers.makedirs(tempdir)
    new_gap = os.path.join(tempdir, 'gap.txt.gz')
    with helpers.getFileHandle(gap_table, 'rt') as reader, helpers.getFileHandle(new_gap, 'wt') as writer:
        for line in reader:
            line = line.strip().split('\t')
            if chr_prefix == 'chr':
                line[1] = line[1] if line[1].startswith('chr') else 'chr' + line[1]
            else:
                line[1] = line[1][3:] if line[1].startswith('chr') else line[1]
            line = '\t'.join(line) + '\n'
            writer.write(line)

    config = {
        'genome_fai_template': reference_fai,
        'gap_table_template': new_gap,
        'chromosomes': chromosomes,
        'chr_name_prefix': chr_prefix
    }
    ref_data_dir = os.path.dirname(reference_fai)
    remixt.analysis.segment.create_segments(
        output, config, ref_data_dir
    )


@cli.command()
@click.option('--segments', required=True)
@click.option('--seqdata', required=True)
@click.option('--haplotypes', required=True)
@click.option('--output', required=True)
@click.option('--tempdir', required=True)
@click.option('--skip_header', is_flag=True, default=False)
def haplotype_allele_readcount(segments, seqdata, haplotypes, output, tempdir, skip_header):
    helpers.makedirs(tempdir)
    tempout = os.path.join(tempdir, 'temp.tsv')
    remixt.analysis.readcount.haplotype_allele_readcount(
        tempout, segments, seqdata, haplotypes, {}
    )
    mondrianutils.haplotypes.finalize_tsv(tempout, output, seqdata, skip_header=skip_header)


@cli.command()
@click.option('--files')
@click.option('--metadata_yaml_files', multiple=True)
@click.option('--samples', multiple=True)
@click.option('--metadata_output')
def generate_metadata(files, metadata_yaml_files, samples, metadata_output):
    mondrianutils.haplotypes.generate_metadata(
        files, metadata_yaml_files, samples,
        metadata_output
    )


@cli.command()
@click.option('--csv', required=True)
@click.option('--yaml', required=True)
@click.option('--metadata_yaml')
@click.option('--sample')
@click.option('--metadata_output')
def generate_infer_haps_metadata(csv, yaml, metadata_yaml, sample, metadata_output):
    mondrianutils.haplotypes.generate_infer_haps_metadata(
        csv, yaml, metadata_yaml, sample, metadata_output
    )

@cli.command()
@click.option('--csv', required=True)
@click.option('--yaml', required=True)
@click.option('--barcodes', required=True)
@click.option('--variants', required=True)
@click.option('--ref_counts', required=True)
@click.option('--alt_counts', required=True)
@click.option('--metadata_yaml')
@click.option('--sample')
@click.option('--metadata_output')
def generate_count_haps_metadata(
        csv, yaml, barcodes, variants, ref_counts, alt_counts,
        metadata_yaml, sample, metadata_output
):
    mondrianutils.haplotypes.generate_count_haps_metadata(
        csv, yaml, barcodes, variants, ref_counts, alt_counts,
        metadata_yaml, sample, metadata_output
    )


@cli.command()
@click.option('--input_bcf_file', required=True)
@click.option('--genetic_map', required=True)
@click.option('--regions_file', required=True)
@click.option('--chromosome', required=True)
@click.option('--tempdir', required=True)
@click.option('--output', required=True)
@click.option('--phased_chromosomes', multiple=True, default=tuple([f'chr{v}' for v in range(1, 23)] + ['chrX']))
@click.option('--is_female', is_flag=True, default=False)
@click.option('--phased_chromosome_x', default='chrX')
@click.option('--shapeit_num_samples', default=100, type=int)
@click.option('--shapeit_confidence_threshold', default=0.95, type=float)
def run_shapeit(input_bcf_file, genetic_map, regions_file, chromosome, tempdir, output, phased_chromosomes,
                is_female, phased_chromosome_x, shapeit_num_samples, shapeit_confidence_threshold):
    mondrianutils.haplotypes.run_shapeit4(
        input_bcf_file, genetic_map, regions_file, chromosome,
        tempdir, output,
        phased_chromosomes=phased_chromosomes,
        is_female=is_female,
        phased_chromx=phased_chromosome_x,
        shapeit_num_samples=shapeit_num_samples,
        shapeit_confidence_threshold=shapeit_confidence_threshold
    )


@cli.command()
@click.option('--input', required=True)
@click.option('--output', required=True)
def convert_haplotypes_csv_to_tsv(input, output):
    mondrianutils.haplotypes.convert_csv_to_tsv(
        input, output
    )


if __name__ == '__main__':
    cli()
