import click
import mondrianutils.variant_calling
@click.group()
def cli():
    pass


@cli.command()
@click.option('--reference', help='specify reference fasta')
@click.option('--chromosomes', multiple=True, default=list(map(str, range(1, 23))) + ['X', 'Y'],
              help='specify target chromosomes')
@click.option('--size', default=1000000, type=int, help='specify interval size')
def generate_intervals(reference, chromosomes, size):
    mondrianutils.variant_calling.generate_intervals(reference, chromosomes, size)


@cli.command()
@click.option('--interval', help='specify reference fasta')
@click.option('--num_splits', type=int, help='specify target chromosomes')
def split_interval(interval, num_splits):
    mondrianutils.variant_calling.split_interval(interval, num_splits)


@cli.command()
@click.option('--inputs', multiple=True, help='vcf files to merge')
@click.option('--output', help='merged output vcf')
def merge_vcf_files(inputs, output):
    mondrianutils.variant_calling.merge_vcf_files(inputs, output)


@cli.command()
@click.option('--reference', help='specify reference fasta')
@click.option('--chromosomes', multiple=True, default=list(map(str, range(1, 23))) + ['X', 'Y'],
              help='specify target chromosomes')
def genome_size(reference, chromosomes):
    mondrianutils.variant_calling.get_genome_size(reference, chromosomes)


@cli.command()
@click.option('--inputs', multiple=True, help='specify input')
@click.option('--output', help='specify output')
def merge_chromosome_depths_strelka(inputs, output):
    mondrianutils.variant_calling.merge_chromosome_depths_strelka(inputs, output)


@cli.command()
@click.option('--input', help='specify input bam')
def get_sample_id_bam(input):
    mondrianutils.variant_calling.get_sample_id_bam(input)


@cli.command()
@click.option('--museq_vcf', required=True)
@click.option('--mutect_vcf', required=True)
@click.option('--strelka_indel', required=True)
@click.option('--strelka_snv', required=True)
@click.option('--consensus_output', required=True)
@click.option('--counts_output', required=True)
@click.option('--chromosomes', default=[str(v) for v in range(1, 23)] + ['X', 'Y'], multiple=True)
def consensus(museq_vcf, mutect_vcf, strelka_indel, strelka_snv, consensus_output, counts_output, chromosomes):
    mondrianutils.variant_calling.consensus(museq_vcf, strelka_snv, mutect_vcf, strelka_indel, consensus_output, counts_output, chromosomes)


@cli.command()
@click.option('--input', required=True)
@click.option('--output', required=True)
@click.option('--tumour', required=True)
@click.option('--normal', required=True)
@click.option('--vcf_normal_id', required=True)
@click.option('--vcf_tumour_id', required=True)
def vcf_reheader_id(input, output, tumour, normal, vcf_normal_id, vcf_tumour_id):
    mondrianutils.variant_calling.vcf_reheader_id(input, output, tumour, normal, vcf_normal_id, vcf_tumour_id)


@cli.command()
@click.option('--input', required=True)
@click.option('--output', required=True)
@click.option('--tumour_bam')
@click.option('--normal_bam')
def update_maf_ids(input, output, tumour_bam, normal_bam):
    mondrianutils.variant_calling.update_maf_ids(input, output, tumour_bam, normal_bam)


@cli.command()
@click.option('--input', required=True)
@click.option('--counts', required=True)
@click.option('--output', required=True)
def update_maf_counts(input, counts, output):
    mondrianutils.variant_calling.update_maf_counts(input, counts, output)


@cli.command()
@click.option('--infiles', multiple=True, required=True)
@click.option('--output', required=True)
def merge_mafs(infiles, output):
    mondrianutils.variant_calling.merge_mafs(infiles, output)


@cli.command()
@click.option('--inputs', multiple=True, required=True)
@click.option('--output', required=True)
def concat_csv(inputs, output):
    mondrianutils.variant_calling.concatenate_csv(inputs, output)


@cli.command()
@click.option('--input', required=True)
@click.option('--output', required=True)
def fix_museq_vcf(input, output):
    mondrianutils.variant_calling.fix_museq_vcf(input, output)


@cli.command()
@click.option('--inputs', multiple=True, required=True)
@click.option('--output', required=True)
@click.option('--threads', type=int, required=True)
@click.option('--tempdir', required=True)
def merge_bams(inputs, output, threads, tempdir):
    mondrianutils.helpers.merge_bams(inputs, output, tempdir, threads)


@cli.command()
@click.option('--files')
@click.option('--metadata_yaml_files', multiple=True)
@click.option('--samples', multiple=True)
@click.option('--metadata_output')
def generate_metadata(files, metadata_yaml_files, samples, metadata_output):
    mondrianutils.variant_calling.generate_metadata(files, metadata_yaml_files, samples, metadata_output)


if __name__ == '__main__':
    cli()
