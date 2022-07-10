import os

import remixt.analysis
import remixt.analysis.haplotype
import remixt.analysis.readcount
import remixt.analysis.segment
import yaml
from mondrianutils import helpers


def infer_haps_grch37(
        thousand_genomes_dir, tempdir, output, snp_genotype, chromosome,
        genetic_map_filename_template, haplotypes_filename_template,
        legend_filename_template, sample_filename,
        chr_name_prefix, phased_chromosome_x
):
    chrom = phased_chromosome_x if chromosome.endswith('X') else chromosome

    assert '{chromosome}' in genetic_map_filename_template
    genetic_map = genetic_map_filename_template.replace('{chromosome}', chrom)
    genetic_map = os.path.join(thousand_genomes_dir, genetic_map)

    assert '{chromosome}' in haplotypes_filename_template
    haplotypes = haplotypes_filename_template.replace('{chromosome}', chrom)
    haplotypes = os.path.join(thousand_genomes_dir, haplotypes)

    assert '{chromosome}' in legend_filename_template
    legend = legend_filename_template.replace('{chromosome}', chrom)
    legend = os.path.join(thousand_genomes_dir, legend)

    sample = os.path.join(thousand_genomes_dir, sample_filename)

    config = {
        'genetic_map_template': genetic_map,
        'haplotypes_template': haplotypes,
        'legend_template': legend,
        'sample_template': sample,
        'ensembl_genome_version': 'GRCh37',
        'chr_name_prefix': chr_name_prefix
    }
    remixt.analysis.haplotype.infer_haps(
        output, snp_genotype, chromosome, tempdir,
        config, None
    )


def infer_haps_grch38(
        thousand_genomes_dir, tempdir, output, snp_genotype, chromosome,
        grch38_1kg_bcf_filename_template, grch38_1kg_X_bcf_filename_template,
        genetic_map_grch38_filename_template, chr_name_prefix, grch38_1kg_phased_chromosome_x,
        grch38_1kg_chromosomes, snp_positions
):
    assert '{chromosome}' in grch38_1kg_bcf_filename_template
    grch38_1kg_bcf_filename_template = grch38_1kg_bcf_filename_template.replace('{chromosome}', chromosome)
    grch38_1kg_bcf_filename_template = os.path.join(thousand_genomes_dir, grch38_1kg_bcf_filename_template)

    grch38_1kg_X_bcf_filename_template = os.path.join(thousand_genomes_dir, grch38_1kg_X_bcf_filename_template)

    assert '{chromosome}' in genetic_map_grch38_filename_template
    genetic_map_grch38_filename_template = genetic_map_grch38_filename_template.replace('{chromosome}', chromosome)
    genetic_map_grch38_filename_template = os.path.join(thousand_genomes_dir, genetic_map_grch38_filename_template)

    config = {
        'snp_positions_filename': snp_positions,
        'grch38_1kg_chromosomes': grch38_1kg_chromosomes,
        'ensembl_genome_version': 'GRCh38',
        'chr_name_prefix': chr_name_prefix,
        'grch38_1kg_bcf_filename_template': grch38_1kg_bcf_filename_template,
        'grch38_1kg_X_bcf_filename_template': grch38_1kg_X_bcf_filename_template,
        'grch38_1kg_phased_chromosome_x': grch38_1kg_phased_chromosome_x,
        'genetic_map_grch38_filename_template': genetic_map_grch38_filename_template
    }

    remixt.analysis.haplotype.infer_haps(
        output, snp_genotype, chromosome, tempdir,
        config, None
    )


def load_tar_meta(thousand_genomes_dir):
    tg_tar_yaml = os.path.join(thousand_genomes_dir, 'meta.yaml')
    assert os.path.exists(tg_tar_yaml)

    with open(tg_tar_yaml, 'rt') as reader:
        data = yaml.safe_load(reader)

    return data


def infer_haps(snp_genotype, chromosome, output, thousand_genomes_tar, snp_positions, tempdir):
    thousand_genomes_dir = os.path.join(tempdir, 'thousand_genomes_impute_tar')

    helpers.untar(thousand_genomes_tar, thousand_genomes_dir)

    meta = load_tar_meta(thousand_genomes_dir)

    if meta['ensembl_genome_version'] == 'GRCh37':
        infer_haps_grch37(
            thousand_genomes_dir, tempdir, output, snp_genotype, chromosome,
            meta['genetic_map_filename_template'], meta['haplotypes_filename_template'],
            meta['legend_filename_template'], meta['sample_filename'],
            meta['chr_name_prefix'], meta['phased_chromosome_x']
        )
    elif meta['ensembl_genome_version'] == 'GRCh38':
        infer_haps_grch38(
            thousand_genomes_dir, tempdir, output, snp_genotype, chromosome,
            meta['grch38_1kg_bcf_filename_template'], meta['grch38_1kg_X_bcf_filename_template'],
            meta['genetic_map_grch38_filename_template'], meta['chr_name_prefix'],
            meta['grch38_1kg_phased_chromosome_x'],
            meta['grch38_1kg_chromosomes'], snp_positions
        )
    else:
        raise Exception()
