import numpy as np
import os
import pandas as pd
import pysam
from mondrianutils import helpers
import csverve
from mondrianutils.dtypes import haplotypes as dtypes


def read_bcf_phased_genotypes(bcf_filename):
    """ Read in a shapeit4 generated BCF file and return dataframe of phased alleles.

    Parameters
    ----------
    bcf_filename : str
        BCF file produced by shapeit4

    Returns
    -------
    pandas.DataFrame
        table of phased alleles
    """
    phased_genotypes = []

    for r in pysam.VariantFile(bcf_filename, 'r'):
        for alt in r.alts:
            chromosome = r.chrom
            position = r.pos
            ref = r.ref

            assert len(r.samples) == 1
            gt_infos = r.samples[0].items()

            assert len(gt_infos) == 1
            assert gt_infos[0][0] == 'GT'
            allele1, allele2 = gt_infos[0][1]

            phased_genotypes.append([chromosome, position, ref, alt, allele1, allele2])

    phased_genotypes = pd.DataFrame(
        phased_genotypes,
        columns=['chromosome', 'position', 'ref', 'alt', 'allele1', 'allele2'])

    return phased_genotypes


def read_phasing_samples(bcf_filenames):
    """ Read a set of phasing samples from BCF files

    Parameters
    ----------
    bcf_filenames : list of str
        list of BCF of phased SNPs

    Yields
    ------
    pandas.DataFrame
        allele1 and allele2 (0/1) indexed by chrom, coord, ref, alt
    """
    for bcf_filename in bcf_filenames:
        phasing = read_bcf_phased_genotypes(bcf_filename)
        phasing.set_index(['chromosome', 'position', 'ref', 'alt'], inplace=True)
        yield phasing


def calculate_haplotypes(phasing_samples, changepoint_threshold=0.95):
    """ Calculate haplotype from a set phasing samples.

    Parameters
    ----------
    phasing_samples : list of pandas.Series
        set of phasing samples for a set of SNPs
    changepoint_threshold : float, optional
        threshold on high confidence changepoint calls, by default 0.95

    Returns
    ------
    pandas.DataFrame
        haplotype info with columns:
            chromosome, position, ref, alt, fraction_changepoint, changepoint_confidence,
            is_changepoint, not_confident, chrom_different, hap_label, allele1, allele2
    """

    haplotypes = None
    n_samples = 0

    for phasing in phasing_samples:
        # Select het positions
        phasing = phasing[phasing['allele1'] != phasing['allele2']]

        # Identify changepoints.  A changepoint occurs when the alternate allele
        # of a heterozygous SNP is on a different haplotype allele from the alternate
        # allele of the previous het SNP.
        changepoints = phasing['allele1'].diff().abs().astype(float).fillna(0.0)

        if haplotypes is None:
            haplotypes = changepoints
        else:
            haplotypes += changepoints
        n_samples += 1

    haplotypes /= float(n_samples)

    haplotypes = haplotypes.rename('fraction_changepoint').reset_index()

    # Calculate confidence in either changepoint or no changepoint
    haplotypes['changepoint_confidence'] = np.maximum(haplotypes['fraction_changepoint'],
                                                      1.0 - haplotypes['fraction_changepoint'])

    # Calculate most likely call of changepoint or no changepoint
    haplotypes['is_changepoint'] = haplotypes['fraction_changepoint'].round().astype(int)

    # Threshold confident changepoint calls
    haplotypes['not_confident'] = (haplotypes['changepoint_confidence'] < float(changepoint_threshold))

    # Calculate hap label
    haplotypes['chrom_different'] = haplotypes['chromosome'].ne(haplotypes['chromosome'].shift())
    haplotypes['hap_label'] = (haplotypes['not_confident'] | haplotypes['chrom_different']).cumsum() - 1

    # Calculate most likely alelle1
    haplotypes['allele1'] = haplotypes['is_changepoint'].cumsum().mod(2)
    haplotypes['allele2'] = 1 - haplotypes['allele1']

    return haplotypes


def write_null(haps_filename):
    with open(haps_filename, 'w') as haps_file:
        haps_file.write('chromosome\tposition\tallele\thap_label\tallele_id\n')


def run_shapeit4(
        input_calls_bcf, genetic_map, regions,
        chromosome, temp_directory, haps_filename,
        phased_chromosomes=tuple([f'chr{v}' for v in range(1, 23)] + ['chrX']),
        is_female=True,
        phased_chromx='chrX',
        shapeit_num_samples=100,
        shapeit_confidence_threshold=0.95
):
    """ Infer haplotype blocks for a chromosome using shapeit4 for grch38

    Parameters
    ----------
    input_calls_bcf : bcf file with calls from bcftools call
    genetic_map : reference genetic map file
    regions : phased reference vcf per chromosome
    chromosome : chromosome name
    temp_directory : tempdir
    haps_filename : output haps file
    phased_chromosomes : list of chromosomes that are phased
    is_female : True if sample is from female
    phased_chromx : name of chrx in reference
    shapeit_num_samples : num of samples/iterations
    shapeit_confidence_threshold : confidence threshold to select sample

    Returns
    -------
    The output haps file will contain haplotype blocks for each heterozygous SNP position. The
    file will be TSV format with the following columns:

        'chromosome': het snp chromosome
        'position': het snp position
        'allele': binary indicator for reference (0) vs alternate (1) allele
        'hap_label': label of the haplotype block
        'allele_id': binary indicator of the haplotype allele

    """

    chr_name_prefix = 'chr' if chromosome.startswith('chr') else ''

    chromosome_1kg = f'chr{chromosome}' if chr_name_prefix == '' else chromosome

    # Skip unphased chromosomes
    if str(chromosome_1kg) not in phased_chromosomes:
        write_null(haps_filename)
        return

    # If we are analyzing male data and this is chromosome X
    # then there are no het snps and no haplotypes
    if chromosome == phased_chromx and not is_female:
        write_null(haps_filename)
        return

    # Temporary directory for shapeit files
    helpers.makedirs(temp_directory)

    # Run shapeit to generate phasing graph
    bingraph_filename = os.path.join(temp_directory, 'phasing.bingraph')

    helpers.run_cmd(
        ['shapeit4',
         '--input', input_calls_bcf,
         '--map', genetic_map,
         '--region', chromosome,
         '--reference', regions,
         '--bingraph', bingraph_filename]
    )

    # Run shapeit to sample from phased haplotype graph
    sample_template = os.path.join(temp_directory, 'sampled.{0}.bcf')
    sample_filenames = []
    for s in range(shapeit_num_samples):
        sample_filename = sample_template.format(s)
        sample_filenames.append(sample_filename)
        helpers.run_cmd(
            ['bingraphsample',
             '--input', bingraph_filename,
             '--output', sample_filename,
             '--sample',
             '--seed', str(s)]
        )
        helpers.run_cmd(
            ['bcftools', 'index', '-f', sample_filename]
        )

    haplotypes = calculate_haplotypes(
        read_phasing_samples(sample_filenames),
        changepoint_threshold=shapeit_confidence_threshold
    )

    haplotypes = pd.concat([
        haplotypes.rename(columns={'allele1': 'allele'})[['chromosome', 'position', 'allele', 'hap_label']].assign(
            allele_id=0),
        haplotypes.rename(columns={'allele2': 'allele'})[['chromosome', 'position', 'allele', 'hap_label']].assign(
            allele_id=1),
    ])

    haplotypes = haplotypes[['chromosome', 'position', 'allele', 'hap_label', 'allele_id']]

    csverve.write_dataframe_to_csv_and_yaml(haplotypes, haps_filename, dtypes=dtypes.dtypes())
