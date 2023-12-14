import click
import mondrianutils.helpers as helpers
import numpy as np
import pandas as pd
import yaml


def _read_vcf(vcf_file, info_field_name):
    data = {}
    info_field_name += '='
    with helpers.getFileHandle(vcf_file, 'rt') as reader:
        for line in reader:
            if line.startswith('#'):
                continue
            line = line.strip().split()
            chrom, pos, _, ref, alt = line[:5]
            svtype = [v for v in line[7].split(';') if v.startswith(info_field_name)]
            assert len(svtype) == 1
            svtype = svtype[0][len(info_field_name):]
            data[(chrom, pos, ref, alt)] = svtype
    return data


def _compare_vcf(vcf_file, ref_vcf_file, info_field_name):
    vcfdata = _read_vcf(vcf_file, info_field_name)
    ref_vcfdata = _read_vcf(ref_vcf_file, info_field_name)
    for (chrom, pos, ref, alt) in ref_vcfdata:
        assert (chrom, pos, ref, alt) in vcfdata
        assert vcfdata[(chrom, pos, ref, alt)] == ref_vcfdata[(chrom, pos, ref, alt)]


def _compare_csv(data, ref_data, approx_cols, atol=0.01, sep=','):
    data = pd.read_csv(data, sep=sep)
    ref_data = pd.read_csv(ref_data, sep=sep)

    assert sorted(data.columns.values) == sorted(ref_data.columns.values)
    cols = sorted(data.columns.values)
    data = data.sort_values(by=cols)
    ref_data = ref_data.sort_values(by=cols)
    data = data.reset_index(drop=True)
    ref_data = ref_data.reset_index(drop=True)

    for colname in data.columns.values:
        if colname in approx_cols:
            assert np.allclose(data[colname], ref_data[colname], atol=atol, equal_nan=True), colname
        else:
            assert data[colname].equals(ref_data[colname]), (data[colname], ref_data[colname])


def compare_alignment(metrics, metrics_ref, gc_metrics, gc_metrics_ref):
    metrics = pd.read_csv(metrics)
    gc_metrics = pd.read_csv(gc_metrics)
    ref_metrics = pd.read_csv(metrics_ref)
    ref_gc_metrics = pd.read_csv(gc_metrics_ref)
    for colname in metrics.columns.values:
        if colname == 'tss_enrichment_score':
            continue
        if metrics[colname].dtype == float:
            assert np.allclose(metrics[colname], ref_metrics[colname], atol=0.001)
        else:
            assert metrics[colname].equals(ref_metrics[colname])
    for colname in gc_metrics.columns.values:
        if gc_metrics[colname].dtype == float:
            assert np.allclose(gc_metrics[colname], ref_gc_metrics[colname], atol=0.001)
        else:
            assert gc_metrics[colname].equals(ref_gc_metrics[colname])


def compare_hmmcopy(reads, reads_ref, metrics, metrics_ref):
    approx_cols = [
        'MSRSI_non_integerness', 'MBRSI_dispersion_non_integerness', 'MBRSM_dispersion',
        'autocorrelation_hmmcopy', 'cv_hmmcopy', 'mad_hmmcopy',
        'total_halfiness', 'scaled_halfiness', 'mean_state_mads', 'mean_state_vars', 'mad_neutral_state',
        'mean_copy',
        'log_likelihood', 'true_multiplier', 'quality'
    ]
    _compare_csv(metrics, metrics_ref, approx_cols, atol=0.05)

    approx_cols = ['cor_gc', 'copy', 'modal_curve']
    _compare_csv(reads, reads_ref, approx_cols, atol=0.1)


def compare_variant_calling(
        museq, museq_ref, mutect, mutect_ref,
        strelka_snv, strelka_snv_ref, strelka_indel,
        strelka_indel_ref
):
    _compare_vcf(museq, museq_ref, 'PR')
    _compare_vcf(mutect, mutect_ref, 'AS_FilterStatus')
    _compare_vcf(strelka_snv, strelka_snv_ref, 'QSS')
    _compare_vcf(strelka_indel, strelka_indel_ref, 'QSI')


def compare_breakpoint_calling(
        destruct, destruct_ref, lumpy, lumpy_ref,
        gridss, gridss_ref, svaba, svaba_ref
):
    _compare_csv(destruct, destruct_ref, ['mate_score', 'log_likelihood', 'log_cdf'], sep='\t', atol=0.001)
    _compare_vcf(lumpy, lumpy_ref, 'SVTYPE')
    _compare_vcf(gridss, gridss_ref, 'SVTYPE')
    _compare_vcf(svaba, svaba_ref, 'SVTYPE')


def compare_snv_genotyping(
        genotyper, genotyper_ref, vartrix, vartrix_ref
):
    _compare_csv(genotyper, genotyper_ref, [])
    _compare_csv(vartrix, vartrix_ref, [])


def compare_sv_genotyping(
        genotyper, genotyper_ref
):
    _compare_csv(genotyper, genotyper_ref, [])


def compare_normalizer(cells_yaml):
    with open(cells_yaml, 'rt') as reader:
        data = yaml.safe_load(reader)

        assert data['cells'] == ['SA1090-A96213A-R22-C43']


