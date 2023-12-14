def dtypes(fastqscreen_genomes=['grch37', 'mm10', 'salmon', 'human', 'mouse']):
    metrics = {
        'multiplier': 'int',
        'tss_ennrichment_score': 'float',
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
        'MSRSI_non_integerness': 'float',
        'MBRSI_dispersion_non_integerness': 'float',
        'MBRSM_dispersion': 'float',
        'autocorrelation_hmmcopy': 'float',
        'cv_hmmcopy': 'float',
        'empty_bins_hmmcopy': 'int',
        'mad_hmmcopy': 'float',
        'mean_hmmcopy_reads_per_bin': 'float',
        'median_hmmcopy_reads_per_bin': 'float',
        'std_hmmcopy_reads_per_bin': 'float',
        'total_mapped_reads_hmmcopy': 'int',
        'total_halfiness': 'float',
        'scaled_halfiness': 'float',
        'mean_state_mads': 'float',
        'mean_state_vars': 'float',
        'mad_neutral_state': 'float',
        'breakpoints': 'int',
        'mean_copy': 'float',
        'state_mode': 'int',
        'log_likelihood': 'float',
        'true_multiplier': 'float',
        'column': 'int',
        'img_col': 'int',
        'primer_i7': 'str',
        'index_i5': 'str',
        'sample_type': 'str',
        'primer_i5': 'str',
        'experimental_condition': 'str',
        'cell_call': 'str',
        'index_i7': 'str',
        'clustering_order': 'int',
        'row': 'int',
        'trim': 'bool',
        'is_s_phase_prob': 'float',
        'is_s_phase': 'bool',
        'quality': 'float',
        'total_mapped_reads': 'int',
        'unpaired_mapped_reads': 'int',
        'paired_mapped_reads': 'int',
        'unpaired_duplicate_reads': 'int',
        'paired_duplicate_reads': 'int',
        'unmapped_reads': 'int',
        'percent_duplicate_reads': 'float',
        'estimated_library_size': 'int',
        'total_reads': 'int',
        'total_duplicate_reads': 'int',
        'total_properly_paired': 'int',
        'coverage_breadth': 'float',
        'coverage_depth': 'float',
        'median_insert_size': 'float',
        'mean_insert_size': 'float',
        'standard_deviation_insert_size': 'float',
        'is_contaminated': 'bool',
        'species': 'str',
        'condition': 'str',
        'index_sequence': 'str',
        'pick_met': 'str',
        'is_control': 'bool',
        'aligned': 'float',
        'expected': 'float',
        'overlap_with_all_filters': 'float',
        'overlap_with_all_filters_and_qual': 'float',
        'overlap_with_dups': 'float',
        'overlap_without_dups': 'float',
        'fastqscreen_nohit': 'int',
        'fastqscreen_total_reads': 'int',
        'fastqscreen_nohit_ratio': 'float',
        'is_normal': 'category',
        'aneuploidy_score': 'float'
    }

    for genome in fastqscreen_genomes:
        metrics['fastqscreen_{}'.format(genome)] = 'int'
        metrics['fastqscreen_{}_multihit'.format(genome)] = 'int'
        metrics['fastqscreen_{}_ratio'.format(genome)] = 'float'

    params = {
        'iteration': 'float',
        'state': 'float',
        'parameter': 'str',
        'cell_id': 'category',
        'value': 'float',
    }

    reads = {
        'chr': 'category',
        'start': 'int',
        'end': 'int',
        'width': 'int',
        'reads': 'int',
        'gc': 'float',
        'cor_gc': 'float',
        'cor_map': 'float',
        'copy': 'float',
        'map': 'float',
        'state': 'int',
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
        'valid': 'bool',
        'ideal': 'bool',
        'modal_curve': 'float',
        'modal_quantile': 'float',
        'multiplier': 'int',
        'is_low_mappability': 'bool',
        'fraction_overlapping_reads': 'float'
    }

    segs = {
        'chr': 'category',
        'start': 'int',
        'end': 'int',
        'state': 'int',
        'median': 'float',
        'multiplier': 'int',
        'cell_id': 'category',
    }

    return locals()
