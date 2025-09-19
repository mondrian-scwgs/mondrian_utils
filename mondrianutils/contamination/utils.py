import os
import pandas as pd
import numpy as np
import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import re
from collections import defaultdict
from io import StringIO


def parse_bamstat_old(fname):
    """Parse BAM statistics from samtools stats output file."""
    skip_fields = set(['#', 'CHK', 'GCT', 'GCC', 'GCT', 'FBC', 'FTC', 'LBC', 'LTC', 'BCC', 'CRC', 'OXC', 'RXC', 'GCD'])
    lines = open(fname).readlines()
    result = {}

    summary = {}
    i = _skip_ahead(lines, 0, skip_fields)

    # summary fields
    while lines[i].startswith('SN'):
        tkns = lines[i].split('\t')
        summary[tkns[1]] = float(tkns[2].strip())
        i += 1
    result['summary'] = summary

    i = _skip_ahead(lines, i, skip_fields)
    if summary['reads mapped and paired:'] == 0:
        return result

    ffq = []
    while lines[i].startswith('FFQ'):
        tkns = lines[i].strip().split('\t')
        ffq.append(np.array(tkns[2:]).astype(int))
        i += 1
    ffq = np.array(ffq)
    result['first_fragment_quality'] = ffq

    i = _skip_ahead(lines, i, skip_fields)

    lfq = []
    while lines[i].startswith('LFQ'):
        tkns = lines[i].strip().split('\t')
        lfq.append(np.array(tkns[2:]).astype(int))
        i += 1
    lfq = np.array(lfq)
    result['last_fragment_quality'] = lfq

    i = _skip_ahead(lines, i, skip_fields)

    gcf = []
    while lines[i].startswith('GCF'):
        tkns = lines[i].strip().split('\t')
        gcf.append(np.array(tkns[1:]))
        i += 1
    gcf = np.array(gcf)
    result['first_fragment_gc'] = gcf

    i = _skip_ahead(lines, i, skip_fields)

    gcl = []
    while lines[i].startswith('GCL'):
        tkns = lines[i].strip().split('\t')
        gcl.append(np.array(tkns[1:]))
        i += 1
    gcl = np.array(gcl)
    result['last_fragment_gc'] = gcl

    # skip the next few fields
    i = _skip_ahead(lines, i, skip_fields)

    insert_sizes = []
    while lines[i].startswith('IS'):
        tkns = lines[i].strip().split('\t')
        insert_sizes.append({'insert_size': int(tkns[1]),
                             'pairs_total': int(tkns[2]),
                             'inward_oriented_pairs': int(tkns[3]),
                             'outward_oriented_pairs': int(tkns[4]),
                             'other_pairs': int(tkns[5])
                             })
        i += 1
    insert_sizes = pd.DataFrame(insert_sizes)
    result['insert_sizes'] = insert_sizes

    i = _skip_ahead(lines, i, skip_fields)

    rl = []
    while lines[i].startswith('RL'):
        tkns = lines[i].strip().split('\t')
        rl.append({'read_length': int(tkns[1]),
                   'count': int(tkns[2]),
                   })
        i += 1
    rl = pd.DataFrame(rl)
    result['read_lengths'] = rl

    i = _skip_ahead(lines, i, skip_fields)

    frl = []
    while lines[i].startswith('FRL'):
        tkns = lines[i].strip().split('\t')
        frl.append({'fragment_length': int(tkns[1]),
                    'count': int(tkns[2]),
                    })
        i += 1
    frl = pd.DataFrame(frl)
    result['first_fragment_lengths'] = frl

    i = _skip_ahead(lines, i, skip_fields)

    lrl = []
    while lines[i].startswith('LRL'):
        tkns = lines[i].strip().split('\t')
        lrl.append({'fragment_length': int(tkns[1]),
                    'count': int(tkns[2]),
                    })
        i += 1
    lrl = pd.DataFrame(lrl)
    result['last_fragment_lengths'] = lrl

    i = _skip_ahead(lines, i, skip_fields)

    mapq = []
    while lines[i].startswith('MAPQ'):
        tkns = lines[i].strip().split('\t')
        mapq.append({'mapq': int(tkns[1]),
                     'count': int(tkns[2]),
                     })
        i += 1
    mapq = pd.DataFrame(mapq)
    result['mapq'] = mapq

    i = _skip_ahead(lines, i, skip_fields)

    indels = []
    while lines[i].startswith('ID'):
        tkns = lines[i].strip().split('\t')
        indels.append({'length': int(tkns[1]),
                       'n_insertions': int(tkns[2]),
                       'n_deletions': int(tkns[3])
                       })
        i += 1
    indels = pd.DataFrame(indels)
    result['indels'] = indels

    i = _skip_ahead(lines, i, skip_fields)

    indels_per_cycle = []
    while lines[i].startswith('IC'):
        tkns = lines[i].strip().split('\t')
        indels_per_cycle.append({'cycle': int(tkns[1]),
                                 'n_insertions_fwd': int(tkns[2]),
                                 'n_insertions_rev': int(tkns[3]),
                                 'n_deletions_fwd': int(tkns[4]),
                                 'n_deletions_rev': int(tkns[5]),
                                 })
        i += 1
    indels_per_cycle = pd.DataFrame(indels_per_cycle)
    result['indels_per_cycle'] = indels_per_cycle

    i = _skip_ahead(lines, i, skip_fields)

    cov = []
    while lines[i].startswith('COV'):
        tkns = lines[i].strip().split('\t')
        cov.append({'cov': int(tkns[2]),
                    'n_bases': int(tkns[3]),
                    })
        i += 1
    cov = pd.DataFrame(cov)
    result['coverage'] = cov

    return result


def _skip_ahead(lines, i, skip_fields):
    """Helper function to skip lines in BAM stats file."""
    while i < len(lines) and (lines[i].startswith('#') or lines[i].split('\t')[0] in skip_fields):
        i += 1
    return i


def parse_bamstat(fname):
    """Parse BAM statistics from samtools stats output file using a refactored approach."""
    # Read all lines and organize by section
    with open(fname, 'r') as f:
        lines = f.readlines()
    
    # Parse lines into sections
    sections = defaultdict(list)
    skip_fields = set(['CHK', 'GCT', 'GCC', 'GCT', 'FBC', 'FTC', 'LBC', 'LTC', 'BCC', 'CRC', 'OXC', 'RXC', 'GCD'])
    
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        parts = line.split('\t', maxsplit=1)
        if len(parts) < 2:
            continue

        section_key = parts[0]
        data = parts[1]

        if section_key in skip_fields:
            continue

        sections[section_key].append(data)
    
    result = {}
    
    # Process SN (Summary Numbers) section
    if 'SN' in sections:
        summary = {}
        for line in sections['SN']:
            parts = line.split('\t')
            key = parts[0]
            value = float(parts[1].split('#')[0].strip())  # Remove comments
            summary[key] = value
        result['summary'] = summary
    
    # Early return if no mapped paired reads
    if result.get('summary', {}).get('reads mapped and paired:', 0) == 0:
        return result
    
    # Process FFQ (First Fragment Qualities) - matrix data
    if 'FFQ' in sections:
        ffq_data = []
        for line in sections['FFQ']:
            values = line.split('\t')
            ffq_data.append([int(x) for x in values[1:]])  # Skip cycle number
        result['first_fragment_quality'] = np.array(ffq_data)
    
    # Process LFQ (Last Fragment Qualities) - matrix data
    if 'LFQ' in sections:
        lfq_data = []
        for line in sections['LFQ']:
            values = line.split('\t')
            lfq_data.append([int(x) for x in values[1:]])  # Skip cycle number
        result['last_fragment_quality'] = np.array(lfq_data)
    
    # Process GCF (First Fragment GC Content) - matrix data
    if 'GCF' in sections:
        gcf_data = []
        for line in sections['GCF']:
            values = line.split('\t')
            gcf_data.append(values)  # Keep all values including cycle
        result['first_fragment_gc'] = np.array(gcf_data)
    
    # Process GCL (Last Fragment GC Content) - matrix data
    if 'GCL' in sections:
        gcl_data = []
        for line in sections['GCL']:
            values = line.split('\t')
            gcl_data.append(values)  # Keep all values including cycle
        result['last_fragment_gc'] = np.array(gcl_data)
    
    # Process IS (Insert Sizes) - tabular data
    if 'IS' in sections:
        csv_data = StringIO('\n'.join(sections['IS']))
        insert_sizes = pd.read_csv(csv_data, sep='\t', names=['insert_size', 'pairs_total', 'inward_oriented_pairs', 'outward_oriented_pairs', 'other_pairs'])
        result['insert_sizes'] = insert_sizes
    
    # Process RL (Read Lengths) - tabular data
    if 'RL' in sections:
        csv_data = StringIO('\n'.join(sections['RL']))
        read_lengths = pd.read_csv(csv_data, sep='\t', names=['read_length', 'count'])
        result['read_lengths'] = read_lengths
    
    # Process FRL (First Fragment Lengths) - tabular data
    if 'FRL' in sections:
        csv_data = StringIO('\n'.join(sections['FRL']))
        first_fragment_lengths = pd.read_csv(csv_data, sep='\t', names=['fragment_length', 'count'])
        result['first_fragment_lengths'] = first_fragment_lengths
    
    # Process LRL (Last Fragment Lengths) - tabular data
    if 'LRL' in sections:
        csv_data = StringIO('\n'.join(sections['LRL']))
        last_fragment_lengths = pd.read_csv(csv_data, sep='\t', names=['fragment_length', 'count'])
        result['last_fragment_lengths'] = last_fragment_lengths
    
    # Process MAPQ (Mapping Qualities) - tabular data
    if 'MAPQ' in sections:
        csv_data = StringIO('\n'.join(sections['MAPQ']))
        mapq = pd.read_csv(csv_data, sep='\t', names=['mapq', 'count'])
        result['mapq'] = mapq
    
    # Process ID (Indels) - tabular data
    if 'ID' in sections:
        csv_data = StringIO('\n'.join(sections['ID']))
        indels = pd.read_csv(csv_data, sep='\t', names=['length', 'n_insertions', 'n_deletions'])
        result['indels'] = indels
    
    # Process IC (Indels per Cycle) - tabular data
    if 'IC' in sections:
        csv_data = StringIO('\n'.join(sections['IC']))
        indels_per_cycle = pd.read_csv(csv_data, sep='\t', names=['cycle', 'n_insertions_fwd', 'n_insertions_rev', 'n_deletions_fwd', 'n_deletions_rev'])
        result['indels_per_cycle'] = indels_per_cycle
    
    # Process COV (Coverage) - tabular data
    if 'COV' in sections:
        csv_data = StringIO('\n'.join(sections['COV']))
        coverage = pd.read_csv(csv_data, sep='\t', names=['rg', 'cov', 'n_bases'])
        # Keep only the cov and n_bases columns
        coverage = coverage[['cov', 'n_bases']]
        result['coverage'] = coverage
    
    return result


def parse_kraken_output(kraken_output_file, output_table, output_human, output_nonhuman):
    """Parse Kraken2 output and separate human vs non-human reads."""
    df = pd.read_table(kraken_output_file,
                       names=['is_classified', 'qname', 'result', 'lengths', 'kmers'])

    my_re = re.compile(r'([A-Za-z0-9 -_\.]+) \(taxid ([0-9]+)\)')

    taxon = []
    taxon_id = []
    for t in df.result:
        match = my_re.search(t)
        taxon.append(match[1])
        taxon_id.append(int(match[2]))

    df['taxon'] = taxon
    df['taxon_id'] = taxon_id

    human_reads = df[df.taxon_id == 9606]
    nonhuman_reads = df[df.taxon_id != 9606]

    df.to_csv(output_table, index=False)
    human_reads.qname.to_csv(output_human, index=False, header=False)
    nonhuman_reads.qname.to_csv(output_nonhuman, index=False, header=False)


def generate_contamination_table_figures(
        kraken_report_files,
        all_reads_stats_files,
        human_reads_stats_files,
        nonhuman_reads_stats_files,
        hmmcopy_metrics_filename,
        library_id,
        summary_table_output,
        multipanel_figure_output,
        chip_figure_output,
        control_cells_output,
        nonhuman_percentage_taxon_output,
        nonhuman_percentage_clade_output,
        nonhuman_composition_output,
        contam_by_column_output,
        ncbi_taxonomy_database='/data1/shahs3/users/myersm2/repos/contamination/kraken_db/ncbi_taxonomy/taxa.sqlite',
        min_percent_aggregate=0.0,
        min_percent_show=2.0,
        min_num_taxa_condense=25
    ):
    """Generate contamination analysis tables and figures from pipeline outputs."""
    # this package requires a taxonomy database file -- use mine: /data1/shahs3/users/myersm2/repos/contamination/kraken_db/ncbi_taxonomy/taxa.sqlite
    from ete3 import NCBITaxa
    
    # summarize kraken2 reports and samtools stats results
    print('started')

    df = _get_summary_table_from_files(
        kraken_report_files,
        all_reads_stats_files,
        human_reads_stats_files,
        nonhuman_reads_stats_files,
        hmmcopy_metrics_filename
    )
    df['prop_unmapped_nonhuman'] = df.nonhuman_reads_unmapped / df.unmapped_reads
    df['prop_mapped_nonhuman'] = (df.nonhuman_reads_bamstats - df.nonhuman_reads_unmapped) / df.total_mapped_reads

    print('got summary table')

    df.to_csv(summary_table_output, index=False)

    print('wrote summary table')


    # make 6-pane plot
    _plot_summary_multipanel(df, library_id, fname=multipanel_figure_output)

    print('generated multipanel summary figure')


    # plot on chip
    fig, ax = plt.subplots(figsize=(5,4), dpi = 200)
    sns.scatterplot(data=df, x='col', y='row', hue='kraken2_prop_nonhuman', s = 15, hue_norm=(0,1), ax=ax)
    plt.tight_layout()
    plt.savefig(chip_figure_output)
    
    print('generated chip figure')

    ## load detailed kraken2 results
    cell_list = sorted(df.cell_id.unique())
    pct = _aggregate_field_from_files(kraken_report_files, cell_list, 'percentage_taxon', id_vars=['scientific_name'], min_val=min_percent_aggregate)
    print(pct.shape)
    print('loaded pct')

    pcc = _aggregate_field_from_files(kraken_report_files, cell_list, 'percentage_clade', id_vars=['scientific_name'], min_val=min_percent_aggregate)
    print(pcc.shape)
    print('loaded pcc')

    metrics = pd.read_csv(hmmcopy_metrics_filename)
    control_cells = metrics[metrics.is_control]
    negative_control_cells = control_cells[control_cells.condition.str.lower().str.contains('negative')].cell_id
    positive_control_cells = control_cells[~control_cells.condition.str.lower().str.contains('negative')].cell_id

    print(metrics.condition.value_counts())

    fig, axes = plt.subplots(1, 2, figsize=(8,4), dpi = 150, sharex='all', sharey = 'all')
    # plot negative control cells
    axes[0].scatter(pct[negative_control_cells].loc['Homo sapiens'],
                    pct[negative_control_cells].loc['unclassified'], s = 10)
    axes[0].set_xlabel("Percent human")
    axes[0].set_ylabel("Percent unclassified")
    axes[0].set_title(f'Negative controls (n={len(negative_control_cells)})')

    # plot positive control cells
    axes[1].scatter(pct[positive_control_cells].loc['Homo sapiens'],
                    pct[positive_control_cells].loc['unclassified'], s = 10)
    axes[1].set_xlabel("Percent human")
    axes[1].set_ylabel("Percent unclassified")
    axes[1].set_title(f'Positive controls (n={len(positive_control_cells)})')
    plt.suptitle(library_id)
    plt.tight_layout()
    plt.savefig(control_cells_output, dpi=200)
    print('plotted control cells')


    ## plot decomposition of non-human reads
    human_taxa = ['Eukaryota', 'Opisthokonta', 'Metazoa', 'Eumetazoa', 'Bilateria', 'Deuterostomia', 'Chordata', 'Craniata', 'Vertebrata', 'Gnathostomata', 'Teleostomi', 'Euteleostomi', 'Sarcopterygii', 'Dipnotetrapodomorpha', 'Tetrapoda', 'Amniota', 'Mammalia', 'Theria', 'Eutheria', 'Boreoeutheria', 'Euarchontoglires', 'Primates', 'Haplorrhini', 'Simiiformes', 'Catarrhini', 'Hominoidea', 'Hominidae', 'Homininae', 'Homo', 'Homo sapiens']
    exclude_taxa = ['unclassified', 'root', 'cellular organisms']

    # restrict to non-human reads
    bac_pcc = pcc[~pcc.index.isin(exclude_taxa + human_taxa)]
    bac_pct = pct[~pct.index.isin(exclude_taxa + human_taxa)]
    
    # make sure bac_pcc and bac_pct are normalized to 100%
    bac_pcc = bac_pcc * (100/bac_pcc.loc['Bacteria'])
    bac_pct = 100 * bac_pct / bac_pct.sum(axis = 0)


    if pct.shape[0] >= min_num_taxa_condense:
        # load NCBI taxonomy database
        ncbi = NCBITaxa(dbfile=ncbi_taxonomy_database)

        # get mapping from taxon ID (int used in database) to scientific name (str) from Kraken2 output
        cr = _aggregate_field_from_files(kraken_report_files, cell_list, 'number_clade', id_vars=['scientific_name', 'taxon_id', 'taxon_rank'])
        sciname2taxid = {r.scientific_name:r.taxon_id for _, r in cr.iloc[:, :0].reset_index().iterrows()}

        condensed = _condense_table(bac_pcc, bac_pct, min_percent_show, ncbi, sciname2taxid)
    else:
        condensed = bac_pct.copy()
    bac_pct.to_csv(nonhuman_percentage_taxon_output)
    bac_pcc.to_csv(nonhuman_percentage_clade_output)

    _plot_triple_cell_composition(condensed, pcc, df.set_index('cell_id').kraken2_total_classified,
                                 fname=nonhuman_composition_output)

    print('plotted cell composition figure')

    plt.figure(dpi = 300, figsize=(6,8))
    plt.subplot(2, 1, 1)
    sns.boxplot(data=df, x='col', y ='kraken2_prop_nonhuman')
    plt.xticks(*plt.xticks(), rotation=90)
    plt.title(library_id)

    plt.subplot(2, 1, 2)
    sns.scatterplot(df, x = 'total_reads_hmmcopy', y='kraken2_prop_nonhuman',
                hue='is_control', s = 10)
    plt.xticks(*plt.xticks(), rotation=90)
    plt.tight_layout()
    plt.savefig(contam_by_column_output, dpi=200)


# Helper functions for contamination analysis (originally from generate_table_figures.py)
def _read_kraken_report(filename):
    """Read and parse Kraken2 report file."""
    report = pd.read_table(filename, names=['percentage_clade', 'number_clade', 'number_taxon', 'taxon_rank', 'taxon_id', 'scientific_name'])
    report.scientific_name = report.scientific_name.str.strip()
    r = report[report.scientific_name == 'root'].iloc[0]
    ratio = r.percentage_clade / r.number_clade
    report['percentage_taxon'] = np.round(report.number_taxon * ratio, 2)
    return report


def _extract_cell_id_from_filename(filepath, suffix):
    """Extract cell ID from filename by removing the suffix."""
    filename = os.path.basename(filepath)
    if not filename.endswith(suffix):
        raise ValueError(f"File {filename} does not end with expected suffix {suffix}")
    return filename[:-len(suffix)]


def _create_file_mappings(file_lists_and_suffixes):
    """Create cell_id -> filepath mappings from lists of files with expected suffixes."""
    mappings = {}
    
    for file_list, suffix, mapping_name in file_lists_and_suffixes:
        cell_to_file = {}
        for filepath in file_list:
            try:
                cell_id = _extract_cell_id_from_filename(filepath, suffix)
                cell_to_file[cell_id] = filepath
            except ValueError as e:
                print(f"Warning: Skipping file due to naming issue: {e}")
                continue
        mappings[mapping_name] = cell_to_file
    
    return mappings


def _get_complete_cells(file_mappings):
    """Get list of cell IDs that have all required file types."""
    if not file_mappings:
        return []
    
    # Find intersection of all cell ID sets
    cell_sets = [set(mapping.keys()) for mapping in file_mappings.values()]
    complete_cells = set.intersection(*cell_sets)
    
    return sorted(list(complete_cells))


def _get_summary_table_from_files(kraken_report_files, all_reads_stats_files, human_reads_stats_files, nonhuman_reads_stats_files, hmmcopy_metrics_filename):
    """Build a summary table for the given library from explicit file lists."""
    
    # Define file types and their expected suffixes
    file_specs = [
        (kraken_report_files, '_report.txt', 'reports'),
        (all_reads_stats_files, '_all_reads_stats.txt', 'all_stats'),
        (human_reads_stats_files, '_human_reads_stats.txt', 'human_stats'),
        (nonhuman_reads_stats_files, '_nonhuman_reads_stats.txt', 'nonhuman_stats')
    ]
    
    # Create file mappings
    file_mappings = _create_file_mappings(file_specs)
    
    # Get cells that have all required files
    complete_cells = _get_complete_cells(file_mappings)
    
    if not complete_cells:
        print("Warning: No cells found with complete file sets")
        return pd.DataFrame()
    
    # Process each cell
    summary_rows = []
    for cell_id in tqdm.tqdm(complete_cells):
        try:
            # Process Kraken data
            kraken_data = _process_kraken_data(file_mappings['reports'][cell_id])
            
            # Process BAM stats data
            bam_data = _process_bam_stats_data(
                cell_id,
                file_mappings['all_stats'][cell_id],
                file_mappings['human_stats'][cell_id],
                file_mappings['nonhuman_stats'][cell_id]
            )
            
            # Combine data for this cell
            summary_rows.append({**kraken_data, **bam_data})
            
        except Exception as e:
            print(f"Error processing cell {cell_id}: {e}")
            continue
    
    if not summary_rows:
        print("Warning: No cells were successfully processed")
        return pd.DataFrame()
    
    # Create main DataFrame and add computed columns
    df = pd.DataFrame(summary_rows)
    df = _add_computed_columns(df)
    
    # Merge with HMMcopy metrics
    df = _merge_with_hmmcopy_metrics(df, hmmcopy_metrics_filename)
    
    # Add spatial coordinates (row/col from cell_id)
    df = _add_spatial_coordinates(df)
    
    return df


def _process_kraken_data(kraken_report_file):
    """Process Kraken report file and extract key metrics."""
    kraken_report = _read_kraken_report(kraken_report_file).set_index('scientific_name')
    kraken_total_reads = kraken_report.number_taxon.sum()
    
    if 'unclassified' in kraken_report.index:
        kraken_total_classified = kraken_total_reads - kraken_report.loc['unclassified', 'number_taxon']
    else:
        kraken_total_classified = kraken_total_reads
    
    kraken_total_human = kraken_report.loc['Homo sapiens', 'number_taxon']
    kraken_prop_human = kraken_total_human / kraken_total_classified if kraken_total_classified > 0 else 0
    kraken_prop_nonhuman = (kraken_total_classified - kraken_total_human) / kraken_total_classified if kraken_total_classified > 0 else 0
    
    return {
        'kraken2_total_fragments': kraken_total_reads,
        'kraken2_total_human': kraken_total_human,
        'kraken2_total_classified': kraken_total_classified,
        'kraken2_prop_human': kraken_prop_human,
        'kraken2_prop_nonhuman': kraken_prop_nonhuman,
    }


def _process_bam_stats_data(cell_id, all_stats_file, human_stats_file, nonhuman_stats_file):
    """Process BAM stats files and extract key metrics."""
    # Parse BAM stats files
    all_bam_stats = parse_bamstat(all_stats_file)
    human_bam_stats = parse_bamstat(human_stats_file)
    nonhuman_bam_stats = parse_bamstat(nonhuman_stats_file)
    
    # Define fields to extract from summary
    useful_fields = {
        'reads unmapped:': int, 'reads MQ0:': int, 'reads properly paired:': int, 'reads duplicated:': int,
        'insert size average:': float, 'insert size standard deviation:': float, 'mismatches:': float, 'error rate:': float,
        'average length:': float, 'percentage of properly paired reads (%):': float
    }
    
    # Extract basic metrics
    result = {'cell_id': cell_id}
    for field, dtype in useful_fields.items():
        key = field.strip(':').replace(' ', '_')
        result[f'human_{key}'] = dtype(human_bam_stats['summary'][field])
        result[f'nonhuman_{key}'] = dtype(nonhuman_bam_stats['summary'][field])
    
    # Add dummy data for missing sections
    _add_dummy_data_if_missing(human_bam_stats, nonhuman_bam_stats)
    
    # Calculate derived metrics
    result.update(_calculate_derived_metrics(all_bam_stats, human_bam_stats, nonhuman_bam_stats))
    
    return result


def _add_dummy_data_if_missing(human_bam_stats, nonhuman_bam_stats):
    """Add dummy data for missing indels/mapq sections."""
    dummy_indel_table = pd.DataFrame({'length': [1], 'n_insertions': [0], 'n_deletions': [0]})
    dummy_mapq_table = pd.DataFrame({'mapq': [0], 'count': [1]})
    
    for stats in [human_bam_stats, nonhuman_bam_stats]:
        if 'indels' not in stats or len(stats['indels']) == 0:
            stats['indels'] = dummy_indel_table.copy()
        if 'mapq' not in stats or len(stats['mapq']) == 0:
            stats['mapq'] = dummy_mapq_table.copy()


def _calculate_derived_metrics(all_bam_stats, human_bam_stats, nonhuman_bam_stats):
    """Calculate derived metrics from BAM stats."""
    result = {}
    
    # MAPQ averages
    result['human_average_mapq'] = np.prod(human_bam_stats['mapq'], axis=1).sum() / human_bam_stats['mapq']['count'].sum()
    result['nonhuman_average_mapq'] = np.prod(nonhuman_bam_stats['mapq'], axis=1).sum() / nonhuman_bam_stats['mapq']['count'].sum()
    
    # Indel metrics
    result['human_n_indels'] = human_bam_stats['indels'].n_insertions.sum() + human_bam_stats['indels'].n_deletions.sum()
    result['nonhuman_n_indels'] = nonhuman_bam_stats['indels'].n_insertions.sum() + nonhuman_bam_stats['indels'].n_deletions.sum()
    
    # Average indel lengths
    human_indel_length = (
        np.prod(human_bam_stats['indels'][['length', 'n_insertions']], axis=1).sum() +
        np.prod(human_bam_stats['indels'][['length', 'n_deletions']], axis=1).sum()
    ) / np.maximum(1, result['human_n_indels'])
    
    nonhuman_indel_length = (
        np.prod(nonhuman_bam_stats['indels'][['length', 'n_insertions']], axis=1).sum() +
        np.prod(nonhuman_bam_stats['indels'][['length', 'n_deletions']], axis=1).sum()
    ) / np.maximum(1, result['nonhuman_n_indels'])
    
    result['human_average_indel_length'] = human_indel_length
    result['nonhuman_average_indel_length'] = nonhuman_indel_length
    
    # Read counts
    result['all_reads_bamstats'] = int(all_bam_stats['summary']['raw total sequences:'])
    result['human_reads_bamstats'] = int(human_bam_stats['summary']['raw total sequences:'])
    result['nonhuman_reads_bamstats'] = int(nonhuman_bam_stats['summary']['raw total sequences:'])
    
    # Mapping metrics
    result['total_mapped_reads'] = all_bam_stats['summary']['reads unmapped:']
    result['prop_reads_unmapped'] = result['total_mapped_reads'] / np.maximum(result['all_reads_bamstats'], 1)
    
    return result


def _add_computed_columns(df):
    """Add computed columns to the main DataFrame."""
    df['total_reads_bamstats'] = df.human_reads_bamstats + df.nonhuman_reads_bamstats
    df['prop_nonhuman_bamstats'] = df.nonhuman_reads_bamstats / df.total_reads_bamstats
    return df


def _merge_with_hmmcopy_metrics(df, hmmcopy_metrics_filename):
    """Merge DataFrame with HMMcopy metrics."""
    hmmcopy_metrics = pd.read_csv(hmmcopy_metrics_filename)
    hmmcopy_metrics['prop_reads_mapped'] = hmmcopy_metrics.total_mapped_reads / hmmcopy_metrics.total_reads
    hmmcopy_metrics['prop_reads_unmapped'] = 1 - hmmcopy_metrics.prop_reads_mapped
    hmmcopy_metrics['total_reads_hmmcopy'] = hmmcopy_metrics.total_reads
    
    return df.merge(hmmcopy_metrics, on=['cell_id'], how='left', suffixes=('_bamstats', ''))


def _add_spatial_coordinates(df):
    """Add row/col coordinates extracted from cell_id."""
    df['row'] = df.cell_id.str.split('-', expand=True)[1].str.slice(1).astype(int)
    df['col'] = df.cell_id.str.split('-', expand=True)[2].str.slice(1).astype(int)
    return df


def _plot_summary_multipanel(df, title, figsize=(8,6), dpi=300, fname=None):
    """Plot some useful summary stats."""
    plt.figure(dpi = 300, figsize = (12, 8))
    plt.subplot(2, 3, 1)
    sns.scatterplot(data=df, x ='total_reads', y ='kraken2_prop_nonhuman', hue = 'is_control')
    plt.ylim(-0.02, 1.02)

    plt.subplot(2, 3, 2)
    sns.scatterplot(data=df, x ='prop_reads_unmapped', y ='kraken2_prop_nonhuman', hue = 'is_control')
    plt.xlim(-0.02, 1.02)
    plt.ylim(-0.02, 1.02)
    
    plt.subplot(2, 3, 3)
    sns.scatterplot(data=df, x ='quality', y ='kraken2_prop_nonhuman', hue = 'is_control')
    plt.xlim(-0.02, 1.02)
    plt.ylim(-0.02, 1.02)
    
    plt.subplot(2, 3, 4)
    sns.scatterplot(data=df, x='prop_reads_unmapped', y='prop_mapped_nonhuman', s = 15, hue='quality')
    plt.ylabel("proportion of mapped reads\nthat map to nonhuman taxa")

    plt.subplot(2, 3, 5)
    sns.scatterplot(data=df, hue='quality', y='prop_unmapped_nonhuman', s = 15, x='human_insert_size_average')
    plt.ylabel("proportion of UNmapped reads\nthat map to nonhuman taxa")
    
    plt.subplot(2, 3, 6)
    sns.scatterplot(data=df, x ='human_insert_size_average', y ='nonhuman_insert_size_average', hue = 'quality')
    
    plt.suptitle(title)
    plt.tight_layout()
    if fname is not None:
        plt.savefig(fname, dpi=dpi)


def _aggregate_field_from_files(kraken_report_files, cell_list, field, id_vars=['scientific_name'], min_val=0):
    """Aggregate field across cells using explicit file paths."""
    # Create mapping from cell_id to report file
    cell_to_file = {}
    for file_path in kraken_report_files:
        # Extract cell_id from filename (assuming format: {cell_id}_report.txt)
        filename = os.path.basename(file_path)
        if filename.endswith('_report.txt'):
            cell_id = filename.replace('_report.txt', '')
            cell_to_file[cell_id] = file_path

    cell_list = sorted(cell_list)

    combined_report = {}
    for cell_id in cell_list:
        report = _read_kraken_report(cell_to_file[cell_id])
        report = report.groupby(id_vars)[field].sum()
        report = report[report > min_val]
        combined_report[cell_id] = report
    combined_report = pd.DataFrame(combined_report).fillna(0)

    return combined_report


def _condense_table(bac_pcc, bac_pct, cut_threshold, ncbi, sciname2taxid):
    """Condense taxonomy table by grouping low-prevalence taxa."""
    taxid2sciname = {v:k for k,v in sciname2taxid.items()}

    mean_clade_rep = np.mean(bac_pcc, axis = 1)
    to_cut = set(bac_pcc.index[np.where(mean_clade_rep < cut_threshold)[0]])
    
    valid_bac = set(bac_pct.index) - to_cut

    catch_all = np.zeros(bac_pct.shape[1])
    new_table = bac_pct.copy()
    to_drop = to_cut.copy()
    for c in to_cut:
        if c not in new_table.index:
            # no mass was assigned specifically to this clade
            to_drop.remove(c)
            continue
        
        # find the lowest clade on the path from root to this node that is represented in "valid_bac"
        taxid = sciname2taxid[c]
        taxid_path = ncbi.get_lineage(taxid)
        sciname_path = [taxid2sciname[x] for x in taxid_path if x in taxid2sciname]
        for candidate in sciname_path[::-1]:
            if candidate in valid_bac:
                # add mass from c to new 
                new_table.loc[candidate] += new_table.loc[c]
                new_table.loc[c] = 0
                break
        if new_table.loc[c].sum() > 0:
            print('found no parent for', c)
            catch_all += new_table.loc[c]
            new_table.loc[c] = 0
    if np.any(catch_all > 0):
        new_table.loc['Dummy non-human root'] = catch_all
    new_table = new_table.drop(labels=list(to_drop))
    return new_table


def _plot_triple_cell_composition(bacteria_df, pcc, total_classified, figsize=(12,12), dpi=200, fname=None):
    """Plot triple panel cell composition figure."""
    plotdf = bacteria_df
    
    plt.figure(figsize=figsize, dpi=dpi)
    plt.subplot(3, 1, 1)
    xs = np.arange(plotdf.shape[1])
    
    my_bottom = 0
    for i in range(len(plotdf)):
        my_height = plotdf.iloc[i]
    
        plt.bar(xs, height=my_height, bottom=my_bottom, label=plotdf.iloc[i].name, facecolor = plt.get_cmap('tab20')(i%20))
        my_bottom += my_height
    
    plt.legend()
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1,1))
    plt.xlabel("")
    plt.ylabel("Proportion of non-human reads")
    
    plotdf = bacteria_df.copy()
    plotdf = plotdf * (pcc.loc['Bacteria'] / 100)
    
    plt.subplot(3, 1, 2)
    xs = np.arange(plotdf.shape[1])
    my_bottom = 0
    for i in range(len(plotdf)):
        my_height = plotdf.iloc[i]
    
        plt.bar(xs, height=my_height, bottom=my_bottom, label=plotdf.iloc[i].name, facecolor = plt.get_cmap('tab20')(i%20))
        my_bottom += my_height
    
    plt.xlabel("")
    plt.ylabel("Proportion of classified reads")

    plt.subplot(3, 1, 3)
    xs = np.arange(plotdf.shape[1])
    my_bottom = 0
    for i in range(len(plotdf)):
        my_height = plotdf.iloc[i] * total_classified.loc[plotdf.columns]
    
        plt.bar(xs, height=my_height, bottom=my_bottom, label=plotdf.iloc[i].name, facecolor = plt.get_cmap('tab20')(i%20))
        my_bottom += my_height
    
    plt.xlabel("Cell")
    plt.ylabel("Number of classified reads")

    ## identify which cells to include
    # always include the first cell in each row
    celldf = plotdf.columns.to_series().str.split('-', expand=True).reset_index(names=['cell_id'])
    celldf = celldf.rename(columns={celldf.columns[-2]:'row', celldf.columns[-1]:'col'})
    celldf['include'] = np.concatenate([[True], celldf.row.iloc[:-1].values != celldf.row.iloc[1:].values])
    
    # also include every other cell s.t. no two adjacent cells are shown
    # i.e., flip the middle entry of each 'False-False-False' triplet
    for i in range(2, len(celldf)):
        if not celldf.include.iloc[i-2:i+1].any():
            celldf.iloc[i-1, -1] = True
    plt.xticks(celldf[celldf.include].index, celldf[celldf.include].cell_id, rotation=90, fontsize='x-small')

    if fname is not None:
        plt.savefig(fname, dpi=dpi, bbox_inches='tight')
