import click
import pandas as pd
import numpy as np
import pickle

def skip_ahead(lines, i, skip_fields):
    while  i < len(lines) and (lines[i].startswith('#') or lines[i].split('\t')[0] in skip_fields):
        i += 1
    return i

def parse_bamstat(fname):
    skip_fields = set(['#', 'CHK', 'GCT', 'GCC', 'GCT', 'FBC', 'FTC', 'LBC', 'LTC', 'BCC', 'CRC', 'OXC', 'RXC', 'GCD'])
    lines = open(fname).readlines()
    result = {}

    summary = {}
    i = skip_ahead(lines, 0, skip_fields)
    
    # summary fields
    while lines[i].startswith('SN'):
        tkns = lines[i].split('\t')
        summary[tkns[1]] = float(tkns[2].strip())
        i += 1
    result['summary'] = summary

    i = skip_ahead(lines, i, skip_fields)
    if summary['reads mapped and paired:'] == 0:
        return result
    
    ffq = []
    while lines[i].startswith('FFQ'):
        tkns = lines[i].strip().split('\t')
        ffq.append(np.array(tkns[2:]).astype(int))
        i += 1
    ffq = np.array(ffq)
    result['first_fragment_quality'] = ffq
    
    i = skip_ahead(lines, i, skip_fields)
    
    lfq = []
    while lines[i].startswith('LFQ'):
        tkns = lines[i].strip().split('\t')
        lfq.append(np.array(tkns[2:]).astype(int))
        i += 1
    lfq = np.array(lfq)
    result['last_fragment_quality'] = lfq
    
    i = skip_ahead(lines, i, skip_fields)
    
    gcf = []
    while lines[i].startswith('GCF'):
        tkns = lines[i].strip().split('\t')
        gcf.append(np.array(tkns[1:]))
        i += 1
    gcf = np.array(gcf)
    result['first_fragment_gc'] = gcf

    i = skip_ahead(lines, i, skip_fields)
    
    gcl = []
    while lines[i].startswith('GCL'):
        tkns = lines[i].strip().split('\t')
        gcl.append(np.array(tkns[1:]))
        i += 1
    gcl = np.array(gcl)
    result['last_fragment_gc'] = gcl
   
    # skip the next few fields
    i = skip_ahead(lines, i, skip_fields)
    
    insert_sizes = []
    while lines[i].startswith('IS'):
        tkns = lines[i].strip().split('\t')
        insert_sizes.append({'insert_size':int(tkns[1]),
                             'pairs_total':int(tkns[2]),
                             'inward_oriented_pairs':int(tkns[3]),
                             'outward_oriented_pairs':int(tkns[4]),
                             'other_pairs':int(tkns[5])
                            })
        i += 1
    insert_sizes = pd.DataFrame(insert_sizes)
    result['insert_sizes'] = insert_sizes
    
    i = skip_ahead(lines, i, skip_fields)
    
    rl = []
    while lines[i].startswith('RL'):
        tkns = lines[i].strip().split('\t')
        rl.append({'read_length':int(tkns[1]),
                             'count':int(tkns[2]),
                            })
        i += 1
    rl = pd.DataFrame(rl)
    result['read_lengths'] = rl

    i = skip_ahead(lines, i, skip_fields)
    
    frl = []
    while lines[i].startswith('FRL'):
        tkns = lines[i].strip().split('\t')
        frl.append({'fragment_length':int(tkns[1]),
                             'count':int(tkns[2]),
                            })
        i += 1
    frl = pd.DataFrame(frl)
    result['first_fragment_lengths'] = frl

    i = skip_ahead(lines, i, skip_fields)
    
    lrl = []
    while lines[i].startswith('LRL'):
        tkns = lines[i].strip().split('\t')
        lrl.append({'fragment_length':int(tkns[1]),
                             'count':int(tkns[2]),
                            })
        i += 1
    lrl = pd.DataFrame(lrl)
    result['last_fragment_lengths'] = lrl

    i = skip_ahead(lines, i, skip_fields)
    
    mapq = []
    while lines[i].startswith('MAPQ'):
        tkns = lines[i].strip().split('\t')
        mapq.append({'mapq':int(tkns[1]),
                             'count':int(tkns[2]),
                            })
        i += 1
    mapq = pd.DataFrame(mapq)
    result['mapq'] = mapq

    i = skip_ahead(lines, i, skip_fields)
    
    indels = []
    while lines[i].startswith('ID'):
        tkns = lines[i].strip().split('\t')
        indels.append({'length':int(tkns[1]),
                             'n_insertions':int(tkns[2]),
                             'n_deletions':int(tkns[3])
                            })
        i += 1
    indels = pd.DataFrame(indels)
    result['indels'] = indels

    i = skip_ahead(lines, i, skip_fields)
    
    indels_per_cycle = []
    while lines[i].startswith('IC'):
        tkns = lines[i].strip().split('\t')
        indels_per_cycle.append({'cycle':int(tkns[1]),
                             'n_insertions_fwd':int(tkns[2]),
                             'n_insertions_rev':int(tkns[3]),
                             'n_deletions_fwd':int(tkns[4]),
                             'n_deletions_rev':int(tkns[5]),
                            })
        i += 1
    indels_per_cycle = pd.DataFrame(indels_per_cycle)
    result['indels_per_cycle'] = indels_per_cycle

    
    i = skip_ahead(lines, i, skip_fields)
    
    cov = []
    while lines[i].startswith('COV'):
        tkns = lines[i].strip().split('\t')
        cov.append({'cov':int(tkns[2]),
                             'n_bases':int(tkns[3]),
                            })
        i += 1
    cov = pd.DataFrame(cov)
    result['coverage'] = cov

    return result

@click.command()
@click.argument('stats_file')
@click.argument('output_file')
def main(
        stats_file,
        output_file,
    
    ):
    result = parse_bamstat(stats_file)
    with open(output_file, 'wb') as f:
        pickle.dump(result, f)


if __name__ == "__main__":
    main()
