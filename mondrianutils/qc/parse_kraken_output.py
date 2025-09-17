import click
import os
import re
import pandas as pd

@click.command()
@click.argument('results_dir')
@click.argument('cell_id')
@click.argument('output_table')
@click.argument('output_human')
@click.argument('output_nonhuman')
def main(
        results_dir,
        cell_id,
        output_table,
        output_human,
        output_nonhuman,
    ):
    df = pd.read_table(os.path.join(results_dir, cell_id, cell_id + '_output.txt'), 
                        names = ['is_classified', 'qname', 'result', 'lengths', 'kmers'])

    my_re = re.compile('([A-Za-z0-9 -_\.]+) \(taxid ([0-9]+)\)')

    taxon = []
    taxon_id = []
    for t in df.result:
        match = my_re.search(t)
        taxon.append(match[1])
        taxon_id.append(int(match[2]))

    df['taxon'] = taxon
    df['taxon_id'] = taxon_id

    human_reads = df[df.taxon_id == 9606]
    nonhuman_reads =  df[df.taxon_id != 9606]

    df.to_csv(output_table, index=False)
    human_reads.qname.to_csv(output_human, index=False, header=False)
    nonhuman_reads.qname.to_csv(output_nonhuman, index=False, header=False)

if __name__ == "__main__":
    main()
