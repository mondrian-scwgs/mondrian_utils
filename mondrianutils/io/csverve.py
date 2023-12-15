import click
import csverve.api as csverve
from mondrianutils.dtypes import hmmcopy_metrics
from mondrianutils.dtypes import hmmcopy_params
from mondrianutils.dtypes import hmmcopy_reads
from mondrianutils.dtypes import hmmcopy_segs
from mondrianutils.dtypes import haplotypes


def rewrite_csv(infile, outfile, dtypes):
    if dtypes == 'hmmcopy_reads':
        dtypes = hmmcopy_reads.dtypes()
    elif dtypes == 'hmmcopy_metrics':
        dtypes = hmmcopy_metrics.dtypes()
    elif dtypes == 'hmmcopy_params':
        dtypes = hmmcopy_params.dtypes()
    elif dtypes == 'hmmcopy_segs':
        dtypes = hmmcopy_segs.dtypes()
    elif dtypes == 'haplotypes':
        dtypes = haplotypes.dtypes()
    else:
        raise Exception()

    csverve.rewrite_csv_file(infile, outfile, dtypes=dtypes)



@click.group()
def cli():
    pass

@cli.command()
@click.option('--infile', required=True, help='Path to the input CSV file')
@click.option('--dtypes', required=True, help='Path to the file containing data types')
@click.option('--outfile', required=True, help='Path to the output CSV file')
def rewrite_csv_cmd(infile, dtypes, outfile):
    rewrite_csv(infile, outfile, dtypes)

if __name__ == '__main__':
    cli()
