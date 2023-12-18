import csverve.api as csverve
from mondrianutils.dtypes.hmmcopy import dtypes as hmmcopy_dtypes
from mondrianutils.dtypes import haplotypes


def rewrite_csv(infile, outfile, dtypes):
    if dtypes == 'hmmcopy_reads':
        dtypes = hmmcopy_dtypes()['reads']
    elif dtypes == 'hmmcopy_metrics':
        dtypes = hmmcopy_dtypes()['metrics']
    elif dtypes == 'hmmcopy_params':
        dtypes = hmmcopy_dtypes()['params']
    elif dtypes == 'hmmcopy_segs':
        dtypes = hmmcopy_dtypes()['segs']
    elif dtypes == 'haplotypes':
        dtypes = haplotypes.dtypes()['haplotypes']
    else:
        raise Exception()

    csverve.rewrite_csv_file(infile, outfile, dtypes=dtypes)
