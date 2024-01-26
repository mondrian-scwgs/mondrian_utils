import os

import csverve.api as csverve
import mondrianutils.helpers as helpers
from mondrianutils.hmmcopy.correct_read_count import CorrectReadCount
from mondrianutils.hmmcopy.plot_hmmcopy import plot_hmmcopy

from mondrianutils.hmmcopy import add_mappability
from mondrianutils.hmmcopy import add_quality

from mondrianutils.hmmcopy.complete_hmmcopy import run_hmmcopy
from mondrianutils.hmmcopy import readcounter
from glob import glob

from mondrianutils.io import overlapping_fraction_per_bin
import pysam


def cell_hmmcopy(
        cellbam, gc_wig_file, map_wig_file, alignment_metrics, exclude_list,
        chromosomes, metrics, params, reads, segments,
        output_tarball, reference, segments_output, bias_output, tempdir,
        quality_classifier_training_data, quality_classifier_model=None,
        mappability_cutoff=0.9, binsize=500000, mapping_quality_threshold=20
):
    cell_id = helpers.get_cells(pysam.AlignmentFile(cellbam, 'rb'))
    assert len(cell_id) == 1, cell_id
    cell_id = cell_id[0]

    helpers.makedirs(tempdir)

    readcount_dir = os.path.join(tempdir, 'readcounter')
    readcount_tmpdir = os.path.join(tempdir, 'readcounter_tmp')

    print('running readcounter')
    readcounter(
        cellbam,
        readcount_dir,
        readcount_tmpdir,
        chromosomes,
        exclude_list=exclude_list,
        ncores=1,
        mapping_quality_threshold=mapping_quality_threshold,
        window_size=binsize
    )
    readcount_wig = glob(f"{readcount_dir}/*.wig")
    assert len(readcount_wig) == 1
    readcount_wig = readcount_wig[0]

    print('running correction')
    corrected_reads = os.path.join(tempdir, 'corrected_reads.csv')
    correction = CorrectReadCount(
        gc_wig_file, map_wig_file, readcount_wig, corrected_reads,
        cell_id, mappability=mappability_cutoff
    )
    correction.main()

    hmmcopy_tempdir = os.path.join(tempdir, 'hmmcopy')
    helpers.makedirs(hmmcopy_tempdir)
    raw_reads = os.path.join(hmmcopy_tempdir, 'reads.csv.gz')
    raw_metrics = os.path.join(hmmcopy_tempdir, 'metrics.csv.gz')
    print('running hmmcopy')
    run_hmmcopy(
        corrected_reads, hmmcopy_tempdir, raw_metrics,
        params, raw_reads, segments, output_tarball
    )

    mappability_reads = os.path.join(hmmcopy_tempdir, 'mappability_reads.csv.gz')
    add_mappability(raw_reads, mappability_reads)

    overlapping_fraction = os.path.join(hmmcopy_tempdir, 'overlapping_fraction.csv.gz')
    overlapping_fraction_tmpdir = os.path.join(hmmcopy_tempdir, 'overlapping_fraction')
    print('running overlapping fraction')
    overlapping_fraction_per_bin(
        cellbam, overlapping_fraction, overlapping_fraction_tmpdir,
        chromosomes=chromosomes,
        binsize=binsize,
        mapping_quality=mapping_quality_threshold, ncores=1
    )
    print('merge fraction and reads')
    csverve.merge_csv(
        [mappability_reads, overlapping_fraction], reads, how='outer',
        on=['chr', 'start', 'end', 'cell_id']
    )

    print('plot hmmcopy')
    plot_hmmcopy(
        reads, segments, params, raw_metrics, reference,
        segments_output, bias_output
    )

    merged_metrics = os.path.join(hmmcopy_tempdir, 'alignment_merged_metrics.csv.gz')
    csverve.merge_csv(
        [raw_metrics, alignment_metrics], merged_metrics, how='inner', on='cell_id', write_header=True, lenient=True
    )

    print('add quality')
    add_quality(merged_metrics, metrics, quality_classifier_training_data,
                joblib_model=quality_classifier_model)
