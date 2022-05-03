import os

import mondrianutils.alignment.coverage_metrics as coverage_metrics
import mondrianutils.helpers as helpers

from . import utils as utils


def test_filtering_unpaired(tmpdir):
    sam = utils.get_all_fail_combinations_sam()
    samfile = os.path.join(tmpdir, 'out.sam')
    with open(samfile, 'wt') as writer:
        writer.write(sam)

    bamfile = os.path.join(tmpdir, 'out.bam')
    helpers.run_cmd(['samtools', 'view', '-bSh', samfile], output=bamfile)
    helpers.run_cmd(['samtools', 'index', bamfile])

    with coverage_metrics.CoverageMetrics(
            bamfile, filter_unpaired=True, filter_duplicates=True,
            filter_supplementary=True, filter_secondary=True
    ) as cov:
        for read in cov._bam_reader.fetch():
            if cov._filter_reads(read) is not True:
                assert False, str(read)


def test_filtering_duplicate(tmpdir):
    sam = utils.get_duplicate_sam()
    samfile = os.path.join(tmpdir, 'out.sam')
    with open(samfile, 'wt') as writer:
        writer.write(sam)

    bamfile = os.path.join(tmpdir, 'out.bam')
    helpers.run_cmd(['samtools', 'view', '-bSh', samfile], output=bamfile)
    helpers.run_cmd(['samtools', 'index', bamfile])

    with coverage_metrics.CoverageMetrics(
            bamfile, filter_duplicates=True
    ) as cov:
        for read in cov._bam_reader.fetch():
            if cov._filter_reads(read) is not True:
                assert False, str(read)


def test_filtering_secondary(tmpdir):
    sam = utils.get_secondary_sam()
    samfile = os.path.join(tmpdir, 'out.sam')
    with open(samfile, 'wt') as writer:
        writer.write(sam)

    bamfile = os.path.join(tmpdir, 'out.bam')
    helpers.run_cmd(['samtools', 'view', '-bSh', samfile], output=bamfile)
    helpers.run_cmd(['samtools', 'index', bamfile])

    with coverage_metrics.CoverageMetrics(
            bamfile, filter_secondary=True
    ) as cov:
        for read in cov._bam_reader.fetch():
            if cov._filter_reads(read) is not True:
                assert False, str(read)


def test_filtering_supplementary(tmpdir):
    sam = utils.get_supplementary_sam()
    samfile = os.path.join(tmpdir, 'out.sam')
    with open(samfile, 'wt') as writer:
        writer.write(sam)

    bamfile = os.path.join(tmpdir, 'out.bam')
    helpers.run_cmd(['samtools', 'view', '-bSh', samfile], output=bamfile)
    helpers.run_cmd(['samtools', 'index', bamfile])

    with coverage_metrics.CoverageMetrics(
            bamfile, filter_supplementary=True
    ) as cov:
        for read in cov._bam_reader.fetch():
            if cov._filter_reads(read) is not True:
                assert False, str(read)


def test_filtering_unpaired(tmpdir):
    sam = utils.get_unpaired_sam()
    samfile = os.path.join(tmpdir, 'out.sam')
    with open(samfile, 'wt') as writer:
        writer.write(sam)

    bamfile = os.path.join(tmpdir, 'out.bam')
    helpers.run_cmd(['samtools', 'view', '-bSh', samfile], output=bamfile)
    helpers.run_cmd(['samtools', 'index', bamfile])

    with coverage_metrics.CoverageMetrics(
            bamfile, filter_unpaired=True
    ) as cov:
        for read in cov._bam_reader.fetch():
            if cov._filter_reads(read) is not True:
                assert False, str(read)
