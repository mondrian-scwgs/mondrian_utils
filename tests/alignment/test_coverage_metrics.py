import os
from itertools import product
from unittest.mock import MagicMock

import mondrianutils.helpers as helpers
import pysam
from mondrianutils.alignment.coverage_metrics import CoverageMetrics


class CoverageTestData(object):

    def __init__(self, output_filename):
        self.output = output_filename
        self._read_idx = 0

    def __enter__(self):
        self.writer = open(self.output, 'wt')
        self.writer.write(self._get_bam_header())
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.writer.close()

    def _get_bam_header(self):
        header = [
            "@HD\tVN:1.3\tSO:coordinate",
            "@SQ\tSN:20\tLN:63025520",
            "@SQ\tSN:21\tLN:48129895",
            "@SQ\tSN:22\tLN:51304566",
            "@RG\tID:aligned.bam_SA1090_SA1090-A96213A-R22-C43\tCN:FL001\tLB:SA1090\tPL:ILLUMINA\tSM:aligned.bam",
            "@CO\tCB:SA1090-A96213A-R22-C43"
        ]
        header = '\n'.join(header) + '\n'
        return header

    def _get_bam_flag(
            self,
            unpaired=False,
            unmapped=False,
            secondary=False,
            duplicate=False,
            supplementary=False
    ):
        # good read, return flag for paired,mapped,firt_in_pair
        if not any((unpaired, unmapped, supplementary, secondary, duplicate)):
            return 67

        # 1 is flag for paired
        code = 0 if unpaired else 1
        code += 4 if unmapped else 0
        code += 256 if secondary else 0
        code += 1024 if duplicate else 0
        code += 2048 if supplementary else 0

        return code

    def _get_read(self, flag, mapq=60):
        read1 = [
            'HISEQ101_{}:5:1206:10358:49496'.format(self._read_idx),
            flag, '20', '32014074', mapq, '27M1D121M', '=', '32014080', '157',
            'GCGAAACCCCGTCTCTACTAAAATTACAAAAAATTAGCCAGGTGTGGTGGCACACGCCTGTAATCCCAGCTA'
            'CTTGGGAGACTGAGGCAGGAGAATGGCTTGAACCCAGGAGGCAGAGGTTGCAGTGAGCCAAGATTGCGCCACTGCA',
            'AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ'
            'JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJ',
            'CB:Z:SA1090-A96213A-R22-C43', 'MC:Z:21M1D129M', 'MD:Z:1T7T13A3^A53G15C8T26G3C2C8', 'PG:Z:MarkDuplicates',
            'RG:Z:aligned.bam_SA1090_SA1090-A96213A-R22-C43', 'NM:i:10 AS:i:99 FS:Z:grch37_2,mm10_2,salmon_0', 'XS:i:94'
        ]
        read1 = '\t'.join([str(v) for v in read1]) + '\n'

        read2 = [
            'HISEQ101_{}:5:1206:10358:49496'.format(self._read_idx),
            flag, '20', '32014074', mapq, '27M1D121M', '=', '32014080', '157',
            'GCGAAACCCCGTCTCTACTAAAATTACAAAAAATTAGCCAGGTGTGGTGGCACACGCCTGTAATCCCAGCTACT'
            'TGGGAGACTGAGGCAGGAGAATGGCTTGAACCCAGGAGGCAGAGGTTGCAGTGAGCCAAGATTGCGCCACTGCA',
            'AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ'
            'JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJ',
            'CB:Z:SA1090-A96213A-R20-C28', 'MC:Z:21M1D129M', 'MD:Z:1T7T13A3^A53G15C8T26G3C2C8',
            'PG:Z:MarkDuplicates.1', 'RG:Z:aligned.bam_SA1090_SA1090-A96213A-R20-C28',
            'NM:i:10 AS:i:99', 'FS:Z:grch37_2,mm10_2,salmon_0', 'XS:i:94'
        ]
        read2 = '\t'.join([str(v) for v in read2]) + '\n'

        self._read_idx += 1

        return read1 + read2

    def add_all_fail_combinations(self):
        flag_combos = [v for v in product([True, False], repeat=5)]
        flag_combos = [v for v in flag_combos if any(v)]
        flags = [self._get_bam_flag(*v) for v in flag_combos]

        reads = [self._get_read(flag, read_idx=i) for i, flag in enumerate(flags)]
        reads = ''.join(reads)

        header = self._get_bam_header()

        sam = header + reads

        self.writer.write(sam)

    def add_read_by_flag(
            self,
            unpaired=False,
            unmapped=False,
            secondary=False,
            duplicate=False,
            supplementary=False
    ):
        flag = self._get_bam_flag(
            unpaired=unpaired,
            unmapped=unmapped,
            secondary=secondary,
            duplicate=duplicate,
            supplementary=supplementary
        )

        sam = self._get_bam_header() + self._get_read(flag)

        self.writer.write(sam)


def test_get_bam_reader(tmpdir):
    samfile = os.path.join(tmpdir, 'out.sam')
    with CoverageTestData(samfile) as testdata:
        testdata.add_read_by_flag(unpaired=True)

    coverage_metrics = CoverageMetrics(samfile)
    bam_reader = coverage_metrics._get_bam_reader()

    assert isinstance(bam_reader, pysam.AlignmentFile)
    assert bam_reader.filename.decode('utf-8') == samfile
    assert bam_reader.mode.decode('utf-8') == 'rb'


def test_genome_length():
    samfile = "path/to/bamfile.bam"
    header = {'SQ': [{'LN': 100}, {'LN': 200}, {'LN': 150}]}
    bam_reader = MagicMock()
    bam_reader.header = header
    coverage_metrics = CoverageMetrics(samfile)
    coverage_metrics._bam_reader = bam_reader

    assert coverage_metrics.genome_length == 450


def test_get_read_intervals():
    samfile = "path/to/bamfile.bam"
    read = MagicMock()
    read.reference_start = 10
    read.reference_end = 20
    read.query_alignment_qualities = [30, 40, 50, 60, 70]
    read.is_reverse = False

    coverage_metrics = CoverageMetrics(samfile, min_base_qual=40)
    regions = coverage_metrics._get_read_intervals(read)

    assert regions == [(11, 14)]


def test_filtering_unpaired(tmpdir):
    samfile = os.path.join(tmpdir, 'out.sam')
    with CoverageTestData(samfile) as testdata:
        testdata.add_all_fail_combinations()

    bamfile = os.path.join(tmpdir, 'out.bam')
    helpers.run_cmd(['samtools', 'view', '-bSh', samfile], output=bamfile)
    helpers.run_cmd(['samtools', 'index', bamfile])

    with CoverageMetrics(
            bamfile, filter_unpaired=True, filter_duplicates=True,
            filter_supplementary=True, filter_secondary=True
    ) as cov:
        for read in cov._bam_reader.fetch():
            if cov._filter_reads(read) is not True:
                assert False, str(read)


def test_filtering_duplicate(tmpdir):
    samfile = os.path.join(tmpdir, 'out.sam')
    with CoverageTestData(samfile) as testdata:
        testdata.add_read_by_flag(duplicate=True)

    bamfile = os.path.join(tmpdir, 'out.bam')
    helpers.run_cmd(['samtools', 'view', '-bSh', samfile], output=bamfile)
    helpers.run_cmd(['samtools', 'index', bamfile])

    with CoverageMetrics(
            bamfile, filter_duplicates=True
    ) as cov:
        for read in cov._bam_reader.fetch():
            if cov._filter_reads(read) is not True:
                assert False, str(read)


def test_filtering_secondary(tmpdir):
    samfile = os.path.join(tmpdir, 'out.sam')
    with CoverageTestData(samfile) as testdata:
        testdata.add_read_by_flag(secondary=True)

    bamfile = os.path.join(tmpdir, 'out.bam')
    helpers.run_cmd(['samtools', 'view', '-bSh', samfile], output=bamfile)
    helpers.run_cmd(['samtools', 'index', bamfile])

    with CoverageMetrics(
            bamfile, filter_secondary=True
    ) as cov:
        for read in cov._bam_reader.fetch():
            if cov._filter_reads(read) is not True:
                assert False, str(read)


def test_filtering_supplementary(tmpdir):
    samfile = os.path.join(tmpdir, 'out.sam')
    with CoverageTestData(samfile) as testdata:
        testdata.add_read_by_flag(supplementary=True)

    bamfile = os.path.join(tmpdir, 'out.bam')
    helpers.run_cmd(['samtools', 'view', '-bSh', samfile], output=bamfile)
    helpers.run_cmd(['samtools', 'index', bamfile])

    with CoverageMetrics(
            bamfile, filter_supplementary=True
    ) as cov:
        for read in cov._bam_reader.fetch():
            if cov._filter_reads(read) is not True:
                assert False, str(read)


def test_filtering_unpaired(tmpdir):
    samfile = os.path.join(tmpdir, 'out.sam')
    with CoverageTestData(samfile) as testdata:
        testdata.add_read_by_flag(unpaired=True)

    bamfile = os.path.join(tmpdir, 'out.bam')
    helpers.run_cmd(['samtools', 'view', '-bSh', samfile], output=bamfile)
    helpers.run_cmd(['samtools', 'index', bamfile])

    with CoverageMetrics(
            bamfile, filter_unpaired=True
    ) as cov:
        for read in cov._bam_reader.fetch():
            if cov._filter_reads(read) is not True:
                assert False, str(read)


def test_get_read_intervals_base_quality_threshold(tmpdir):
    samfile = os.path.join(tmpdir, 'out.sam')
    with CoverageTestData(samfile) as testdata:
        testdata.add_read_by_flag()

    with CoverageMetrics(samfile, min_base_qual=10) as covmet:
        read = MagicMock(
            reference_start=0,
            reference_end=10,
            query_alignment_qualities=[30, 20, 10, 5, 5, 10, 20, 30, 40, 50, 50],
            is_reverse=False
        )
        regions = covmet._get_read_intervals(read)

        assert regions == [(0, 2), (5, 10)]

    with CoverageMetrics(samfile, min_base_qual=0) as covmet:
        read = MagicMock(
            reference_start=0,
            reference_end=10,
            query_alignment_qualities=[30, 20, 10, 5, 5, 10, 20, 30, 40, 50, 50],
            is_reverse=False
        )
        regions = covmet._get_read_intervals(read)

        assert regions == [(0, 10)]


def test_get_read_intervals_base_quality_threshold_reversed(tmpdir):
    samfile = os.path.join(tmpdir, 'out.sam')
    with CoverageTestData(samfile) as testdata:
        testdata.add_read_by_flag()

    with CoverageMetrics(samfile, min_base_qual=10) as covmet:
        read = MagicMock(
            reference_start=0,
            reference_end=10,
            query_alignment_qualities=[30, 20, 10, 5, 5, 10, 20, 30, 40, 50, 50],
            is_reverse=True
        )
        regions = covmet._get_read_intervals(read)

        assert regions == [(0, 5), (8, 10)]

    with CoverageMetrics(samfile, min_base_qual=0) as covmet:
        read = MagicMock(
            reference_start=0,
            reference_end=10,
            query_alignment_qualities=[30, 20, 10, 5, 5, 10, 20, 30, 40, 50, 50],
            is_reverse=True
        )
        regions = covmet._get_read_intervals(read)

        assert regions == [(0, 10)]


def test_merge_overlapping_intervals(tmpdir):
    samfile = os.path.join(tmpdir, 'out.sam')
    with CoverageTestData(samfile) as testdata:
        testdata.add_read_by_flag()

    with CoverageMetrics(samfile, min_base_qual=10) as covmet:
        merged = covmet._merge_overlapping_intervals([[0, 10], [20, 30]])
        assert merged == [[0, 10], [20, 30]]

        merged = covmet._merge_overlapping_intervals([[0, 10], [5, 15], [20, 30]])
        assert merged == [[0, 15], [20, 30]]

        merged = covmet._merge_overlapping_intervals([[5, 15], [0, 10], [20, 30]])
        assert merged == [[0, 15], [20, 30]]


def test_merge_overlapping_intervals(tmpdir):
    samfile = os.path.join(tmpdir, 'out.sam')
    with CoverageTestData(samfile) as testdata:
        testdata.add_read_by_flag()

    with CoverageMetrics(samfile, min_base_qual=10) as covmet:
        covmet._bam_reader = MagicMock(
            header={'SQ': [{'LN': 1000}]}
        )

        cov = covmet.get_coverage(
            {'RID1': {'1': [[1, 10], [5, 15], [20, 30]]}}
        )
        assert cov == 0.024

        cov = covmet.get_coverage(
            {'RID1': {'1': [[5, 15], [1, 10], [20, 30]]}}
        )
        assert cov == 0.024

        cov = covmet.get_coverage(
            {'RID1': {'1': [[1, 10], [20, 30]]}}
        )
        assert cov == 0.019
