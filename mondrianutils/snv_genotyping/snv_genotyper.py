from collections import defaultdict

import csverve.api as csverve
import mondrianutils.helpers as helpers
import pandas as pd
import pysam
import vcf


class SnvGenotyper(object):
    def __init__(
            self,
            bamfile,
            targets,
            output,
            cell_barcodes=False,
            interval=None,
            count_duplicates=False,
            min_mqual=20,
            sparse=False,
            ignore_untagged_reads=False
    ):
        self.bam = self._get_bam_reader(bamfile)

        self.cell_barcodes = cell_barcodes
        self.all_cells = self.get_cells()

        self.chrom = None
        self.begin = None
        self.end = None
        if interval:
            self.chrom, self.begin, self.end = self._parse_interval(interval)

        self.targets = targets
        self.output = output

        self.count_duplicates = count_duplicates
        self.min_mqual = min_mqual
        self.sparse = sparse
        self.ignore_untagged_reads = ignore_untagged_reads

    @property
    def dtypes(self):
        dtypes = {
            'chrom': 'str',
            'pos': int,
            'ref': 'str',
            'alt': 'str',
            'cell_id': 'str',
            'ref_counts': int,
            'alt_counts': int
        }
        return dtypes

    def __get_bam_header(self):
        return self.bam.header

    def _get_cells_from_header(self):
        header = self.__get_bam_header()
        cells = []
        for line in str(header).split('\n'):
            if not line.startswith("@CO"):
                continue
            line = line.strip().split()
            cb = line[1]
            cell = cb.split(':')[1]
            cells.append(cell)
        return cells

    def _get_cells_from_barcodes(self):
        cells = []
        with helpers.getFileHandle(self.cell_barcodes, 'rt') as reader:
            for line in reader:
                cells.append(line.strip())
            return cells

    def get_cells(self):
        if self.cell_barcodes:
            return self._get_cells_from_barcodes()
        else:
            cells = self._get_cells_from_header()
            if len(cells) == 0:
                raise Exception('No cells ids found in bam header and no cell barcodes file provided')
            return cells

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass
        # clean up output if there are any exceptions
        # if exc_type and os.path.exists(self.output):
        #     os.remove(self.output)

    def _get_bam_reader(self, bamfile):
        """returns pysam bam object
        :returns pysam bam object
        """
        return pysam.AlignmentFile(bamfile, 'rb')

    def _fetch(self, chrom, start, end):
        """returns iterator over reads in the specified region
        :param chrom: chromosome name (str)
        :param start: bin starting pos (int)
        :param end: bin end pos (int)
        :returns iterator over reads
        """
        return self.bam.fetch(chrom, start, end)

    @staticmethod
    def _parse_interval(interval):
        """
        allowed formats:
        1. None
        2. just chrom ex: 1
        3. chrom, start, end. ex: 1:1-1000

        Parameters
        ----------
        region :

        Returns
        -------

        """

        if interval is None:
            return None, None, None

        if ':' not in interval:
            return interval, None, None

        chrom, coords = interval.split(':')

        assert '-' in coords

        beg, end = coords.split('-')

        beg = int(beg) - 1
        end = int(end)

        return chrom, beg, end

    def _fetch_vcf_reader(self, vcf_file):
        vcf_reader = vcf.Reader(filename=vcf_file)

        if self.chrom is None:
            return vcf_reader

        try:
            vcf_reader = vcf_reader.fetch(self.chrom, start=self.begin, end=self.end)

        except ValueError:
            vcf_reader = ()

        return vcf_reader

    def load_targets(self, vcf_file):
        reader = self._fetch_vcf_reader(vcf_file)

        targets = [(record.CHROM, record.POS, record.REF, record.ALT) for record in reader]

        return targets

    def _check_read(self, read):
        valid = True

        if read.alignment.mapping_quality < self.min_mqual:
            valid = False

        elif read.alignment.is_duplicate and (not self.count_duplicates):
            valid = False

        elif read.alignment.is_unmapped:
            valid = False

        elif read.alignment.is_qcfail:
            valid = False

        elif read.alignment.is_secondary:
            valid = False

        elif self.ignore_untagged_reads:
            try:
                read.alignment.get_tag('CB')
            except KeyError:
                valid = False

        return valid

    def _get_counts_pos(self, chrom, pos, ref, alt):
        ref_count = defaultdict(int)
        alt_count = defaultdict(int)

        for pileupcol in self.bam.pileup(
                str(chrom), int(pos) - 200, int(pos) + 200,
                ignore_overlaps=False, max_depth=1e6
        ):
            for pileupread in pileupcol.pileups:

                if not self._check_read(pileupread):
                    continue

                if pileupcol.pos == pos - 1 and pileupread.query_position is not None:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    cell_id = pileupread.alignment.get_tag('CB')

                    if base == ref:
                        ref_count[cell_id] += 1
                    elif base == alt:
                        alt_count[cell_id] += 1

        return ref_count, alt_count

    def get_counts(self, targets):
        data = []

        for (chrom, pos, ref, alts) in targets:

            for alt in alts:
                ref_count, alt_count = self._get_counts_pos(chrom, pos, ref, alt)

                for cell in self.all_cells:
                    row = [chrom, pos, ref, alt, cell, ref_count[cell], alt_count[cell]]

                    if self.sparse:
                        if ref_count[cell] == 0 and alt_count[cell] == 0:
                            continue

                    data.append(row)

        data = pd.DataFrame(data, columns=['chrom', 'pos', 'ref', 'alt', 'cell_id', 'ref_counts', 'alt_counts'])

        return data

    def genotyping(self):
        targets = self.load_targets(self.targets)

        df = self.get_counts(targets)

        csverve.write_dataframe_to_csv_and_yaml(
            df, self.output, self.dtypes
        )
