'''
Created on Oct 10, 2017

@author: dgrewal
'''
import glob
import os
from collections import defaultdict

import argparse
import numpy as np
import pandas as pd
import pysam
import tqdm


def create_bincounts():
    return defaultdict(int)


def create_cellcounts():
    return defaultdict(create_bincounts)


class CellReadCounts:
    """
    Store and manage read counts per bin for each cell and chromosome.

    :param chromosomes: List of chromosome names
    :type chromosomes: list of str
    :param chr_lengths: Chromosome lengths keyed by chromosome name
    :type chr_lengths: dict
    :param window_size: Bin size used to count reads
    :type window_size: int
    """
    
    def __init__(self, chromosomes, chr_lengths, window_size):
        self.chromosomes = chromosomes
        self.chr_lengths = chr_lengths
        self.window_size = window_size
        self.counts = defaultdict(create_cellcounts)

    def get_observed_cells(self):
        """
        Get a sorted list of all observed cell IDs.

        :return: List of cell IDs
        :rtype: list of str
        """
        cells = set()
        for chrom, data in self.counts.items():
            for cell in data:
                cells.add(cell)
        return list(sorted(cells))

    def get_all_bins(self, chrom):
        """
        Get all bins for a given chromosome.

        :param chrom: Chromosome name
        :type chrom: str
        :return: List of bin tuples (start, end)
        :rtype: list of tuple
        """
        reflen = self.chr_lengths[chrom]

        start = 0
        end = start + self.window_size
        bins = [(start, end)]
        while end < reflen:
            start += self.window_size
            end = min(start + self.window_size, reflen)
            bins.append((start, end))
        return bins

    def write_wig_files(self, output_dir, cells=None):
        """
        Write wiggle format output files for each cell.

        :param output_dir: Directory to write .wig files
        :type output_dir: str
        :param cells: Optional list of cell IDs to write; all if None
        :type cells: list of str or None
        """
        if cells is None:
            cells = self.get_observed_cells()

        for cell in cells:
            with open(os.path.join(output_dir, '{}.wig'.format(cell)), 'wt') as outfile:

                header_str = "track type=wiggle_0 name={}\n".format(cell)
                outfile.write(header_str)

                for chrom in self.chromosomes:
                    track_header_str = "fixedStep chrom={} start=1 step={} span={}\n".format(
                        chrom, self.window_size, self.window_size)
                    outfile.write(track_header_str)

                    bins = self.get_all_bins(chrom)
                    for bin in bins:
                        outfile.write(str(self.counts[chrom][cell][bin]) + '\n')

    def write_dataframe(self, output_filename, cells=None):
        """
        Generate and write read counts as a pandas DataFrame.

        :param output_filename: Output file path for the dataframe
        :type output_filename: str
        :param cells: Optional list of cell IDs to include
        :type cells: list of str or None
        """
        if cells is None:
            cells = self.get_observed_cells()

        data = []
        for cell in cells:
            cell_data = []
            for chrom in self.chromosomes:
                bins = self.get_all_bins(chrom)
                for bin in bins:
                    start, end = bin
                    cell_data.append({
                        'cell_id': cell,
                        'chromosome': chrom,
                        'start': start,
                        'end': end,
                        'count': self.counts[chrom][cell][bin],
                    })
            
            cell_data = pd.DataFrame(cell_data)
            cell_data['cell_id'] = pd.Categorical(cell_data['cell_id'], categories=cells)
            cell_data['chromosome'] = pd.Categorical(cell_data['chromosome'], categories=self.chromosomes)

            data.append(cell_data)

        if len(data) == 0:
            data = pd.DataFrame(columns=['cell_id', 'chromosome', 'start', 'end', 'count'])
        else:
            data = pd.concat(data, ignore_index=True)

        data.to_csv(output_filename, index=False)

    def get_overlapping_bin(self, pos, chrom):
        """
        Get the bin (start, end) that overlaps with a given position.

        :param pos: Genomic position
        :type pos: int
        :param chrom: Chromosome name
        :type chrom: str
        :return: Bin as (start, end)
        :rtype: tuple of int
        """

        ## TODO: refactor.
        # if pos = 100 and binsize is 10. put it in 90,100. not 100,110
        if pos > self.window_size and pos % self.window_size == 0:
            end = (pos // self.window_size) * self.window_size
            start = end - self.window_size
        else:
            start = (pos // self.window_size) * self.window_size
            end = start + self.window_size

        if end > self.chr_lengths[chrom]:
            end = self.chr_lengths[chrom]

        return (start, end)

    def add_read(self, chrom, cell_id, pos):
        """
        Add a read to the count for a given cell, chromosome, and position.

        :param chrom: Chromosome name
        :type chrom: str
        :param cell_id: Cell barcode/tag
        :type cell_id: str
        :param pos: Genomic position of the read
        :type pos: int
        """
        binval = self.get_overlapping_bin(pos, chrom)
        self.counts[chrom][cell_id][binval] += 1


def merge_cell_read_counts(counts_list, chromosomes):
    """
    Merge read count data from multiple CellReadCounts objects.

    :param counts_list: List of CellReadCounts instances
    :type counts_list: list of CellReadCounts
    :param chromosomes: Chromosomes to include in merged data
    :type chromosomes: list of str
    :return: Merged CellReadCounts object
    :rtype: CellReadCounts
    """
    assert len(counts_list) > 0

    window_size = counts_list[0].window_size
    chr_lengths = counts_list[0].chr_lengths
    for counts in counts_list:
        assert window_size == counts.window_size
        assert set(chr_lengths.items()) == set(counts.chr_lengths.items())
        assert set(counts.chromosomes).issubset(chromosomes)

    merged = CellReadCounts(chromosomes, chr_lengths, window_size)
    for counts in counts_list:
        for chrom, chrom_data in sorted(counts.counts.items()):
            for cell, cell_data in chrom_data.items():
                for bin, count in cell_data.items():
                    merged.counts[chrom][cell][bin] += count
    
    return merged


class ReadCounter:
    """
    Count reads per bin from a BAM file, optionally filtering by mapping quality and exclusion regions.

    :param bam: Path to input BAM file
    :type bam: str
    :param window_size: Bin size to use for counting
    :type window_size: int
    :param mapq: Mapping quality threshold
    :type mapq: int
    :param excluded: Optional path to regions to exclude
    :type excluded: str or None
    :param tag_name: Read tag to use for cell barcodes (default: 'CB')
    :type tag_name: str
    :param ignore_missing_tags: Skip reads without the barcode tag
    :type ignore_missing_tags: bool
    """

    def __init__(
            self, bam, window_size, mapq, excluded=None, tag_name='CB',
            ignore_missing_tags=False,
    ):
        self.bamfile = bam
        self.window_size = window_size
        self.mapq_threshold = mapq

        if excluded is not None:
            self.excluded = pd.read_csv(excluded, sep="\t", )
            self.excluded.columns = ["chrom", "start", "end"]
        else:
            self.excluded = None
        
        self.tag_name = tag_name
        self.ignore_missing_tags = ignore_missing_tags

        if self.is_empty():
            self.chr_lengths = None
        else:
            self.chr_lengths = self.__get_chr_lengths()

    def is_empty(self):
        """
        Check if the BAM file is empty or contains 'NO DATA'.

        :return: True if file is empty or contains 'NO DATA'
        :rtype: bool
        """
        if os.path.getsize(self.bamfile) < 10:
            with open(self.bamfile, 'rt') as reader:
                if reader.readline().startswith('NO DATA'):
                    return True

    def __get_chrom_excluded(self, chrom):
        """
        Generate a binary mask of excluded positions for a chromosome.

        :param chrom: Chromosome name
        :type chrom: str
        :return: Numpy array where 1 indicates excluded position
        :rtype: numpy.ndarray
        """
        chrom_length = self.chr_lengths[chrom]
        chrom_excluded = np.zeros(chrom_length + 1, dtype=np.uint8)

        for start, end in self.excluded.loc[self.excluded['chrom'] == chrom, ['start', 'end']].values:
            start = min(start, chrom_length)
            end = min(end, chrom_length)
            chrom_excluded[start:end] = 1

        return chrom_excluded

    def __get_chr_lengths(self):
        """
        Get chromosome lengths from BAM file header.

        :return: Dictionary with chromosome names and lengths
        :rtype: dict
        """
        with pysam.AlignmentFile(self.bamfile, 'rb') as bam:
            names = bam.references
            lengths = bam.lengths
            return {name: length for name, length in zip(names, lengths)}

    def filter(self, pileupobj, chrom_excluded=None):
        """
        Apply filtering to a read: exclude low mapping quality, duplicates, and masked regions.

        :param pileupobj: Pysam read object
        :type pileupobj: pysam.AlignedSegment
        :param chrom_excluded: Optional binary mask of excluded positions
        :type chrom_excluded: numpy.ndarray or None
        :return: True if the read should be excluded
        :rtype: bool
        """
        pos = pileupobj.reference_start
        if chrom_excluded is not None and chrom_excluded[pos]:
            return True

        if pileupobj.is_duplicate:
            return True

        if pileupobj.mapping_quality < self.mapq_threshold:
            return True

        return False

    def get_data(self, chrom):
        """
        Iterate over reads in a chromosome, count reads per bin.

        :param chrom: Chromosome name
        :type chrom: str
        :return: CellReadCounts object with counts for the chromosome
        :rtype: CellReadCounts
        """
        with pysam.AlignmentFile(self.bamfile, 'rb') as bam:

            counts = CellReadCounts([chrom], self.chr_lengths, self.window_size)

            chrom_excluded = None
            if self.excluded is not None:
                chrom_excluded = self.__get_chrom_excluded(chrom)

            for pileupobj in tqdm.tqdm(bam.fetch(chrom)):
                if self.filter(pileupobj, chrom_excluded):
                    continue

                if self.ignore_missing_tags and not pileupobj.has_tag(self.tag_name):
                    continue

                cell_id = pileupobj.get_tag(self.tag_name)

                counts.add_read(chrom, cell_id, pileupobj.pos)
            
            return counts


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('bam',
                        help='specify the path to the input bam file')

    parser.add_argument('outdir',
                        help='specify path to the output dir')

    parser.add_argument('--chromosomes',
                        nargs='*',
                        default=[str(v) for v in range(1, 23)] + ['X', 'Y'],
                        help='specify target chromosomes'
                        )
    parser.add_argument('-w', '--window_size',
                        type=int,
                        default=1000,
                        help='specify bin size')
    parser.add_argument('-m', '--mapping_quality_threshold',
                        type=int,
                        default=0,
                        help='threshold for the mapping quality, reads ' \
                             'with quality lower than threshold will be ignored')

    parser.add_argument('--seg',
                        default=False,
                        action='store_true',
                        help='write the output in seg format')

    parser.add_argument('--exclude_list',
                        default=None,
                        help='regions to skip')

    parser.add_argument('--reference', default=None)

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = parse_args()
    with ReadCounter(args.bam, args.output, args.window_size,
                     args.chromosomes, args.mapping_quality_threshold,
                     args.seg, excluded=args.exclude_list,
                     reference=args.reference) as rcount:
        rcount.main()
