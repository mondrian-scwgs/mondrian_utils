import csverve.api as csverve
from mondrianutils.dtypes.breakpoint import dtypes as breakpoint_dtypes
import mondrianutils.helpers as helpers
import pandas as pd
import pysam


class SvGenotyper(object):
    def __init__(
            self,
            bamfile,
            destruct_reads,
            destruct_table,
            output
    ):
        self.bam = self._get_bam_reader(bamfile)
        self.cells = self._get_cell_ids()

        self.destruct_reads = destruct_reads
        self.destruct_table = destruct_table
        self.output = output

    @property
    def dtypes(self):
        dtypes = breakpoint_dtypes()['genotyping']
        return dtypes

    def __get_bam_header(self):
        return self.bam.header

    def _get_bam_reader(self, bamfile):
        """returns pysam bam object
        :returns pysam bam object
        """
        return pysam.AlignmentFile(bamfile, 'rb')

    def _get_cell_ids(self):
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

    def _fetch(self):
        """returns iterator over reads in the specified region
        :param chrom: chromosome name (str)
        :param start: bin starting pos (int)
        :param end: bin end pos (int)
        :returns iterator over reads
        """
        return self.bam.fetch()

    def load_reads(self, destruct_reads):
        read_ids = {}
        with helpers.getFileHandle(destruct_reads, 'rt') as reader:
            for line in reader:
                line = line.strip().split()
                prediction_id = line[0]
                read_id = line[-2]

                assert read_id.startswith('+')
                read_id = read_id[1:]

                if read_id not in read_ids:
                    read_ids[read_id] = []
                read_ids[read_id].append(prediction_id)
        return read_ids

    def load_table(self):
        df = pd.read_csv(
            self.destruct_table, sep='\t', dtype=str,
            usecols=['prediction_id', 'chromosome_1', 'strand_1', 'position_1', 'chromosome_2', 'strand_2',
                     'position_2']
        )

        df = df.set_index('prediction_id')
        df = df.to_dict('index')

        return df

    def _get_count_data(self, reads_data):

        outputdata = {}

        for read in self.bam.fetch(until_eof=True):
            read_id = read.query_name

            try:
                cell_id = read.get_tag('CB')
            except:
                cell_id = "UNKNOWN"

            if read_id in reads_data:
                prediction_ids = reads_data[read_id]

                for prediction_id in prediction_ids:

                    if prediction_id not in outputdata:
                        outputdata[prediction_id] = {}
                    if cell_id not in outputdata[prediction_id]:
                        outputdata[prediction_id][cell_id] = 0

                    outputdata[prediction_id][cell_id] += 1

        return outputdata

    def annotate_table(self, df, counts):
        output = []

        for prediction_id, prediction_counts in counts.items():
            for cell_id, cell_count in prediction_counts.items():
                tabledata = df[prediction_id]
                chromosome_1 = tabledata['chromosome_1']
                strand_1 = tabledata['strand_1']
                position_1 = tabledata['position_1']
                chromosome_2 = tabledata['chromosome_1']
                strand_2 = tabledata['strand_1']
                position_2 = tabledata['position_1']

                output.append(
                    [
                        prediction_id, chromosome_1, strand_1, position_1,
                        chromosome_2, strand_2, position_2, cell_id, cell_count
                    ]
                )

        df = pd.DataFrame(
            output,
            columns=[
                'prediction_id', 'chromosome_1', 'strand_1', 'position_1', 'chromosome_2',
                'strand_2', 'position_2', 'cell_id', 'read_count'
            ]
        )
        return df

    def main(self):
        reads_data = self.load_reads(self.destruct_reads)

        outdata = self._get_count_data(reads_data)

        table = self.load_table()

        table = self.annotate_table(table, outdata)

        csverve.write_dataframe_to_csv_and_yaml(
            table, self.output, self.dtypes
        )
