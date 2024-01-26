'''
Created on Sep 8, 2015

@author: dgrewal
'''
import sys

import matplotlib

matplotlib.use("Agg")
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import logging
from mondrianutils import helpers
from .clustermap import ClusterMap
import csverve.api as csverve

sys.setrecursionlimit(2000)


class PlotPcolor(object):
    '''
    merges files. no overlap queries, simple concatenation
    since columns are different, select header and insert values at proper
    indices. use N/A for missing.
    '''

    def __init__(
            self,
            infile,
            metrics,
            output,
            column_name='state',
            max_cn=20,
            chromosomes=[str(v) for v in range(1, 23)] + ['X', 'Y'],
            mappability_threshold=0.9,
            sidebar_column='pick_met',
            scale_by_cells=False,
            disable_clustering=False
    ):
        self.input = infile
        self.metrics = metrics
        self.output = output

        self.column_name = column_name
        self.max_cn = max_cn

        self.chromosomes = chromosomes
        self.mappability_threshold = mappability_threshold

        self.scale_by_cells = scale_by_cells
        self.sidebar_column = sidebar_column
        self.disable_clustering = disable_clustering

    def read_segs(self):
        return helpers.load_and_pivot_reads_data(self.input, self.column_name)

    def read_metrics(self):
        metrics = csverve.read_csv(self.metrics)
        metrics = metrics.set_index('cell_id')
        return metrics

    def plot_heatmap(self, data, ccdata, title, lims, pdfout, distance_matrix=None):
        """
        generate heatmap, annotate and save

        """
        ClusterMap(
            data,
            ccdata,
            self.max_cn,
            chromosomes=self.chromosomes,
            scale_by_cells=self.scale_by_cells,
            distance_matrix=distance_matrix
        )

        plt.suptitle(title)

        plt.subplots_adjust(right=0.85)

        pdfout.savefig(pad_inches=0.2)

        plt.close("all")

    def plot_heatmap_by_sep(self, data, sepdata, colordata):
        """
        generate and save plot to output
        """

        def genplot(data, samples):
            pltdata = data.loc[samples]

            title = ' (%s) n=%s/%s' % (sep, len(samples), num_samples)

            if self.disable_clustering:
                allvals = {v: i for i, v in enumerate(set(colordata.values()))}
                dist_mat = [allvals[colordata[v]] for v in pltdata.index]
                dist_mat = np.matrix([[v] * 10 for v in dist_mat])
            else:
                dist_mat = None

            self.plot_heatmap(
                pltdata, colordata, title, lims, pdfout, distance_matrix=dist_mat
            )

        if not self.output:
            return

        sns.set_style('whitegrid')
        sns.set(font_scale=1.5)

        pdfout = PdfPages(self.output)

        if not data.values.size:
            logging.getLogger("single_cell.plot_heatmap").warn("no data to plot")
            return

        vmax = np.nanmax(data.values)
        vmin = np.nanmin(data.values)
        lims = (vmin, vmax)

        for sep, samples in sepdata.items():
            num_samples = len(samples)

            samples = set(samples).intersection(set(data.index))

            if len(samples) < 2:
                continue

            if len(samples) > 3000:
                logging.getLogger("single_cell.plot_heatmap").warn(
                    'The output file will only plot 1000 cells per page'
                )
                samples = sorted(samples)
                # plot in groups of 1000
                sample_sets = [samples[x:x + 1000]
                               for x in range(0, len(samples), 1000)]
                if len(sample_sets[-1]) == 1:
                    sample_sets[-2] += sample_sets[-1]
                    del sample_sets[-1]

                for samples in sample_sets:
                    genplot(data, list(samples))
            else:
                genplot(data, list(samples))
        pdfout.close()

    def main(self):
        '''
        main function
        '''
        data = self.read_segs()

        if data.empty:
            logging.getLogger("single_cell.plot_heatmap").warn("no data to plot")
            open(self.output, "w").close()
            return

        metrics = self.read_metrics()

        sepdata = {'all': list(metrics.index)}
        colordata = metrics[self.sidebar_column].to_dict()

        self.plot_heatmap_by_sep(data, sepdata, colordata)
