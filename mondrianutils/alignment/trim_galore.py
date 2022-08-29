'''
Created on Jan 4, 2017

@author: dgrewal
'''

import logging
import os
import shutil
from shutil import copyfile

from mondrianutils import helpers


class RunTrimGalore(object):
    """
    Trim fastq files with trimgalore
    """

    def __init__(self, seq1, seq2, fq_r1, fq_r2, trimgalore, cutadapt, tempdir,
                 adapter, adapter2, report_r1, report_r2, qc_report_r1,
                 qc_report_r2, qc_zip_r1, qc_zip_r2, num_threads, run_fastqc=False):
        self.seq1 = seq1
        self.seq2 = seq2
        self.trimgalore_path = trimgalore
        self.cutadapt_path = cutadapt
        self.tempdir = tempdir
        self.adapter = adapter
        self.adapter2 = adapter2
        self.fastqc_report_r1 = qc_report_r1
        self.fastqc_report_r2 = qc_report_r2
        self.fastqc_zip_r1 = qc_zip_r1
        self.fastqc_zip_r2 = qc_zip_r2
        self.fastq_r1 = fq_r1
        self.fastq_r2 = fq_r2
        self.report_r1 = report_r1
        self.report_r2 = report_r2
        self.empty = False
        self.num_threads = num_threads
        self.run_fastqc = run_fastqc

        self.check_inputs()

        if not os.path.exists(self.tempdir):
            os.makedirs(self.tempdir)

    def check_file(self, path):
        """
        check
        * if they exist or not,
        * are they too small
        * are they gzipped
        """
        if not os.path.exists(path):
            raise IOError("Couldn't find file: %s" % path)

        if os.path.getsize(path) < 1000:
            logging.getLogger("single_cell.align.trim").warn(
                "extremely small file detected: %s" % path)

            with helpers.getFileHandle(path) as reader:
                num_lines = len(reader.readlines())
            if not self.empty and num_lines == 0:
                self.empty = True
                logging.getLogger("single_cell.align.trim").warn(
                    "empty file: %s, skipping trimming" % path)

        ext = os.path.splitext(path)[1][1:]
        # get the last extension
        ext = ext[-1]

        return ext

    def check_inputs(self):
        """
        check input files
        * if they exist or not,
        * are they too small
        * are they gzipped
        """

        ext1 = self.check_file(self.seq1)
        ext2 = self.check_file(self.seq2)

        assert ext1 == ext2, "input fastqfiles should have the same extension"

    def fake_run(self):
        # copy files to output and rename to match trimgalore output
        output_r1 = os.path.join(self.tempdir, "R1_001_val_1.fq.gz")
        output_r2 = os.path.join(self.tempdir, "R2_001_val_2.fq.gz")

        copyfile(self.seq1, output_r1)
        copyfile(self.seq2, output_r2)

    def run_trimgalore(self):
        """
        launch trimgalore
        """

        if self.empty:
            self.fake_run()
            return

        cmd = [self.trimgalore_path]

        if self.run_fastqc:
            cmd.append('--fastqc')

        cmd = ['--paired',
               '--path_to_cutadapt', self.cutadapt_path,
               '--output_dir', self.tempdir + '/',
               ]

        if self.adapter:
            cmd.extend(['--adapter', self.adapter])

        if self.adapter2:
            cmd.extend(['--adapter2', self.adapter2])

        if not self.num_threads == 1:
            cmd.extend(['--fastqc_args', '--threads {}'.format(self.num_threads)])

        cmd.extend([self.seq1, self.seq2])

        helpers.run_cmd(cmd)

    def move_files(self, fname, outpath):
        """
        move files from from temp dir to the expected path
        """
        dir = os.path.dirname(outpath)


        if dir and not os.path.exists(dir):
            os.makedirs(dir)

        path = os.path.join(self.tempdir, fname)
        shutil.move(path, outpath)
        assert os.path.isfile(outpath)

    def get_file(self, r1_out, r2_out, ext, trimreport=False):
        """
        find the file in outdir and rename to r1_out and r2_out
        """
        outfiles = os.listdir(self.tempdir)

        reps = [v for v in outfiles if ext in v]

        if ext == '.fq.gz':
            assert reps != [], "Couldn't move %s files" % ext

        seq1 = os.path.basename(self.seq1)
        seq2 = os.path.basename(self.seq2)

        if trimreport:
            for rep in reps:
                if seq1 in rep:
                    self.move_files(rep, r1_out)
                elif seq2 in rep:
                    self.move_files(rep, r2_out)
                else:
                    if ext == '.fq.gz':
                        raise Exception("Couldn't move %s files" % ext)
                    else:
                        logging.getLogger("single_cell.align.trim").warn(
                            "Couldn't move %s files" % ext)
            return

        for rep in reps:
            if "val_1" in rep:
                self.move_files(rep, r1_out)
            elif "val_2" in rep:
                self.move_files(rep, r2_out)
            else:
                if ext == '.fq.gz':
                    raise Exception("Couldn't move %s files" % ext)
                else:
                    logging.getLogger("single_cell.align.trim").warn(
                        "Couldn't move %s files" % ext)

    def gather_outputs(self):
        """
        rename the output metrics files.
        """
        # trimming reports
        self.get_file(
            self.report_r1,
            self.report_r2,
            "_trimming_report.txt",
            trimreport=True)

        # fastqc html files
        self.get_file(
            self.fastqc_report_r1,
            self.fastqc_report_r2,
            ".html")

        # fastqc zip files
        self.get_file(self.fastqc_zip_r1, self.fastqc_zip_r2, ".zip")

        # trimmed fastq files
        self.get_file(self.fastq_r1, self.fastq_r2, ".fq.gz")


def trim_galore(fastq1, fastq2, trim_fastq1, trim_fastq2, adapter1, adapter2, tempdir, num_threads):
    """
    run fastqc on both fastq files
    run trimgalore if needed, copy if not.
    """
    with helpers.getFileHandle(fastq1) as reader:
        if not reader.readline():
            shutil.copyfile(fastq1, trim_fastq1)
            shutil.copyfile(fastq2, trim_fastq2)

    reports_dir = os.path.join(tempdir, 'fastqc_reports')
    helpers.makedirs(reports_dir)
    rep1 = os.path.join(reports_dir, 'trimgalore_R1.html')
    rep2 = os.path.join(reports_dir, 'trimgalore_R2.html')
    qcrep1 = os.path.join(reports_dir, 'trimgalore_qc_R1.html')
    qcrep2 = os.path.join(reports_dir, 'trimgalore_qc_R2.html')
    qczip1 = os.path.join(reports_dir, 'trimgalore_qc_R1.zip')
    qczip2 = os.path.join(reports_dir, 'trimgalore_qc_R2.zip')

    run_tg = RunTrimGalore(
        fastq1, fastq2, trim_fastq1, trim_fastq2, 'trim_galore', 'cutadapt', tempdir,
        adapter1, adapter2, rep1, rep2, qcrep1, qcrep2, qczip1, qczip2, num_threads
    )
    run_tg.run_trimgalore()
    run_tg.gather_outputs()
