import os

import mondrianutils.alignment.utils as alignment_utils
import mondrianutils.helpers as helpers


def get_bam_files(dirpath):
    files = os.listdir(dirpath)
    files = [filepath for filepath in files if filepath.endswith('.bam')]
    files = [filepath for filepath in files if 'MT' not in filepath]
    files = [os.path.join(dirpath, filepath) for filepath in files]
    return files


def get_cell_id(bampath):
    cell_id = os.path.basename(bampath)
    cell_id = cell_id.replace('.bam', '')
    return cell_id


def tag_bam_cmd(input_bam, output_bam, cell_id):
    cmd = ['alignment_utils', 'tag_bam_with_cellid',
           '--infile', input_bam, '--outfile', output_bam,
           '--cell_id', cell_id]

    return cmd


def dlp_bams_to_mondrian_bam(dlp_bam_dir, output_bam, tempdir, ncores=8):
    dlp_bam_files = get_bam_files(dlp_bam_dir)

    tagged_bam_temp = os.path.join(tempdir, 'tagged_bams')
    helpers.makedirs(tagged_bam_temp)
    tag_cmds = []
    all_tagged_bams = {}
    for bamfile in dlp_bam_files:
        tagged_bam = os.path.join(tagged_bam_temp, os.path.basename(bamfile))
        cell_id = get_cell_id(bamfile)
        all_tagged_bams[cell_id] = tagged_bam

        cmd = tag_bam_cmd(bamfile, tagged_bam, cell_id)

        tag_cmds.append(cmd)

    scripts_tempdir = os.path.join(tempdir, 'scripts')
    helpers.makedirs(scripts_tempdir)
    helpers.run_in_gnu_parallel(tag_cmds, scripts_tempdir, ncores)

    merge_cells_temp = os.path.join(tempdir, 'merge_cells')
    helpers.makedirs(merge_cells_temp)
    alignment_utils.merge_cells(all_tagged_bams, merge_cells_temp, ncores, output_bam)
