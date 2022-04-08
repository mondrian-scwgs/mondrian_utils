import argparse

from mondrianutils.dlp_utils.dlp_bams_to_mondrian_bam import dlp_bams_to_mondrian_bam


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()

    dlp_bams_to_mondrian_bam = subparsers.add_parser('dlp_bams_to_mondrian_bam')
    dlp_bams_to_mondrian_bam.set_defaults(which='dlp_bams_to_mondrian_bam')
    dlp_bams_to_mondrian_bam.add_argument('--dlp_bam_dir', required=True)
    dlp_bams_to_mondrian_bam.add_argument('--output', required=True)
    dlp_bams_to_mondrian_bam.add_argument('--tempdir', required=True)
    dlp_bams_to_mondrian_bam.add_argument('--cores', default=8)

    args = vars(parser.parse_args())

    return args


def utils():
    args = parse_args()

    if args['which'] == 'dlp_bams_to_mondrian_bam':
        dlp_bams_to_mondrian_bam(args['dlp_bam_dir'], args['output'], args['tempdir'], ncores=args['cores'])
    else:
        raise Exception()
