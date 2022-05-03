from itertools import product


def _get_bam_flag(
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


def _get_bam_header():
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


def _get_read(flag, mapq=60, read_idx=0):
    read1 = [
        'HISEQ101_{}:5:1206:10358:49496'.format(read_idx),
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
        'HISEQ101_{}:5:1206:10358:49496'.format(read_idx),
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

    return read1 + read2


def get_all_fail_combinations_sam():
    flag_combos = [v for v in product([True, False], repeat=5)]
    flag_combos = [v for v in flag_combos if any(v)]
    flags = [_get_bam_flag(*v) for v in flag_combos]

    reads = [_get_read(flag, read_idx=i) for i, flag in enumerate(flags)]
    reads = ''.join(reads)

    header = _get_bam_header()

    sam = header + reads

    return sam


def get_duplicate_sam():
    flag = _get_bam_flag(duplicate=True)
    return _get_bam_header() + _get_read(flag, read_idx=0)

def get_unpaired_sam():
    flag = _get_bam_flag(unpaired=True)
    return _get_bam_header() + _get_read(flag, read_idx=0)

def get_unmapped_sam():
    flag = _get_bam_flag(unmapped=True)
    return _get_bam_header() + _get_read(flag, read_idx=0)

def get_secondary_sam():
    flag = _get_bam_flag(secondary=True)
    return _get_bam_header() + _get_read(flag, read_idx=0)

def get_supplementary_sam():
    flag = _get_bam_flag(supplementary=True)
    return _get_bam_header() + _get_read(flag, read_idx=0)
