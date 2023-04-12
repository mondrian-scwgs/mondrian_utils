import pysam


def get_vcf_header(reference_fasta, sample_id):
    fasta = pysam.FastaFile(filename=reference_fasta)

    header = ['##fileformat=VCFv4.2\n']
    for reference in fasta.references:
        header.append(f'##contig=<ID={reference},length={fasta.get_reference_length(reference)}>\n')

    header.append('##ALT=<ID=INS,Description="Insertion">\n')
    header.append('##ALT=<ID=DEL,Description="Deletion">\n')
    header.append('##ALT=<ID=DUP,Description="Duplication">\n')
    header.append('##ALT=<ID=INV,Description="Inversion">\n')
    header.append('##ALT=<ID=BND,Description="Breakend; Translocation">\n')
    header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    header.append('##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">\n')
    header.append(
        '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">\n')
    header.append(
        '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Structural variation with imprecise breakpoints">\n'
    )
    header.append('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">\n')
    header.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">\n')
    header.append('##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromosome for BND SVs">\n')
    header.append(
        '##INFO=<ID=STRAND,Number=1,Type=String,Description="Strands of supporting reads for structural variant">\n'
    )
    header.append(f'#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  {sample_id}\n')

    return header


def read_destruct_calls(destruct_csv):
    with open(destruct_csv, 'rt') as reader:
        header = reader.readline().strip().split('\t')
        header = {v: i for i, v in enumerate(header)}

        for line in reader:
            line = line.strip().split('\t')

            data = dict(
                chromosome_1=line[header['chromosome_1']],
                position_1=int(line[header['position_1']]),
                strand_1=line[header['strand_1']],
                chromosome_2=line[header['chromosome_2']],
                position_2=int(line[header['position_2']]),
                strand_2=line[header['strand_2']],
                prediction_id=line[header['prediction_id']],
                breakpoint_type=line[header['type']],
                rearrangement_type=line[header['rearrangement_type']],
                read_count=line[header['num_reads']],
            )

            yield data


def reclassify_breakpoints(breakpoint_type):
    mapping = {
        'deletion': 'DEL',
        'duplication': 'DUP',
        'inversion': 'BND',
        'foldback': 'BND',
        'translocation': 'BND'
    }

    return mapping[breakpoint_type]


def get_alt(strand_1, strand_2, chrom, pos, breakpoint_type):
    if not breakpoint_type == 'BND':
        return f'<{breakpoint_type}>'

    strands = f'{strand_1}{strand_2}'

    if strands == '+-':
        alt = f'N[{chrom}:{pos}['
    elif strands == '-+':
        alt = f']{chrom}:{pos}]N'
    elif strands == '++':
        alt = f'N]{chrom}:{pos}]'
    elif strands == '--':
        alt = f'[{chrom}:{pos}[N'
    else:
        raise Exception()

    return alt


def get_svlen(position_1, position_2, breakpoint_type):
    svlen = abs(position_1 - position_2)

    if breakpoint_type == 'DEL':
        return -svlen
    elif breakpoint_type == 'DUP':
        return svlen
    else:
        return


def process_single_line_data(call):
    assert call['chromosome_1'] == call['chromosome_2']

    svtype = reclassify_breakpoints(call['breakpoint_type'])

    if call['position_1'] > call['position_2']:
        position_1 = call['position_2']
        strand_1 = call['strand_2']
        position_2 = call['position_1']
        strand_2 = call['strand_1']
    else:
        position_1 = call['position_1']
        strand_1 = call['strand_1']
        position_2 = call['position_2']
        strand_2 = call['strand_2']

    svlen = get_svlen(position_1, position_2, svtype)

    alt = get_alt(strand_1, strand_2, call['chromosome_1'], position_1, svtype)

    yield {
        'chrom': call['chromosome_1'],
        'pos': position_1,
        'end': position_2,
        'id': call['prediction_id'],
        'alt': alt,
        'svtype': svtype,
        'breakpoint_type': call['breakpoint_type'],
        'rearrangement_type': call['rearrangement_type'],
        'svlen': svlen,
        'support': call['read_count'],
        'destruct_id': call['prediction_id'],
        'strand': f'{strand_1}{strand_2}'
    }


def process_double_line_data(call):
    svtype = reclassify_breakpoints(call['breakpoint_type'])

    alt_1 = get_alt(
        call['strand_1'],
        call['strand_2'],
        call['chromosome_2'],
        call['position_2'],
        svtype
    )

    call_1 = {
        'chrom': call['chromosome_1'],
        'pos': call['position_1'],
        'chrom2': call['chromosome_2'],
        'end': call['position_2'],
        'id': call['prediction_id'] + '_A',
        'alt': alt_1,
        'svtype': svtype,
        'breakpoint_type': call['breakpoint_type'],
        'rearrangement_type': call['rearrangement_type'],
        'support': call['read_count'],
        'destruct_id': call['prediction_id'],
        'strand': f'{call["strand_1"]}{call["strand_2"]}'

    }

    yield call_1

    alt_2 = get_alt(
        call['strand_2'],
        call['strand_1'],
        call['chromosome_1'],
        call['position_1'],
        svtype
    )

    call_2 = {
        'chrom': call['chromosome_2'],
        'pos': call['position_2'],
        'chrom2': call['chromosome_1'],
        'end': call['position_1'],
        'id': call['prediction_id'] + '_B',
        'alt': alt_2,
        'svtype': svtype,
        'breakpoint_type': call['breakpoint_type'],
        'rearrangement_type': call['rearrangement_type'],
        'support': call['read_count'],
        'destruct_id': call['prediction_id'],
        'strand': f'{call["strand_2"]}{call["strand_1"]}'
    }

    yield call_2


def process_destruct_vcf_data(destruct_data):
    for call in destruct_data:
        breakpoint_type = reclassify_breakpoints(call['breakpoint_type'])

        if breakpoint_type == 'DUP' or breakpoint_type == 'DEL':
            for v in process_single_line_data(call):
                yield v
        else:
            for v in process_double_line_data(call):
                yield v


def write_vcf(calls, reference, output, sample_id):
    header = get_vcf_header(reference, sample_id)

    with open(output, 'wt') as writer:
        for line in header:
            writer.write(line)

        for call in calls:

            info = [
                'IMPRECISE',
                f'SVTYPE={call["svtype"]}',
                'CIPOS=-100,100', 'CIEND=-100,100',
                f'END={call["end"]}',
                f'SUPPORT={call["support"]}',
                f'BREAKPOINT_TYPE={call["breakpoint_type"]}',
                f'rearrangement_TYPE={call["rearrangement_type"]}',
                f'DESTRUCT_ID={call["destruct_id"]}',
                f'STRAND={call["strand"]}'
            ]
            if 'chrom2' in call:
                info.append(f'CHR2={call["chrom2"]}')
                info.append(f'MATEID={call["id"]}')
            else:
                info.append(f'SVLEN={call["svlen"]}')

            info = ';'.join(info)

            outstr = [
                call['chrom'],
                call['pos'],
                call['id'],
                'N',
                call['alt'],
                60,
                'PASS',
                info,
                'GT:DV',
                './.:2'
            ]

            outstr = '\t'.join([str(v) for v in outstr]) + '\n'

            writer.write(outstr)


def destruct_csv_to_vcf(infile, outfile, reference, sample_id):
    destruct_calls = read_destruct_calls(infile)
    destruct_calls = process_destruct_vcf_data(destruct_calls)

    write_vcf(destruct_calls, reference, outfile, sample_id)
