import os

from mondrianutils.breakpoint_calling.destruct_csv_to_vcf import (
    reclassify_breakpoints,
    get_alt,
    process_single_line_data,
    get_svlen,
    process_destruct_vcf_data,
    process_double_line_data,
    get_vcf_header,
    write_vcf
)


def reference_fasta(fasta_path):
    with open(fasta_path + '.fai', 'wt') as writer:
        writer.write('1\t249250621\t75\t60\t61\n')
        writer.write('2\t243199373\t253404949\t60\t61\n')

    open(fasta_path, 'wt').close()


def test_get_vcf_header(tmpdir):
    reference = os.path.join(tmpdir, 'reference.fasta')
    reference_fasta(reference)

    header = get_vcf_header(reference, 'SA123')
    assert isinstance(header, list)
    assert len(header) > 0
    assert header[0] == "##fileformat=VCFv4.2\n"

    header = [v.strip() for v in header]
    assert '##contig=<ID=1,length=249250621>' in header
    assert '##contig=<ID=2,length=243199373>' in header

    assert header[-1].split() == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SA123']


def test_reclassify_breakpoints():
    assert reclassify_breakpoints("deletion") == "DEL"
    assert reclassify_breakpoints("duplication") == "DUP"
    assert reclassify_breakpoints("inversion") == "BND"
    assert reclassify_breakpoints("foldback") == "BND"
    assert reclassify_breakpoints("translocation") == "BND"


def test_get_alt():
    assert get_alt("+", "-", "chr1", 100, "BND") == "N[chr1:100["
    assert get_alt("-", "+", "chr2", 200, "BND") == "]chr2:200]N"
    assert get_alt("+", "+", "chr3", 300, "BND") == "N]chr3:300]"
    assert get_alt("-", "-", "chr4", 400, "BND") == "[chr4:400[N"


def test_get_svlen():
    assert get_svlen(100, 200, "DEL") == -100
    assert get_svlen(200, 100, "DEL") == -100
    assert get_svlen(100, 200, "DUP") == 100
    assert get_svlen(200, 100, "DUP") == 100
    assert get_svlen(100, 200, "BND") is None


def test_process_double_line_data():
    call = {
        'chromosome_1': 'chr1',
        'position_1': 100,
        'strand_1': '+',
        'chromosome_2': 'chr2',
        'position_2': 200,
        'strand_2': '-',
        'prediction_id': 'call1',
        'breakpoint_type': 'translocation',
        'rearrangement_type': 'type1',
        'read_count': 2
    }
    result = list(process_double_line_data(call))
    assert len(result) == 2
    assert isinstance(result[0], dict)
    assert isinstance(result[1], dict)

    assert result[0]['chrom'] == 'chr1'
    assert result[1]['chrom'] == 'chr2'
    assert result[0]['pos'] == 100
    assert result[1]['pos'] == 200
    assert result[0]['chrom2'] == 'chr2'
    assert result[1]['chrom2'] == 'chr1'
    assert result[0]['end'] == 200
    assert result[1]['end'] == 100
    assert result[0]['id'] == 'call1_A'
    assert result[1]['id'] == 'call1_B'
    assert result[0]['alt'] == 'N[chr2:200['
    assert result[1]['alt'] == ']chr1:100]N'
    assert result[0]['svtype'] == result[1]['svtype'] == 'BND'
    assert result[0]['breakpoint_type'] == result[1]['breakpoint_type'] == 'translocation'
    assert result[0]['rearrangement_type'] == result[1]['rearrangement_type'] == 'type1'
    assert result[0]['support'] == result[1]['support'] == 2


def test_process_single_line_data():
    call = {
        'chromosome_1': 'chr1',
        'position_1': 100,
        'strand_1': '+',
        'chromosome_2': 'chr1',
        'position_2': 200,
        'strand_2': '-',
        'prediction_id': 'call1',
        'breakpoint_type': 'deletion',
        'rearrangement_type': 'type1',
        'read_count': 2
    }
    result = list(process_single_line_data(call))
    assert len(result) == 1
    assert isinstance(result[0], dict)


def test_process_destruct_vcf_data():
    destruct_data = [
        {
            'chromosome_1': 'chr1',
            'position_1': 100,
            'strand_1': '+',
            'chromosome_2': 'chr1',
            'position_2': 200,
            'strand_2': '-',
            'prediction_id': 'call1',
            'breakpoint_type': 'deletion',
            'rearrangement_type': 'type1',
            'read_count': 2
        },
        {
            'chromosome_1': 'chr1',
            'position_1': 100,
            'strand_1': '+',
            'chromosome_2': 'chr2',
            'position_2': 200,
            'strand_2': '-',
            'prediction_id': 'call2',
            'breakpoint_type': 'translocation',
            'rearrangement_type': 'type2',
            'read_count': 3
        }
    ]
    result = list(process_destruct_vcf_data(destruct_data))
    assert len(result) == 3
    assert isinstance(result[0], dict)
    assert isinstance(result[1], dict)
    assert isinstance(result[2], dict)


def test_write_vcf_header(tmpdir):
    reference = os.path.join(tmpdir, 'reference.fasta')
    reference_fasta(reference)
    sample_id = "SA123"
    output_vcf = os.path.join(tmpdir, 'output.vcf')

    write_vcf([], reference, output_vcf, sample_id)

    with open(output_vcf, 'rt') as reader:
        data = [v.strip() for v in reader]
        assert '##contig=<ID=1,length=249250621>' in data
        assert '##contig=<ID=2,length=243199373>' in data
        assert data[-1].split() == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SA123']


def test_write_vcf_del(tmpdir):
    reference = os.path.join(tmpdir, 'reference.fasta')
    reference_fasta(reference)
    sample_id = "SA123"
    output_vcf = os.path.join(tmpdir, 'output.vcf')

    calls = [
        {
            'chrom': 'chr1',
            'pos': 100,
            'end': 200,
            'id': 'call1',
            'alt': 'DEL',
            'svtype': 'deletion',
            'breakpoint_type': 'deletion',
            'rearrangement_type': 'type1',
            'svlen': -100,
            'support': 2,
            'destruct_id': 'call1',
            'strand': '+-'
        }
    ]
    write_vcf(calls, reference, output_vcf, sample_id)

    with open(output_vcf, 'rt') as reader:
        data = [v.strip() for v in reader]

        info = 'IMPRECISE;SVTYPE=deletion;CIPOS=-100,100;CIEND=-100,100;END=200;SUPPORT=2;' \
               'BREAKPOINT_TYPE=deletion;rearrangement_TYPE=type1;DESTRUCT_ID=call1;STRAND=+-;SVLEN=-100'
        assert data[-1].split('\t') == ['chr1', '100', 'call1', 'N', 'DEL', '60', 'PASS', info, 'GT:DV', './.:2']


def test_write_vcf_bnd(tmpdir):
    reference = os.path.join(tmpdir, 'reference.fasta')
    reference_fasta(reference)
    sample_id = "SA123"
    output_vcf = os.path.join(tmpdir, 'output.vcf')

    calls = [
        {'chrom': 'chr1', 'pos': 100, 'chrom2': 'chr2',
         'end': 200, 'id': 'call1_A', 'alt': 'N[chr2:200[',
         'svtype': 'BND', 'breakpoint_type': 'translocation',
         'rearrangement_type': 'type1', 'support': 2,
         'destruct_id': 'call1', 'strand': '+-'
         },
        {'chrom': 'chr2', 'pos': 200, 'chrom2': 'chr1',
         'end': 100, 'id': 'call1_B', 'alt': ']chr1:100]N',
         'svtype': 'BND', 'breakpoint_type': 'translocation',
         'rearrangement_type': 'type1', 'support': 2,
         'destruct_id': 'call1', 'strand': '-+'
         }
    ]
    write_vcf(calls, reference, output_vcf, sample_id)

    with open(output_vcf, 'rt') as reader:
        data = [v.strip() for v in reader]

        info = 'IMPRECISE;SVTYPE=BND;CIPOS=-100,100;CIEND=-100,100;END=200;SUPPORT=2;BREAKPOINT_TYPE=translocation;' \
               'rearrangement_TYPE=type1;DESTRUCT_ID=call1;STRAND=+-;CHR2=chr2;MATEID=call1_A'
        assert data[-2].split('\t') == ['chr1', '100', 'call1_A', 'N', 'N[chr2:200[', '60', 'PASS', info, 'GT:DV',
                                        './.:2']

        info = 'IMPRECISE;SVTYPE=BND;CIPOS=-100,100;CIEND=-100,100;END=100;SUPPORT=2;BREAKPOINT_TYPE=translocation;' \
               'rearrangement_TYPE=type1;DESTRUCT_ID=call1;STRAND=-+;CHR2=chr1;MATEID=call1_B'
        assert data[-1].split('\t') == ['chr2', '200', 'call1_B', 'N', ']chr1:100]N', '60', 'PASS', info, 'GT:DV',
                                        './.:2']
