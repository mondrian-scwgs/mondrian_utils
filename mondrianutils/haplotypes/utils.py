import csverve.api as csverve
import json
import os
import pandas as pd
import yaml
from mondrianutils import helpers
from mondrianutils.dtypes.haplotypes import dtypes
from mondrianutils import __version__


def add_cell_id_to_seqdata(seqdata, cellid):
    with pd.HDFStore(seqdata) as store:
        store.put('cell_id', pd.Series([cellid]))


def get_cell_id_from_seqdata(seqdata):
    with pd.HDFStore(seqdata) as store:
        if '/cell_id' not in store.keys():
            return
        cellid = store.get('cell_id')
        cellid = list(cellid)
        assert len(cellid) == 1
        return cellid[0]


def infer_type(files):
    with open(files, 'rt') as reader:
        files = json.load(reader)

    filetypes = sorted(set([v['left'] for v in files]))

    # more than one wf
    if 'haplotype_counts' in filetypes and 'infer_haplotype' in filetypes:
        return 'haplotype_calling'
    elif 'haplotype_counts' in filetypes:
        return 'haplotype_counting'
    elif 'infer_haplotype' in filetypes:
        return 'infer_haplotype'
    else:
        raise Exception()


def generate_metadata(
        files, metadata_yaml_files, samples, metadata_output
):
    wf_type = infer_type(files)
    data = helpers.metadata_helper(files, metadata_yaml_files, samples, wf_type)

    with open(metadata_output, 'wt') as writer:
        yaml.dump(data, writer, default_flow_style=False)


def generate_infer_haps_metadata(
        csvfile, yamlfile, metadata_input, sample, metadata_output
):
    with open(metadata_input, 'rt') as reader:
        data = yaml.safe_load(reader)

    out_data = dict()
    out_data['meta'] = dict(
        type='haplotype_infer',
        version=__version__,
        lane_ids=data['meta']['lane_ids'],
        sample_ids=data['meta']['sample_ids'],
        library_ids=data['meta']['library_ids'],
        cell_ids=data['meta']['cell_ids'],
        sample=sample
    )

    files = {
        os.path.basename(csvfile): {
            'result_type': 'infer_haplotype',
            'auxiliary': helpers.get_auxiliary_files(csvfile)
        },
        os.path.basename(yamlfile): {
            'result_type': 'infer_haplotype',
            'auxiliary': helpers.get_auxiliary_files(yamlfile)
        }
    }

    with open(metadata_output, 'wt') as writer:
        yaml.dump(files, writer, default_flow_style=False)


def generate_count_haps_metadata(
        csvfile, yamlfile, barcodes, variants, ref_counts, alt_counts,
        metadata_input, sample, metadata_output
):
    with open(metadata_input, 'rt') as reader:
        data = yaml.safe_load(reader)

    out_data = dict()
    out_data['meta'] = dict(
        type='haplotype_count',
        version=__version__,
        lane_ids=data['meta']['lane_ids'],
        sample_ids=data['meta']['sample_ids'],
        library_ids=data['meta']['library_ids'],
        cell_ids=data['meta']['cell_ids'],
        sample=sample
    )

    files = {
        os.path.basename(csvfile): {
            'result_type': 'infer_haplotype',
            'auxiliary': helpers.get_auxiliary_files(csvfile)
        },
        os.path.basename(yamlfile): {
            'result_type': 'infer_haplotype',
            'auxiliary': helpers.get_auxiliary_files(yamlfile)
        },
        os.path.basename(csvfile): {
            'result_type': 'infer_haplotype',
            'auxiliary': helpers.get_auxiliary_files(barcodes)
        },
        os.path.basename(yamlfile): {
            'result_type': 'infer_haplotype',
            'auxiliary': helpers.get_auxiliary_files(variants)
        },
        os.path.basename(csvfile): {
            'result_type': 'infer_haplotype',
            'auxiliary': helpers.get_auxiliary_files(ref_counts)
        },
        os.path.basename(yamlfile): {
            'result_type': 'infer_haplotype',
            'auxiliary': helpers.get_auxiliary_files(alt_counts)
        }
    }

    with open(metadata_output, 'wt') as writer:
        yaml.dump(files, writer, default_flow_style=False)


def finalize_tsv(infile, outfile, seqdata, skip_header=False):
    df = pd.read_csv(infile, sep='\t')

    cellid = get_cell_id_from_seqdata(seqdata)
    if cellid:
        df['cell_id'] = cellid

    csverve.write_dataframe_to_csv_and_yaml(
        df, outfile, skip_header=skip_header, dtypes=dtypes()
    )


def annotate_haps(haps_csv, thousand_genomes, tempdir, output_csv):
    helpers.makedirs(tempdir)
    temp_output = os.path.join(tempdir, 'output.csv')

    annotation_data = {}

    with helpers.getFileHandle(thousand_genomes, 'rt') as db:
        for line in db:
            line = line.strip().split('\t')

            chrom, pos, ref, alt = line

            annotation_data[(chrom, pos)] = (ref, alt)

    with helpers.getFileHandle(haps_csv, 'rt') as reader, helpers.getFileHandle(temp_output, 'wt') as writer:

        header = reader.readline().strip().split('\t')
        header.extend(['ref', 'alt'])
        header = ','.join(header) + '\n'
        writer.write(header)

        for line in reader:
            line = line.strip().split('\t')

            chrom = line[0]
            pos = line[1]

            if (chrom, pos) in annotation_data:
                ref, alt = annotation_data[(chrom, pos)]
            else:
                ref = 'NA'
                alt = 'NA'

            line.extend([ref, alt])
            line = ','.join(line) + '\n'

            writer.write(line)

    csverve.rewrite_csv_file(
        temp_output, output_csv,
        dtypes=dtypes()
    )


def convert_csv_to_tsv(csv_infile, tsv_outfile):
    df = csverve.read_csv(csv_infile)
    df.to_csv(tsv_outfile, sep='\t', index=False)
