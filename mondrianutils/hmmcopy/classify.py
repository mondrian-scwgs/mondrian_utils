'''
Created on Jun 26, 2018

@author: dgrewal
'''
import joblib
import csverve.api as csverve
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from mondrianutils import helpers
from mondrianutils.dtypes.hmmcopy import dtypes as hmmcopy_dtypes


def read_from_h5(filename, tablename):
    with pd.HDFStore(filename) as h5store:
        data = h5store[tablename]
    return data


def read_data(filename, tablename):
    if filename.endswith('.h5'):
        data = read_from_h5(filename, tablename)
    elif filename.endswith('.csv'):
        data = csverve.read_csv(filename)
    elif filename.endswith('.csv.gz'):
        data = csverve.read_csv(filename)
    else:
        raise Exception("unknown file format")

    return data


def train_classifier(filename, joblib_model=None):
    if joblib_model:
        try:
            return joblib.load(joblib_model)
        except:
            print("unable to load joblib model, loading from h5")

    training_data = read_from_h5(filename, '/training_data')

    labels = training_data["label"]

    del training_data["label"]

    clf = RandomForestClassifier(n_estimators=500)

    model = clf.fit(training_data, labels)

    features = training_data.columns.values
    model.feature_names_ = features

    return model


def load_data(metrics, colnames):
    metrics = csverve.read_csv(metrics)

    metrics = metrics.set_index('cell_id')

    data = []
    for colname in colnames:
        coldata = metrics[colname]

        if colname == 'scaled_halfiness':
            # haploid poison adds inf, replace with big number since 0 is considered good
            # and we want to score to decrease
            coldata = coldata.replace(np.inf, 1e10)
        data.append(coldata)

    data = pd.concat(data, axis=1)
    data = data.replace(-np.inf, np.nan)
    data = data.fillna(0)

    return data


def classify(model, data):
    predictions = model.predict_proba(data)

    index_good_proba = np.where(model.classes_ == 1)

    assert len(index_good_proba) == 1

    index_good_proba = index_good_proba[0]

    predictions = predictions[:, index_good_proba]

    predictions = dict(zip(data.index, predictions))

    return predictions


def write_to_hdf(output, hmmcopy_tablename, data):
    with pd.HDFStore(output, 'a', complevel=9, complib='blosc') as outstore:
        outstore.put(hmmcopy_tablename, data, format='table')


def write_to_csv(output, data, gzipped=False):
    compression = 'gzip' if gzipped else None
    data.to_csv(output, index=False, compression=compression, na_rep='NA')


def write_to_output(hmmcopy_filename, output, predictions):
    data = csverve.read_csv(hmmcopy_filename)

    data['cell_id'] = data['cell_id'].astype('str')

    data['quality'] = data['cell_id'].map(predictions)
    data.quality = data.quality.astype(float)

    organisms = [v for v in data.columns.values if v.startswith('fastqscreen_')]
    organisms = sorted(set([v.split('_')[1] for v in organisms]))
    organisms = [v for v in organisms if v not in ['nohit', 'total']]

    csverve.write_dataframe_to_csv_and_yaml(
        data, output,
        dtypes=hmmcopy_dtypes(fastqscreen_genomes=organisms)['metrics'],
        write_header=True
    )


def add_quality(metrics, output, training_data_h5, joblib_model=None):
    model = train_classifier(training_data_h5, joblib_model=joblib_model)

    feature_names = model.feature_names_

    data = load_data(metrics, feature_names)

    predictions = classify(model, data)

    write_to_output(
        metrics,
        output,
        predictions)
