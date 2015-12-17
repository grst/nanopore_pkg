#!/usr/bin/env python3

import argparse
from npcaller import fast5
import h5py
import re
import pandas
import numpy as np
import pickle


def get_models(file):
    """
    Extract the Metrichor-models for the HMM-basecaller from a fast5 file.

    Args:
        file: fast5-file

    Returns:
        Dict of Name => Dataframe, each containing a model
    """
    f5 = h5py.File(file, 'r')

    models = []

    # collect model names
    try:
        log = f5['/Analyses/Basecall_2D_000/Log']
        log = bytes(log[...]).decode("utf-8").split("\n")
        model = [x for x in log if x.find(".model") >= 0]
        model = [re.search(r"\"(.*)\"", x).group(1) for x in model]
        models.extend(model)
    except KeyError:
        raise fast5.Fast5PathException("Error opening fast5 file or unsupported file structure.")

    # collect models
    model_params = {}
    for t in ["template", "complement"]:
        try:
            m_name = [x for x in model if x.find(t) >= 0][0]
            m_name = m_name.split("/")[-1].split(".")[0] # get model name from file path
            f5_path = fast5.identify_events_path(f5)
            params = pandas.DataFrame(np.array(f5[(f5_path + '/Model').format(t)]))
            params["kmer"] = params["kmer"].map(lambda x: x.decode("utf-8"))
            model_params[m_name] = params
            print(">> Model found: {0}".format(m_name))
            print("Mean {0}\t Median {1}".format(params["level_mean"].median(), params["level_mean"].mean()))
        except IndexError:
            """model not available in this file"""
            print("{0}-model is specified in log, but not found in the file. Skipped.".format(t))

    return model_params


def dump_models(model_params, basename):
    """
    Dump models to one file each.

    Args:
        model_params: output of get_models()
        basename: will be appended with model name.

    """

    for name, model in model_params.items():
        pickle.dump(model, open(basename + "." + name + ".pickle", 'wb'))


if __name__ == "__main__":
    argp = argparse.ArgumentParser("Extracts all models contained in a given fast5 file")
    argp.add_argument("-f", "--fast5", required=True, type=str, help="Path to fast5 file")
    argp.add_argument("-o", "--output", required=True, type=str,
                      help="model basename. /path/to/my_model will result in e.g. /path/to/my_model.template.pickle")

    args = argp.parse_args()

    model_params = get_models(args.fast5)
    dump_models(model_params, args.output)