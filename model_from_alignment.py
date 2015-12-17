#!/usr/bin/env python3

import argparse



if __name__ == "__main__":
    argp = argparse.ArgumentParser("Align processed reads from metrichor to reference; "
                                   "extract correct kmers and calculate mean and stdv for the model. "
                                   "Only data from 2D reads is used."
                                   ""
                                   "Warning. This script loads all fast5-files to RAM. Ensure that you have enough "
                                   "or do not use all the data. ")
    argp.add_argument("-f", "--fast5dir", required=True, type=str,
                      help="directory with fast 5 files.")
    argp.add_argument("-r", "--reference", required=True, type=argparse.FileType('r'),
                      help="fasta file with reference sequence")
    argp.add_argument("-o", "--output", required=True, type=str,
                      help="model basename. /path/to/my_model will result in e.g. /path/to/my_model.template.pickle")
    argp.add_argument("-n", "--ncores", required=False, type=int,
                      help="#CPU cores")
    argp.add_argument("-k", "--kmer", required=False, type=int,
                      help="length of kmer", default=6)


    args = argp.parse_args()