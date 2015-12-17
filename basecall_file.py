#!/usr/bin/env python2

# the hmm has to be precomputed, otherwise really bad performance.

if __name__ == "__main__":
    argp = argparse.ArgumentParser("Read fast5-files from folder perform the basecalling according to a model. ")
    argp.add_argument("-f", "--fast5dir", required=True, type=str,
                      help="directory with fast 5 input-files.")
    argp.add_argument("-o", "--output", required=True, type=argparse.FileType('w'),
                      help="fasta file with the called reads")
    argp.add_argument("-n", "--ncores", required=False, type=int,
                      help="#CPU cores")

    args = argp.parse_args()