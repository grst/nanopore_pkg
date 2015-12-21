#!/usr/bin/env python3
import argparse

if __name__ == "__main__":
    argp = argparse.ArgumentParser("Map reads to reference and calculate statistics.")
    argp.add_argument("-f", "--fasta", required=True, type=argparse.FileType('r'),
                      help="directory with fast 5 files.")
    argp.add_argument("-r", "--reference", required=True, type=argparse.FileType('r'),
                      help="fasta file with reference sequence")
    argp.add_argument("-g", "--graphmap", required=False, type=str, default="",
                      help="Path to graphmap alignment tool")
    argp.add_argument("-n", "--ncores", required=False, type=int,
                      help="#CPU cores")

    args = argp.parse_args()