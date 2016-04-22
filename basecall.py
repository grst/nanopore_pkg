#!/usr/bin/env python2
import sys
import argparse
from npcaller.basecaller import Basecaller


def main(filelist, output, template, complement, ncores):
    bc = Basecaller(ncores=ncores, template_model=template, complement_model=complement)
    files = [l.strip() for l in filelist.readlines()]
    bc.process_files(files)
    bc.write_to_fasta(output)


if __name__ == "__main__":
    argp = argparse.ArgumentParser("Read fast5-files from folder perform the basecalling according to a model. ")
    argp.add_argument("-f", "--filelist", required=True, type=argparse.FileType('r'),
                      help="a list of fast5 files to be called, one per line")
    argp.add_argument("-o", "--output", required=True, type=argparse.FileType('w'),
                      help="fasta file with the called reads")
    argp.add_argument("-t", "--template", required=False, type=argparse.FileType('rb'),
                      help="path to template model")
    argp.add_argument("-c", "--complement", required=False, type=argparse.FileType('rb'),
                      help="path to complement model")
    argp.add_argument("-n", "--ncores", required=False, type=int, default=None,
                      help="#CPU cores")

    args = argp.parse_args()
    if args.template is None and args.complement is None:
        print("Template and/or Complement model must be specified")
        sys.exit(1)

    main(args.filelist, args.output, args.template, args.complement, args.ncores)