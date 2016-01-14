#!/usr/bin/env python3
import argparse
from npcaller.fasta import FastaReader, FastaWriter
from npcaller.fast5 import Fast5File
from npcaller.validator import align_to_reference, sam_to_bam
from tempfile import NamedTemporaryFile
from pysam import AlignmentFile
import logging
from itertools import chain
logger = logging.getLogger("alignment")
logger.setLevel(logging.INFO)


class AlignmentEndException(Exception):
    pass


class ModelMaker(object):
    """
    Make a model from alignment of events to reference.
    """
    def __init__(self, filelist, ref_file, model_basename, ncores, k, graphmap_bin):
        """

        Args:
            filelist:
            ref_file:
            model_basename:
            ncores:
            k:

        Returns:

        """
        self.ref_file = ref_file
        self.model_basename = model_basename
        self.ncores = ncores
        self.k = k
        self.graphmap_bin = graphmap_bin

        fr = FastaReader(ref_file)
        _, self.ref_seq = next(fr.get_entries())

        self.f5files = {}
        filelist = [l.strip() for l in filelist.readlines()]
        for file in filelist:
            f5file = Fast5File(file)
            self.f5files[f5file.get_id(), f5file]

        self.alignment = self.make_bam()

        self.find_correct_kmers()

        self.make_model()

    def make_bam(self):
        with NamedTemporaryFile() as fasta_file:
            samfile = self.model_basename + ".sam"
            bamfile = self.model_basename
            fw = FastaWriter(fasta_file)
            for f5file in self.f5files.values():
                fw.write_entry(f5file.get_id(), Fast5File.get_seq("template"))
            fw.flush()
            align_to_reference(fasta_file, self.ref_file, samfile,
                               graphmap_bin=self.graphmap_bin, ncores=self.ncores)
            sam_to_bam(samfile, bamfile)

        return AlignmentFile(bamfile)

    def find_correct_kmers(self):
        total_events = 0
        result = list()
        reads = [x for x in self.alignment.fetch()]
        logger.info("{0} reads found in alignment".format(len(reads)))
        for read in reads:
            f5file = self.f5files[read.query_name]
            pairs = [list(t) for t in zip(*read.get_aligned_pairs())]
            assert(pairs[0][0] == 0), "alignment is not null-indexed."
            correct, total = self._process_events(f5file, pairs)
            total_events += total
            result.append(correct)

        true_events = list(chain.from_iterable([r.get() for r in result]))
        logger.info("Identified {0} correct kmers of {1} total kmers. That's {2:%}".format(
            true_events, total_events, true_events/total_events
        ))

    def make_model(self):
        pass

    def _process_events(self, f5file, pairs):
        i_seq = 0
        correct = list()
        total = 0
        strand = "template"
        called_seq = f5file.get_seq(strand)
        for ev in f5file.get_corrected_events(strand):
            total += 1
            i_seq += self._gapmove(ev["move"], pairs[0], i_seq)
            try:
                ev_index = self._event_indexes(pairs[0], i_seq)
            except AlignmentEndException:
                # not the whole read is aligned
                break
            read_kmer = self._get_nt_kmer(ev_index, pairs[0], called_seq)
            assert(read_kmer == ev["kmer"]), (i_seq, ev, read_kmer, ev_index)
            if self._is_correct_kmer(ev_index, pairs, called_seq):
                ev["ref_position"] = pairs[1][ev_index[0]] #first position of kmer in reference
                correct.append(ev)
        return correct, total

    def _event_indexes(self, pairing_seq, offset):
        """ get the next entries from the pairing array
        such that k non-gap characters are contained"""
        count = 0
        kmer = []
        for i in range(offset, len(pairing_seq)):
            if count == self.k: break
            if pairing_seq[i] is not None:
                count += 1
                kmer.append(i)
        if len(kmer) != self.k:
            raise AlignmentEndException
        return kmer

    @staticmethod
    def _gapmove(to_move, seq, offset):
        """move by 'move' (from metrichor) in the aligned sequence.
        additionally increase index to compensate for each gap
        """
        move = to_move
        for i in seq[offset:]:
            if i is None:
                move += 1
            else:
                to_move -= 1
                if to_move <= 0:
                    return move

    @staticmethod
    def _get_nt_kmer(index, pairs, seq):
        """convert sequence indexes into the corresponding nucleotides.
        gaps are converted into ''
        """
        seq_index = [pairs[x] for x in index]
        nt_kmer = [seq[x] for x in seq_index]
        return "".join(nt_kmer)

    @staticmethod
    def _is_consecutive_seq(seq):
        """check if the sequence 'seq' consists of consecutive numbers"""
        return len(set(list(map(lambda ix:ix[1]-ix[0], enumerate(seq))))) <= 1

    def _is_correct_kmer(self, ev_index, pairs, read):
        """check if a kmer corresponds completely wit the reference.
        This is the case if:
            * the read positions are consecutive (no indels)
            * the ref positions are consecutive (no indels)
            * the nucleotides are idential (no substitutions)
        """
        assert(len(ev_index) == self.k), "invalid event index"
        read_index = [pairs[0][x] for x in ev_index]
        ref_index = [pairs[1][x] for x in ev_index]

        if None in read_index or not ModelMaker._is_consecutive_seq(read_index):
            """indel in read"""
            return False

        if None in ref_index or not ModelMaker._is_consecutive_seq(ref_index):
            """indel in ref"""
            return False

        read_seq = [read[x] for x in read_index]
        ref_seq = [self.ref_seq[x] for x in ref_index]
        if read_seq == ref_seq:
            """full_match"""
            return True
        else:
            """substitution"""
            return False


if __name__ == "__main__":
    argp = argparse.ArgumentParser("Align processed reads from metrichor to reference; "
                                   "extract correct kmers and calculate mean and stdv for the model. "
                                   "Only data from 2D reads is used."
                                   ""
                                   "Warning. This script loads all fast5-files to RAM. Ensure that you have enough "
                                   "or do not use all the data. ")
    argp.add_argument("-f", "--filelist", required=True, type=argparse.FileType('r'),
                      help="a list of fast5 files to be aligned, one per line")
    argp.add_argument("-r", "--reference", required=True, type=argparse.FileType('r'),
                      help="fasta file with reference sequence")
    argp.add_argument("-o", "--output", required=True, type=str,
                      help="model basename. /path/to/my_model will result in e.g. /path/to/my_model.template.pickle")
    argp.add_argument("-n", "--ncores", required=False, type=int,
                      help="#CPU cores", default=None)
    argp.add_argument("-k", "--kmer", required=False, type=int,
                      help="length of kmer", default=6)
    argp.add_argument("-g", "--graphmap", required=False, type=str, default="",
                       help="Path to graphmap alignment tool")

    args = argp.parse_args()
    ModelMaker(args.filelist, args.reference, args.output, args.ncores, args.kmer, args.graphmap)