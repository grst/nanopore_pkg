#!/usr/bin/env python3
import argparse
from npcaller.fasta import FastaReader, FastaWriter
from npcaller.fast5 import Fast5File
from npcaller.validator import align_to_reference, sam_to_bam
from tempfile import NamedTemporaryFile
from pysam import AlignmentFile
import numpy as np
import pandas
import logging
from multiprocessing import Pool
import itertools
import pickle
logging.basicConfig(level=logging.INFO)


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

        self.logger = logging.getLogger("alignment")
        self.logger.setLevel(logging.NOTSET)

        fr = FastaReader(ref_file)
        _, self.ref_seq = next(fr.get_entries())

        self.f5files = {}
        filelist = [l.strip() for l in filelist.readlines()]
        for file in filelist:
            f5file = Fast5File(file)
            self.f5files[f5file.get_id()] = f5file

    def make_all_models(self):
        """
        Make the models for both template and complement.
        """
        for strand in ["template", "complement"]:
            self.make_model(strand)

    def make_model(self, strand):
        """
        Runner function: make the model for either template or complement strand.
        Pickles the model to {basename}.strand.pickle.
        Args:
            strand:

        Returns:

        """
        self.logger.info("Making model for {0}-strand".format(strand))
        alignment = self._make_bam(strand)
        correct_kmers = self._find_correct_kmers(alignment, strand)
        model = self._make_stats(correct_kmers)
        pickle.dump(model, open("{0}.{1}.{2}".format(self.model_basename, strand, "pickle"), 'wb'), protocol=2)

    def _make_bam(self, strand):
        """
        Extract sequences from fast5 files and map them to the reference sequence

        Args:
            strand: either template or complement

        Returns:
            AlignmentFile object of the generated bam-Alignment

        """
        # with NamedTemporaryFile() as fasta_file:
        with NamedTemporaryFile('w') as fasta_file:
            samfile = self.model_basename + ".sam"
            bamfile = self.model_basename
            fw = FastaWriter(fasta_file)
            for f5file in self.f5files.values():
                header = f5file.get_id()
                seq = f5file.get_seq(strand)
                fw.write_entry(header, seq)
            fw.flush()
            align_to_reference(fasta_file, self.ref_file, samfile,
                               graphmap_bin=self.graphmap_bin, ncores=self.ncores)
            sam_to_bam(samfile, bamfile)

        return AlignmentFile(bamfile + ".bam")

    def _find_correct_kmers(self, alignment, strand):
        """

        Args:
            alignment (AlignmentFile): Pysam Object of the bam-Alignment
            strand (str): either template or complement

        Returns:
            list of correctly mapped events.
        """
        total_events = 0
        result = list()
        reads = [x for x in alignment.fetch()]
        self.logger.info("{0} reads found in alignment".format(len(reads)))
        for read in reads:
            f5file = self.f5files[read.query_name]
            pairs = [list(t) for t in zip(*read.get_aligned_pairs())]
            assert(pairs[0][0] == 0), "alignment is not null-indexed."
            correct, total = self._process_events(f5file, pairs, strand)
            total_events += total
            result.append(correct)

#        true_events = list(chain.from_iterable([r.get() for r in result]))
        correct_events = [x for x in itertools.chain.from_iterable(result)]
        self.logger.info("Identified {0} correct kmers of {1} total kmers. That's {2:%}".format(
            len(correct_events), total_events, len(correct_events)/total_events
        ))
        return correct_events

    def _make_stats(self, correct_events):
        """
        sort the events into their respective kmer-buckets and calculate the target
        statistics (mean, sd) for the model

        Args:
            correct_events (list): list of correctly mapped events

        Returns:
            Pandas Dataframe containing the model

        """
        self.logger.info("started calculating statistics")
        all_kmers = ["".join(i) for i in itertools.product("ACGT", repeat=self.k)]
        stat_map = {}
        for attr in ["mean", "stdv"]:
            stat_map[attr] = {kmer: [] for kmer in all_kmers}
            for ev in correct_events:
                stat_map[attr][ev["kmer"]].append(ev[attr])

        # make model file
        model = []
        for kmer in all_kmers:
            model.append({
                "kmer": kmer,
                "level_mean": np.mean(stat_map["mean"][kmer]),
                "level_stdv": np.std(stat_map["mean"][kmer]),
                "sd_mean": np.mean(stat_map["stdv"][kmer]),
                "sd_stdv": np.std(stat_map["stdv"][kmer]),
                "weight": 1000.0  # not implemented in this model, use neutral value
            })
        return pandas.DataFrame(model)

    def _process_events(self, f5file, pairs, strand):
        """
        Helper function which processes the events per f5file.

        Args:
            f5file (Fast5File):
            pairs (list): list of pairs (read nt <-> ref nt)
            strand: either template or complement

        Returns:
            correctly mapped events of the given file.

        """
        i_seq = 0
        correct = list()
        total = 0
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
                ev["ref_position"] = pairs[1][ev_index[0]]  # first position of kmer in reference
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
                                   "Only data from 2D reads is used.")
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
    mm = ModelMaker(args.filelist, args.reference, args.output, args.ncores, args.kmer, args.graphmap)
    mm.make_all_models()