from __future__ import print_function
import pickle
import ghmm
from npcaller.fast5 import Fast5File
from multiprocessing import Pool
import sys
from npcaller.fasta import FastaWriter
import numpy as np

# The GHMM object cannot be pickled (SWIG object)
# To use the multiprocessing library, a workaround is to save the objects in
# global variables. This is ugly, but works for this application.
template_model = None
complement_model = None


def call_file(file_path):
    """
    Execute basecalling on the given file.

    This is the method called by the multiprocessing.Pool.map() and executed in
    one worker thread. Should be part of the Basecaller class, but the
    multiprocessing module does not allow that (Pickling error...)

    Args:
        file_path: path to fast5 file

    Returns:
        Tuple (file_path, template_sequence, complement_sequence)

    """
    template, complement = None, None
    f5file = Fast5File(file_path)
    template_events = f5file.get_corrected_events("template")
    complement_events = f5file.get_corrected_events("complement")
    if template_model and template_events:
        means = [ev["mean"] for ev in template_events]
        template = template_model.predict(means)
    if complement_model and complement_events:
        means = [ev["mean"] for ev in complement_events]
        complement = complement_model.predict(means)
    return file_path, template, complement


class ModelException(Exception):
    pass


class HmmModel(object):
    TRANSMAT_CONST_MOVE = 0
    TRANSMAT_MOVE_MAX_2 = 2
    emission_domain = ghmm.Float()      # emission domain for HMM
    model = None

    def __init__(self, model_file, which_transmat=TRANSMAT_MOVE_MAX_2):
        self.model_params = pickle.load(model_file)
        self.nmers = len(self.model_params.ix[0, "kmer"])
        self.all_kmers = self.model_params["kmer"].tolist()
        assert 4 ** self.nmers == len(self.model_params), "# model params do not match length of kmers"
        self.nstates = len(self.model_params)
        self.init_model(which_transmat)

    def init_model(self, which_transmat):
        """
        Build a HMM according to the parameter File.

        The model is only based on the means of events. Stdv is not taken into account.

        Args:
            which_transmat: HmmModel.TRANSMAT_CONST_MOVE or HmmModel.TRANSMAT_MOVE_MAX_2

        Returns: GaussianEmissionHMM
        """
        mat_gen = {
            self.TRANSMAT_MOVE_MAX_2: self.mk_transmat2,
            self.TRANSMAT_CONST_MOVE: self.mk_transmat1
        }
        A = mat_gen[which_transmat]()
        B = self.model_params[["level_mean", "level_stdv"]].values.tolist()  # mu, std of each state
        pi = [1/float(self.nstates)] * self.nstates   # initial probabilities per state
        # generate model from parameters
        self.model = ghmm.HMMFromMatrices(self.emission_domain,
                                          ghmm.GaussianDistribution(self.emission_domain), A, B, pi)

    def mk_transmat1(self):
        """make a transition matrix assuming move=1"""
        transmat = np.empty((self.nstates, self.nstates))
        for j, from_kmer in enumerate(self.all_kmers):
            for i, to_kmer in enumerate(self.all_kmers):
                p = 1/4. if from_kmer[-(self.nmers-1):] == to_kmer[:(self.nmers-1)] else 0.
                transmat[j, i] = p

        return transmat.tolist()

    def mk_transmat2(self):
        """
        make a transition matrix assuming move=0 or move=1 or move=2

        The fractions are chosen by trial and error. Performance can probably be gained by choosing
        them systematically.
        """
        transmat = np.empty((self.nstates, self.nstates))
        for j, from_kmer in enumerate(self.all_kmers):
            for i, to_kmer in enumerate(self.all_kmers):
                p = 0
                if from_kmer[-(self.nmers-2):] == to_kmer[:(self.nmers-2)]:
                    """move=2"""
                    p = (2/50.) * (1/16.)
                elif from_kmer[-(self.nmers-1):] == to_kmer[:(self.nmers-1)]:
                    """move=1"""
                    p = (47/50.) * (1/4.)
                elif from_kmer == to_kmer:
                    """move=0"""
                    p = (1/50.) * 1
                transmat[j, i] = p

        return transmat.tolist()

    def path_to_seq(self, states):
        """

        Args:
            states: list of states (result of viterbi algorithm)

        Returns:
            nucleotide sequence

        """
        kmers = [self.model_params.ix[x, "kmer"] for x in states]
        seq = [kmer[0] for kmer in kmers] + [kmers[-1][1:]]
        return "".join(seq)

    def predict(self, means):
        """
        Predict the most likely path of states for the given means (viterbi)

        Args:
            means: list of means

        Returns:
            seq: called nucleotide sequence

        """
        seq = ghmm.EmissionSequence(self.emission_domain, means)
        result = self.model.viterbi(seq)
        states = result[0]
        return self.path_to_seq(states)


class Basecaller(object):
    def __init__(self, ncores=None, template_model=None, complement_model=None):
        """
        Args:
            ncores: number of CPU cores used. If ncores is None, then cpu_count() is used.
            template_model: template model file
            complement_model: complement model file

        """
        if template_model is None and complement_model is None:
            raise ModelException("no model specified")

        self.ncores = ncores
        self.template_model = HmmModel(template_model) if template_model else None
        self.complement_model = HmmModel(complement_model) if complement_model else None
        self.results = []

    def process_files(self, files):
        """
        Execute multi-threaded basecalling on a list of fast5-files

        Args:
            files: list of paths to fast5 files

        """
        # Global variables as workaround to share unpicklable objects with the worker threads.
        global template_model
        global complement_model
        template_model = self.template_model
        complement_model = self.complement_model

        p = Pool(self.ncores)
        try:
            print("Calling Files: ")
            for i, res in enumerate(p.imap_unordered(call_file, files), 1):
                self.results.append(res)
                sys.stdout.write('\rdone {0:%}'.format(i/float(len(files))))
                sys.stdout.flush()
            p.close()
            p.join()
        except KeyboardInterrupt:
            p.terminate()
            raise KeyboardInterrupt

    def write_to_fasta(self, output_file):
        """
        write the results to a fasta file

        Args:
            output_file (file): file object to write results

        """
        fw = FastaWriter(output_file)
        for file_path, template, complement in self.results:
            if template:
                header = "{0} {1}".format(file_path, "template")
                fw.write_entry(header, template)
            if complement:
                header = "{0} {1}".format(file_path, "complement")
                fw.write_entry(header, complement)
