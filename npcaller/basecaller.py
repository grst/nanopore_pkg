import pickle
import ghmm
from npcaller.fast5 import Fast5File
from multiprocessing import Pool
import sys
from npcaller.fasta import FastaWriter
import numpy as np

class ModelException:
    pass

class HmmModel(object):
    TRANSMAT_CONST_MOVE = HmmModel.mk_transmat1
    TRANSMAT_MOVE_MAX_2 = HmmModel.mk_transmat2
    emission_domain = ghmm.Float() # emission domain for HMM
    model = None

    def __init__(self, model_file, f_transmat=HmmModel.TRANSMAT_MOVE_MAX_2):
        self.model_params = pickle.load(model_file)
        self.nmers = len(self.model_params[["kmer"]][0])
        assert 4 ** self.nmers == len(self.model_params), "# model params do not match length of kmers"
        self.nstates = len(self.model_params)
        self.init_model(f_transmat)

    def init_model(self, f_transmat):
        """
        Build a HMM according to the parameter File.

        The model is only based on the means of events. Stdv is not taken into account.

        Args:
            f_transmat: function to construct the transition matrix

        Returns: GaussianEmissionHMM
        """
        A = self.f_transmat()
        B = self.model_params[["level_mean", "level_stdv"]].values.tolist() #mu, std of each state
        pi = [1/float(self.nstates)] * self.nstates   # initial probabilities per state
        # generate model from parameters
        self.model = ghmm.HMMFromMatrices(self.emission_domain,
                                          ghmm.GaussianDistribution(self.emission_domain), A, B, pi)


    def mk_transmat1(self):
        """make a transition matrix assuming move=1"""
        transmat = np.empty((self.nstates, self.nstates))
        for j, from_kmer in enumerate(self.model_params[["kmer"]]):
            for i, to_kmer in enumerate(self.model_params[["kmer"]]):
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
        for j, from_kmer in enumerate(self.model_params[["kmer"]]):
            for i, to_kmer in enumerate(self.model_params[["kmer"]]):
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
        kmers = [self.model_params[["kmer"]][x] for x in states]
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
            f5_files: list of fast5 files
            ncores: number of CPU cores used. If ncores is None, then cpu_count() is used.
            template_model: template model file
            complement_model: complement model file

        Returns:

        """
        if template_model is None and complement_model is None:
            raise ModelException("no model specified")

        self.ncores = ncores
        self.template_model = HmmModel(template_model) if template_model else None
        self.complement_model = HmmModel(complement_model) if complement_model else None
        self.results = []

    def process_files(self, files):
        def file_worker(file_path):
            template, complement = None, None
            f5file = Fast5File(file_path)
            template_events = f5file.get_events("template")
            complement_events = f5file.get_events("complement")
            if self.template_model and template_events:
                means = [ev["mean"] for ev in template_events]
                template = self.template_model.predict(means)
            if self.complement_model and complement_events:
                means = [ev["mean"] for ev in complement_events]
                complement = self.complement_model.predict(means)
            return file_path, template, complement

        p = Pool(self.ncores)
        try:
            print("Calling Files: ")
            for i, res in enumerate(p.imap_unordered(file_worker, files), 1):
                self.results.append(res)
                sys.stdout.write('\rdone {0:%}'.format(i/float(len(files))))
            p.close()
            p.join()
        except KeyboardInterrupt:
            p.terminate()

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



