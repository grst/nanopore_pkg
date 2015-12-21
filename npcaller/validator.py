from fasta import FastaReader, FastaWriter
import random
import pysam
from subprocess import call
import math
import csv
import itertools
from multiprocessing import Pool

def randomize_fasta(fasta_in, fasta_out):
    """
    read fasta_in, keep the sequence lengths, but mutate all nucleotides randomly, assuming a
    uniform distribution.

    Args:
        fasta_in: path or file handle to fasta input file
        fasta_out: path or file handle to fasta output file

    """
    fr = FastaReader(fasta_in)
    fw = FastaWriter(fasta_out)
    for header, seq in fr.get_entries():
        header += " random"
        seq = "".join(random.choice("ACGT") for _ in seq)
        fw.write_entry(header, seq)
    fr.close()
    fw.close()


def align_to_reference(fasta_file, ref_file, sam_file, graphmap_bin="graphmap", ncores="-1"):
    """
    Align reads to the reference using graphmap.

    Args:
        fasta_file: fasta file with reads to align
        ref_file: fasta file with reference (genome) sequence
        sam_file: output file
        graphmap_bin: path to graphmap binary
        ncores: number of cpu cores to use for alignment. If ncores=-1, ncores equals min(24, cpu_count()/2)

    """
    call([graphmap_bin,
          "-r", ref_file.name,
          "-d", fasta_file.name,
          "-o", sam_file.name,
          "-t", ncores])


def _mk_read_stats(samfile_path, read_id, total_reads, ref_seq, chunk):
    """
    make statistics for a chunk of reads. This function is the worker function for multiprocessing.Pool.map().

    Since the Pysam AlignedRead is not thread-safe and induces a deadlock, every worker has to open its own
    AlignmentFile Object and skip all reads until the read of interest is reached.

    This should be part of AlignmentStats, but the native multiprocessing library does not allow that (pickling error)

    Args:
        read: pysam.AlignedRead
    """
    samfile = pysam.AlignmentFile(samfile_path)
    # since reads can only be accessed though the generator, we have to skip those we are not interested in:
    generator = samfile.fetch()
    try:
        read = None
        for i in range(read_id):            # skip previous reads
            read = generator.next()

        for i in range(chunk):              # process the reads of interest
            read = generator.next()

            s = {
                "read_len": read.get_tag("ZQ"),
                "mapping_quality": read.mapping_quality,
                "aln_len": read.alen,
                "aln_score": read.get_tag("AS"),
                "aln_evalue": read.get_tag("ZE"),
                "aln_editdistance": read.get_tag("NM"),
                "mapped_nts": sum([l for op, l in read.cigartuples if op == 0]),
                "ins": sum(l for op, l in read.cigartuples if op == 1),
                "del": sum(l for op, l in read.cigartuples if op == 2),
                "subst": sum(1 for query, ref in read.get_aligned_pairs(matches_only=True)
                             if read.seq[query] != ref_seq[ref]),
                "is_significant": AlignmentStats.e2p(read.get_tag("ZE"))
                                  < AlignmentStats.SIGNIFICANCE / total_reads   # bonferroni correction
            }
            assert s["read_len"] >= s["mapped_nts"]
            assert s["aln_len"] >= s["mapped_nts"]
            assert len(read.query_sequence) == s["read_len"]
            yield s

    except StopIteration:
        pass


class AlignmentStats(object):
    """
    Generate statistics for a SAM-file.
    """

    SIGNIFICANCE = .05

    def __init__(self,  sam_file, ref_seq, ncores=-1):
        """
        Instantiate object and calculate statistics.

        Args:
            sam_file (File): File handle of SAM-file
            ref_seq (str): reference sequence
            ncores (int): # cpu cores
        """
        self.ncores = ncores
        self.sam_file = sam_file
        self.pysam_file = pysam.AlignmentFile(sam_file.name)
        self.ref_seq = ref_seq
        self.file_stats = {}
        self.read_stats = []
        self.summary = []
        self.mk_stats()
        self.mk_summary()

    @staticmethod
    def e2p(e):
        """
        convert Evalue to Pvalue (of Alignment)
        """
        return 1-np.exp(-e)

    def mk_stats(self):
        """
        Do the actual calculations...
        """
        samtools_result = call(["samtools view {0} | grep -v '^#' | cut -f17 | cut -d':' -f3".format(
            self.sam_file.name)], shell=True) # use samtools and shell tools to the the #lines and #total nts
        self.file_stats["total_reads"] =  len(samtools_result),   #for some reason pysam.unmapped does not work...
        self.file_stats["mapped_reads"] =  self.pysam_file.mapped,
        self.file_stats["total_nts"] = sum(int(x) for x in samtools_result)

        chunk = int(math.ceil(self.file_stats["mapped_reads"]/float(self.ncores))) #chunk size for multiprocessing
        p = Pool(self.ncores)
        try:
            read_stats = p.map(_mk_read_stats,
                               itertools.repeat(self.sam_file.name),
                               range(0, self.file_stats["mapped_reads"], chunk),
                               itertools.repeat(self.file_stats["total_reads"]),
                               itertools.repeat(self.ref_seq),
                               itertools.repeat(chunk))
            p.close()
        except KeyboardInterrupt:
            p.terminate()

        self.read_stats = list(itertools.chain(*read_stats))

     def sumstat(self, stat):
        """
        sum of the stat <stat> for all reads.

        Args:
            stat: statistics to sum up.

        """
        return sum(r[stat] for r in self.read_stats)

    def mk_summary(self):
        lines = [
            ["mapped_reads/total_reads", self.file_stats["mapped_reads"], self.file_stats["total_reads"]],
            ["significant_reads/total_reads", self.sumstat("is_significant"), self.file_stats["total_reads"]],
            ["alignment_score/alignment_length", self.sumstat("aln_score"), self.sumstat("aln_len")],
            ["mapped_nts/total_nts", self.sumstat("mapped_nts"), self.file_stats["total_nts"]],
            ["editdistance/alignment_length", self.sumstat("aln_editdistance"), self.sumstat("aln_len")],
            ["SNPs/mapped_nts", self.sumstat("subst"), self.sumstat("mapped_nts")],
            ["ins/mapped_nts", self.sumstat("ins"), self.sumstat("mapped_nts")],
            ["del/mapped_nts", self.sumstat("del"), self.sumstat("mapped_nts")],
        ]

        def process_line(line):
            return [str(x) for x in (line + ["{0:%}".format(float(line[1])/line[2])])]

        self.summary = [process_line(line) for line in lines]

    def to_tsv(self, output_file):
        """
        write summary to tsv file.

        Args:
            output_file (File 'w'):

        """
        csvw = csv.writer(output_file, delimiter="\t", quoting=csv.QUOTE_NONE)
        for row in self.summary:
            csvw.writerow(row)
