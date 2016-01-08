from npcaller.fasta import FastaReader, FastaWriter
import random
import pysam
from subprocess import call, check_output
import math
import numpy as np
import csv
import itertools
import tempfile
from multiprocessing import Pool, cpu_count


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


def filter_fasta(fasta_in, fasta_out, keyword):
    """
    read fasta_in, keep all sequences that contain the keyword in the header

    Args:
        fasta_in: path or file handle to fasta input file
        fasta_out: path or file handle to fasta output file
        keyword: filter criterion

    """
    fr = FastaReader(fasta_in)
    fw = FastaWriter(fasta_out)
    for header, seq in fr.get_entries():
        if header.find(keyword) >= 0:
            fw.write_entry(header, seq)
    fr.close()
    fw.close()


def align_to_reference(fasta_file, ref_file, samfile_name, graphmap_bin="graphmap", ncores="-1"):
    """
    Align reads to the reference using graphmap.
    Create SAM and BAM file. Sort and Index Samfiles.

    Args:
        fasta_file: fasta file with reads to align
        ref_file: fasta file with reference (genome) sequence
        samfile_name: name of sam file
        graphmap_bin: path to graphmap binary
        ncores: number of cpu cores to use for alignment. If ncores=-1, ncores equals min(24, cpu_count()/2)

    """
    call([graphmap_bin,
          "-r", ref_file.name,
          "-d", fasta_file.name,
          "-o", samfile_name,
          "-t", str(ncores)])


def sam_to_bam(samfile, bamfile, samtools_bin="samtools"):
    """
    Convert SAM to sorted and indexed BAM.

    Args:
        samfile: Input Sam
        bamfile: Output Bam
        samtools_bin: path to samtools binary

    """
    with tempfile.NamedTemporaryFile() as tmp_bam:
        # samtools takes only basename in 'sort' but full filename in 'index'
        call(["{2} view -S -b {0} > {1}".format(samfile, tmp_bam.name, samtools_bin)], shell=True)
        call([samtools_bin, "sort", tmp_bam.name, bamfile])
        call([samtools_bin, "index", bamfile + ".bam"])


def _mk_read_stats(samfile_path, read_id, total_reads, ref_seq, chunk):
    """
    make statistics for a chunk of reads. This function is the worker function for multiprocessing.Pool.map().

    Since the Pysam AlignedRead is not thread-safe and induces a deadlock, every worker has to open its own
    AlignmentFile Object and skip all reads until the read of interest is reached.

    This should be part of AlignmentStats, but the native multiprocessing library does not allow that (pickling error)

    Args:

    """
    result = []
    samfile = pysam.AlignmentFile(samfile_path)
    # since reads can only be accessed though the generator, we have to skip those we are not interested in:
    generator = samfile.fetch()
    try:
        read = None
        for i in range(read_id):            # skip previous reads
            read = next(generator)

        for i in range(chunk):              # process the reads of interest
            read = next(generator)

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
            result.append(s)

    except StopIteration:
        pass

    return result


class AlignmentStats(object):
    """
    Generate statistics for a SAM-file.
    """

    SIGNIFICANCE = .05

    def __init__(self, bam_file, ref_seq, ncores=None):
        """
        Instantiate object and calculate statistics.

        Args:
            bam_file (File): File handle of BAM-file
            ref_seq (str): reference sequence
            ncores (int): # cpu cores
        """
        self.ncores = cpu_count() if not ncores else ncores
        self.bam_file = bam_file
        self.pysam_file = pysam.AlignmentFile(bam_file.name)
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

        Args:
            e: evalue
        """
        return 1-np.exp(-e)

    def mk_stats(self):
        """
        Do the actual calculations...
        """
        samtools_result = check_output(["samtools view {0} | grep -v '^#' | cut -f17 | cut -d':' -f3".format(
            self.bam_file.name)], shell=True)   # use samtools and shell tools to the the #lines and #total nts
        samtools_result = samtools_result.decode("utf-8").strip().split("\n")
        self.file_stats["total_reads"] = len(samtools_result)   # for some reason pysam.unmapped does not work...
        self.file_stats["mapped_reads"] = self.pysam_file.mapped
        self.file_stats["total_nts"] = sum(int(x) for x in samtools_result)

        chunk = int(math.ceil(self.file_stats["mapped_reads"]/float(self.ncores)))  # chunk size for multiprocessing
        p = Pool(self.ncores)
        try:
            read_stats = p.starmap(_mk_read_stats, zip(
                               itertools.repeat(self.bam_file.name),
                               range(0, self.file_stats["mapped_reads"], chunk),
                               itertools.repeat(self.file_stats["total_reads"]),
                               itertools.repeat(self.ref_seq),
                               itertools.repeat(chunk)))
            p.close()
            self.read_stats = list(itertools.chain(*read_stats))

        except KeyboardInterrupt:
            p.terminate()

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
