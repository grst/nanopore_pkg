"""
toolbox for handling fast5 (hdf5) files.
"""

import h5py
import re

class lazy_property(object):
    '''
    meant to be used for lazy evaluation of an object attribute.
    property should represent non-mutable data, as it replaces itself.
    '''

    def __init__(self,fget):
        self.fget = fget
        self.func_name = fget.__name__


    def __get__(self,obj,cls):
        if obj is None:
            return None
        value = self.fget(obj)
        setattr(obj,self.func_name,value)
        return value

class Fast5Exception(Exception):
    pass


class Fast5PathException(Fast5Exception):
    pass


def identify_events_path(fast5):
    """
    Newer versions of the a fast5 file have a different path to the event table.
    Identify which type of fast5 file we have and return the path accordingly.

    Args:
        fast5: h5py File object

    Returns:
        corresponding index of the fast5 file
    """

    try:
        _ = fast5["/Analyses/Basecall_1D_000/BaseCalled_template"]
        return "/Analyses/Basecall_1D_000/BaseCalled_{0}"
    except KeyError:
        try:
            _ = fast5["/Analyses/Basecall_2D_000/BaseCalled_template"]
            return "/Analyses/Basecall_2D_000/BaseCalled_{0}"
        except KeyError:
            raise Fast5PathException("Unknown version of fast5 format")


class Fast5File(object):
    def __init__(self, path):
        result = re.search(r'ch(\d+)_file(\d+)_', path)
        self.file_id = int(result.group(2))
        self.channel_id = int(result.group(1))
        self.seqs = {
            "template": None,
            "complement": None
        }
        try:
            self.f5 = h5py.File(path, 'r')
            self.event_path = identify_events_path(self.f5)
        except OSError:
            raise Fast5Exception("Unable to open fast5-file.")

    def get_id(self):
        """
        Returns: unique identifier for f5 file.
        """
        return "ch_{0}_file_{1}".format(self.channel_id, self.file_id)

    def get_attrs(self, strand):
        """
        Get Attributes for template/complement strand

        Args:
            strand: either "template" or "complement"

        Returns:
            {shift: foo, drift: foo, scale: foo} or None if strand not in File

        """
        assert strand in ["template", "complement"]
        attrs = ["shift", "scale", "drift"]
        try:
            all_attrs = self.f5[self.event_path.format(strand) + "/Model"].attrs
            return {k: all_attrs[k] for k in attrs}
        except (KeyError, ValueError):
            return None

    def get_events(self, strand):
        """
        Get events for template/complement strand

        Args:
            strand: either "template" or "complement"

        Returns:
            generator yielding one event-dict at a time or None if strand not in File

        """
        assert strand in ["template", "complement"]
        try:
            events = self.f5[self.event_path.format(strand) + "/Events"]
            for raw_ev in events:
                ev = {}
                ev["mean"] = raw_ev[0]
                ev["start"] = float(raw_ev[1])
                ev["stdv"] = raw_ev[2]
                ev["length"] = float(raw_ev[3])
                ev["kmer"] = bytes(raw_ev[4]).decode('utf-8')
                ev["move"] = int(raw_ev[6])
                yield ev
        except KeyError:
            pass

    def get_corrected_events(self, strand):
        """
        Get events for template/complement strand and apply the shift/scale/drift corrections.
        Unfortunately, these corrections are not documented exactly anywhere. Some information is on
        https://wiki.nanoporetech.com/display/BP/1D+Basecalling+overview.

        Args:
            strand: either "template" or "complement"

        Returns:
            generator yielding one event-dict at a time or None if strand not in File

        """
        assert strand in ["template", "complement"]
        attrs = self.get_attrs(strand)
        if attrs is None:
            raise Fast5Exception("events cannot be shift/scale/drift corrected.")

        # apparently, drift is applied relatively to the start of the read.
        read_start = next(self.get_events(strand))["start"]

        shift = lambda mean: mean - attrs["shift"]
        scale = lambda mean: mean / attrs["scale"]
        drift = lambda mean, start: mean - (start - read_start) * attrs["drift"]

        for ev in self.get_events(strand):
            ev["mean"] = drift(scale(shift(ev["mean"])), ev["start"])
            yield ev

    def get_seq(self, strand):
        """
        get the nt-sequence, based on the kmers and 'move'.
        Evaluation is done in a lazy fashion. If the function
        was called once, the sequence can be accessed in constant time.

        Args:
            strand (str): either template or complement

        Returns:
            (str) nucleotide sequence

        """
        assert strand in ["template", "complement"]
        if self.seqs[strand] is None:
            self.seqs[strand] = self.events2seq(self.get_events(strand))
        return self.seqs[strand]

    @staticmethod
    def events2seq(events):
        """turn events into a nt sequence based on the 'move' column

        Args:
            events: list of events (dictionaries)
        """
        seq = list()
        seq.append(next(events)["kmer"])
        for ev in events:
            move = ev["move"]
            kmer = ev["kmer"]
            seq.append(kmer[len(kmer)-move:])
        return "".join(seq)


