"""
toolbox for handling fast5 (hdf5) files.
"""

import h5py


class Fast5PathException(Exception):
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

