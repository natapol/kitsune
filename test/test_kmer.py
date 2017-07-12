# -*- coding: utf-8 -*-

import ksiga.kmerutil as kmer


def test_convertKmerToLocation():
    """ Produce a correct insert location.
    """
    fn = kmer.create_kmer_loc_fn(7)

    assert fn("AAAAAA" + "A") == 0
    assert fn("AAAAAA" + "C") == 1
    assert fn("AAAAAA" + "G") == 2
    assert fn("AAAAAA" + "T") == 3