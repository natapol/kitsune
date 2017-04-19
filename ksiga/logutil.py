# -*- coding: utf-8 -*-

""" Provide utility functions relate to logging.
"""

import sys
import os

QUIET = False
handle = sys.stderr


def warn(msg):
    """TODO: Docstring for warn.

    Args:
        msg (TODO): TODO

    Returns: TODO

    """
    handle.write(msg)
    handle.write(os.linesep)


def notify(msg):
    """TODO: Docstring for warn.

    Args:
        msg (TODO): TODO

    Returns: TODO

    """
    if QUIET:
        return

    if type(msg) is not str:
        msg = str(msg)
    handle.write(msg)
    handle.write(os.linesep)
