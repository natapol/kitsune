# -*- coding: utf-8 -*-

""" Log-utility
"""

import sys
import os

QUITE = False
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
    if not QUITE:
        if type(msg) is not str:
            msg = str(msg)
        handle.write(msg)
        handle.write(os.linesep)
