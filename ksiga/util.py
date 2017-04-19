# -*- coding: utf-8 -*-

import gzip


def copen(filename, mode="r", encoding="utf-8", **kwargs):
    """Conveniant method for openning file.

    Args:
        filename (str): filename

    Returns: TODO

    """
    if filename.endswith(".gz"):
        return gzip.open(filename, mode=mode, encoding=encoding, **kwargs)
    else:
        return open(filename, mode, encoding=encoding, **kwargs)


def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
