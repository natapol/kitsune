# -*- coding: utf-8 -*-

""" Logging utility"""

import os, errno
import logging
from logging.handlers import RotatingFileHandler

DEFAULT_FORMAT = logging.Formatter('%(asctime)s - '
                              '%(filename)s:%(funcName)s '
                              '%(levelname)s - '
                              '%(lineno)d:\t'
                              '- %(message)s')

DEFAULT_LOGGER = logging.getLogger("ksiga")
DEFAULT_LOGGER.addHandler(logging.StreamHandler())
DEFAULT_LOGGER.setLevel(logging.INFO)

def getLogger(name, level=logging.INFO, logformat=DEFAULT_FORMAT):

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    sh = logging.StreamHandler()
    sh.setLevel(level)
    sh.setFormatter(logformat)

    # try:
        # os.makedirs("logs")
    # except OSError as e:
        # if e.errno == errno.EEXIST and os.path.isdir("logs"):
            # pass
        # else:
            # raise

    # fh = RotatingFileHandler("logs/kali.log", maxBytes=21474836480, backupCount=5)  # 20 MB

    # fh.setLevel(logging.DEBUG)
    # fh.setFormatter(logformat)

    logger.addHandler(sh)
    # logger.addHandler(fh)

    return logger


def conditional_logging(logger, method, message):
    """Log only if the global variable is available (LOGGING=TRUE)"""
    global LOGGING

    if LOGGING is True:
        logging_method = getattr(logger, method)
        logging_method(message)


def debug(message):
    DEFAULT_LOGGER.debu(message)
    pass


def info(message):
    DEFAULT_LOGGER.info(message)


def warning(message):
    DEFAULT_LOGGER.warning(message)


def notify(message):
    DEFAULT_LOGGER.notify(message)