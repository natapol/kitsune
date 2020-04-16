# -*- coding: utf-8 -*-

import re
import os

##### Setup
from setuptools import setup

def readme():
    """
    create logdescription from rst file
    """
    with open('README.rst') as f:
        return f.read()

# The full version, including alpha/beta/rc tags.
release = re.sub('^v', '', os.popen('git describe --tags').read().strip())

setup(
    name='kitsune',
    version=release,
    description='a toolkit for evaluation of the lenght of k-mer in a given genome dataset for alignment-free phylogenimic analysis',
    long_description=readme(),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
    ],
    keywords='kitsune',
    url='https://github.com/natapol/kitsune',
    packages=['kitsune'],
    author='Natapol Pornputtapong',
    author_email='natapol.por@gmail.com',
    scripts=[],
    entry_points={
        'console_scripts': ['kitsune=kitsune.main:main'],
    },
    include_package_data=True,
    install_requires=[
        "numpy >= 1.1.0",
        "scipy>=0.18.1",
        "biopython>=1.68",
        "tqdm>=4.32"
    ],
    zip_safe=False
)
