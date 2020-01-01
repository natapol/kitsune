# -*- coding: utf-8 -*-

##### Setup
from setuptools import setup, find_packages

def readme():
    with open('README.rst') as f:
        return f.read()

setup(
    name='kitsune',
    version='1.2.6',
    description='tools for finding an optimal kmer',
    long_description=readme(),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
    ],
    keywords='kitsune',
    url='https://gitlab.com/natapol.por/ksiga',
    packages=['kitsune'],
    author='Natapol Pornputtapong',
    author_email='natapol.por@gmail.com',
    license="Apache 2.0",
    scripts=[],
    entry_points = {
        'console_scripts': ['kitsune=kitsune.main:main'],
    },
    include_package_data=True,
    install_requires=["numpy >= 1.1.0",
                        "scipy>=0.18.1",
                        "biopython>=1.68",
                        "tqdm>=4.32"
                    ],
    zip_safe=False
)
