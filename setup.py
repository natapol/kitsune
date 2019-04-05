# -*- coding: utf-8 -*-

##### Setup
from setuptools import setup, find_packages

def readme():
    with open('README.rst') as f:
        return f.read()

setup(
    name='kitsune',
    version='0.9.2',
    description='tools for finding an optimal kmer',
    long_description=readme(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
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
    test_suite='nose.collector',
    tests_require=['nose', 'nose-cover3'],
    include_package_data=True,
    install_requires=["numpy >= 1.1.0",
                        "scipy>=0.18.1",
                        "biopython>=1.68",
                        "tqdm"
                    ],
    zip_safe=False
)
