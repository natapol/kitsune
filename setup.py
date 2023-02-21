import sys

import setuptools

from kitsune.kitsune import __version__

if sys.version_info[0] < 3 or (sys.version_info[0] == 3 and sys.version_info[1] < 5):
    sys.stdout.write(
        "kitsune requires Python 3.5 or higher. Your current Python version is {}.{}.{}\n".format(
            sys.version_info[0], sys.version_info[1], sys.version_info[2]
        )
    )

setuptools.setup(
    author='Natapol Pornputtapong',
    author_email='natapol.por@gmail.com',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
    ],
    description=(
        "A toolkit for evaluation of the lenght of k-mer in a given genome dataset for "
        "alignment-free phylogenimic analysis"
    ),
    download_url='https://pypi.org/project/kitsune/',
    entry_points={
        'console_scripts': ['kitsune=kitsune.kitsune:main'],
    },
    install_requires=[
        "numpy>=1.1.0",
        "scipy>=0.18.1",
        "tqdm>=4.32"
    ],
    keywords=[
        "bioinformatics",
        "kitsune"
    ],
    license="GPL-3.0 License",
    license_files=["LICENSE"],
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    name='kitsune',
    packages=setuptools.find_packages(),
    platform=["Linux", "Mac OSX"],
    python_requires=">=3.5",
    url='https://github.com/natapol/kitsune',
    version=__version__,    
    zip_safe=False
)
