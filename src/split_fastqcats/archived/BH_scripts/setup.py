import sys
import os
from setuptools import setup, find_packages

# Check Python version
if sys.version_info < (3, 6):
    raise SystemExit("""Requires Python 3.6 or later.""")

# Read requirements
with open('requirements.txt') as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]

# Read README for long description
with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

# Read version from __init__.py
with open('src/__init__.py', 'r') as f:
    for line in f:
        if line.startswith('__version__'):
            version = line.strip().split('=')[1].strip(' \'"')
            break

setup(
    name='split_fastqcats',
    version=version,
    description='A tool for processing and splitting FastQ reads based on primer sequences',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Adam Cribbs',
    author_email='adam.cribbs@ndorms.ox.ac.uk',
    url='https://github.com/cribbslab/split_fastqcats',
    license='MIT',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'split-fastqcats=split_fastqcats.fastq_splitter:main',
        ],
    },
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS',
    ],
    keywords=['bioinformatics', 'sequencing', 'fastq', 'primer', 'processing'],
    include_package_data=True,
    zip_safe=False,
    test_suite='tests',
)
