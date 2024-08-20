# split_fastqcats
Containing code for splitting out concatenated fastq files


## Installation

### Conda installation - **in progress**

The preferred method for installation is through conda/mamba.  Preferably the
installation should be in a seperate environment::

    mamba env create -f conda/environments/split_fastqcats.yml
    conda activate fastcats
    pip install -e .

    split_fastqcats --help

## Usage

Run the ``split_fastqcats --help`` command view the help documentation.


Then run the pipeline::

    split_fastqcats split tests/input.fastq.gz output.fastq.gz
