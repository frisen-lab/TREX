![CI status](https://github.com/frisen-lab/TREX/workflows/CI/badge.svg)

Overview
--------

TREX is an experimental workflow that enables simultaneous lineage TRacking and EXpression profiling of single cells using RNA-sequencing. 

An essential part of this workflow is presented here: the extraction of genetic barcodes or "cloneIDs" from single-cell transcriptomes and the reconstruction of related cells.

The computational tool uses BAM files of one or multiple sequencing libraries as an input for the generation of cloneID count matrices and identifies clonally related cells based on Jaccard similarity between each pair of cloneID+cells. 

Currently, TREX is compatible with common RNA-sequencing library preparation methods and data formats provided by 10X Chromium, 10X Visium and Smart-Seq3. 

Installation
------------

TREX requires Python 3.6 or newer.
We recommend that you install TREX into a "virtual environment", which can be done by running the following commands:

    python3 -m venv trex-venv
    trex-venv/bin/pip install git+https://github.com/frisen-lab/TREX.git

`trex-venv` is the name of the directory that will be created and which will contain the virtual environment.
You can choose a different name.

Test the installation by running `trex-venv/bin/trex --version`.
You can also choose to run `source trex-venv/bin/activate`, which will make the command directly available, that is, you can then simply type `trex` to run it. (Examples below assume you have done this.)

Running TREX on a minimal test dataset
--------------------------------------

Clone the Git repository or [download it as a ZIP file](https://github.com/frisen-lab/TREX/archive/master.zip) and unpack it.
The directory `tests/data/` contains a test dataset in the proper format. Run TREX on it:

    trex run10x -s 695 -e 724 tests/data/

This will create a folder `trex_run/` in the current directory (use the `-o` option to choose a different folder) with the results.

The input directory for TREX must be a Cell Ranger output directory. See the contents of the `tests/data/outs` directory to learn which are the minimum files necessary. Cell Ranger must have been configured to map against a reference augmented by an extra chromosome that contains the cloneID. By default, that extra chromosome is assumed to be the last in the BAM file (use `--chromosome` to choose a different one). The options the options `-s` and `-e` set where on the extra chromosome the cloneID is located (`-s` gives start and `-e` gives end in 1-based coordinates).

Please also run `trex run10x --help` to see the other available command-line options.
