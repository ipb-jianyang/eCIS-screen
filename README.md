# Reedme for ScEDS-screen
This repository contains the code used for screening putative Self-contained Effector Delivery Systems loci from genomic sequences.

## Note
This pipeline was designed and tested on Linux system only, though it may also work on other platforms with the required Perl environment.

## Requirements
All scripts of this pipeline were custom wrtten in `Perl` (5.10.1) with `bioperl` (1.7).
Additional program required is `HMMER` (3.1b2) available from http://hmmer.org/download.html.

## Install
The script `HMMsearch_genomesII.pl` is the master of the pipeline, you need to be manually edited it to fit your local system environment before it can work properly. The two settings are given in the top of the script, includeing the path to `hmmsearch` program ($hmmsearch) and the path of this pipeline ($script_path).
After update the settings you may test the pipeline using the following command (run within the install directory):

    $ ./HMMsearch_genomesII.pl ScEDS.hmm example_genomes screen_summary.txt

## Usage
The basic usage looks like the following:

    $ HMMsearch_genomesII.pl <query HMM profile> <path to genomes for screen> <output summary file>

- The `ScEDS.hmm` file included is a pre-build HMM profile contains multiple HMMs for all components of ScEDSs.
- The path to genomes is expected to contain sub-directories for each genome, which should contain the GenBank format file named as `*_genomic.gbff.gz`. An example of the genome directory is given as `example_genomes`, which includes two genomes.
- The pipeline will produce intermdiate files (e.g. CDS proteins, HMMER output) in the CURRENT working directory (where you run the script), these files can be helpful to improve the speed of the pipeline for further runs.
- The pipeline may run for a few hours or days depends on the number of genomes for screening. The progress message will be printed to the screen by default.
- The output is in Tab-delimited text format, view in in Excel for easy read.
