Search for selection in MPRA data
====

Smith, McManus, and Fraser 2013 and Agoglia and Fraser 2013 both used data from
a Massively Parallel Reporter Assay (MPRA) to test for selection on enhancers
and non-coding elements. In this project, I'm hoping to find similar signs of
selection in the saturation mutagenesis data from the Ahithuv lab [1].



[1] https://mpra.gs.washington.edu/satMutMPRA/ and https://www.biorxiv.org/content/10.1101/505362v1


## Prerequisites

The goal of the project is to have a basically turnkey solution. To that end,
as many of the dependencies as possible are loaded into a conda environment.
However, you will still need a copy of the RepBase database.  Instructions on
how to load that in to come, but for now, I expect it to be in
Reference/RepBase24.03.fasta/all.fasta

    cat *.ref appendix/*.ref > all.fasta


You will also need to have installed:

    * FigTree (sudo apt-get install figtree)
    * Snakemake and Python
    * The RefSeq Genomic blast database
    * Probably other things I haven't discovered yet (please file an issue!)

