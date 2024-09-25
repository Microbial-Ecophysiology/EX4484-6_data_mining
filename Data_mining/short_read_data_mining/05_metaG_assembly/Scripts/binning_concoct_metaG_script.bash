#!/bin/bash

wdir=$1

# CONCOCT Binning
cut_up_fasta.py $wdir/assembly/final.contigs.fa -c 10000 -o 0 --merge_last -b $wdir/binning/concoct/contigs_10K.bed > $wdir/binning/concoct/contigs_10K.fa

concoct_coverage_table.py $wdir/binning/concoct/contigs_10K.bed $wdir/mapping/*.bam > $wdir/binning/concoct/coverage_table.tsv

concoct --composition_file $wdir/binning/concoct/contigs_10K.fa --coverage_file $wdir/binning/concoct/coverage_table.tsv -b $wdir/binning/concoct/concoct_output/ -l 2500 -t 80

merge_cutup_clustering.py $wdir/binning/concoct/concoct_output/clustering_gt2500.csv > $wdir/binning/concoct/concoct_output/clustering_merged.csv

extract_fasta_bins.py $wdir/assembly/final.contigs.fa $wdir/binning/concoct/concoct_output/clustering_merged.csv --output_path $wdir/binning/concoct/concoct_output/fasta_bins
