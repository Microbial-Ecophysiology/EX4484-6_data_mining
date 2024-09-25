#!/bin/bash

wdir=$1

#Metabat2 Binning

jgi_summarize_bam_contig_depths --outputDepth $wdir/mapping/depth.txt $wdir/mapping/*.bam
metabat2 -i $wdir/assembly/final.contigs.fa -a $wdir/mapping/depth.txt -o $wdir/binning/metabat/PRJNA889212 -m 2500 -v -t 40
