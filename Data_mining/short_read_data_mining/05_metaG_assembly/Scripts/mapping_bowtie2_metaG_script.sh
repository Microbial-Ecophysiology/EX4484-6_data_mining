#!/bin/bash

SID=$1
wdir=$2

bowtie2 --threads 30 -x $wdir/mapping/contigs -1 $wdir/clean_reads/${SID}_1.fastq -2 $wdir/clean_reads/${SID}_2.fastq -S $wdir/mapping/${SID}".sam"

samtools sort -o $wdir/mapping/${SID}".bam" $wdir/mapping/${SID}".sam"

samtools index $wdir/mapping/${SID}".bam"
