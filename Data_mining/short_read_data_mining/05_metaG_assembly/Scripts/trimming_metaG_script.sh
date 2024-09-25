#!/bin/bash

# read input arguments
wdir=$2
SID=$1


# Adapter trimming
bbduk.sh in1=$wdir/downloaded/${SID}"_1.fastq.gz" in2=$wdir/downloaded/${SID}"_2.fastq.gz" out1=$wdir/clean_reads/${SID}"_adapter_trimmed_1.fastq.gz" out2=$wdir/clean_reads/${SID}"_adapter_trimmed_2.fastq.gz">

# Quality trimming
bbduk.sh in1=$wdir/clean_reads/${SID}"_adapter_trimmed_1.fastq.gz" in2=$wdir/clean_reads/${SID}"_adapter_trimmed_2.fastq.gz" out1=$wdir/clean_reads/${SID}"_1.fastq" out2=$wdir/clean_reads/${SID}"_2.fastq" qtr>
