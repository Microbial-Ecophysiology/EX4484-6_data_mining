# Anvio script to refine EX MAGs found in metagenomic studies

# test with one sample: PRJNA713414

## Workflow for anvio 7.1
## This workflow is based on the assumption that you already did your bin reassembly and want to refine your target bins.

# Prepare files
# Create a concatenated assembly by combining all your reassembled bins into one file
WDIR_A="/storage/hdd1/chh/Mara_metaG_data_mining/v1/EX_MAGs"
WDIR_B="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/Thermoplasmata_MGs/metaG_MAGs"

# Change fasta header for anvio while creating a concatenated assembly file 
while read line 
do
# Change fasta header for anvio while creating a concatenated assembly file 
  mkdir -p $WDIR_B/anvio_refine/${line}/input
  ls -1 $WDIR_A/${line}/reassembled_bins/reassembled_bins | xargs -n 1 basename | while read SID
  do
	sed -e '/^>/s/\./_/' -e "/^>/s/[^[:alnum:]>-]/_/g"  $WDIR_A/${line}/reassembled_bins/reassembled_bins/${SID}
  done >> $WDIR_B/anvio_refine/${line}/input/${line}_concat_assembly.fasta

# create a list with all contigs and corresponding bins
  ls -1 $WDIR_A/${line}/reassembled_bins/reassembled_bins  | xargs -n 1 basename | sed 's/\.fa//'  | while read SID
  do
    grep '^>' $WDIR_A/${line}/reassembled_bins/reassembled_bins/${SID}".fa" | sed  -e "/^>/s/[^[:alnum:]>-]/_/g"  -e 's/^>//' -e "s/$/\t${SID}/" -e 's/\./_/g'
  done >> $WDIR_B/anvio_refine/${line}/input/${line}_binning_results.txt
 
done < $WDIR_B/sample_names.txt



# Run mapping of all clean reads on the concatenated bin assembly
conda activate anvio-7.1

while read line
do
  cat $WDIR_A/${line}/clean_reads/*1.fastq.gz > $WDIR_B/anvio_refine/${line}/input/Read_1.fastq.gz
  cat $WDIR_A/${line}/clean_reads/*2.fastq.gz > $WDIR_B/anvio_refine/${line}/input/Read_2.fastq.gz
  gunzip $WDIR_B/anvio_refine/${line}/input/*fastq.gz

  bowtie2-build $WDIR_B/anvio_refine/${line}/input/${line}_concat_assembly.fasta $WDIR_B/anvio_refine/${line}/input/contigs 
  bowtie2 --threads 100 -x $WDIR_B/anvio_refine/${line}/input/contigs  -1 $WDIR_B/anvio_refine/${line}/input/Read_1.fastq -2 $WDIR_B/anvio_refine/${line}/input/Read_2.fastq -S $WDIR_B/anvio_refine/${line}/input/${line}.sam 
  samtools sort -m 3G -o $WDIR_B/anvio_refine/${line}/input/${line}.bam $WDIR_B/anvio_refine/${line}/input/${line}.sam
  samtools index $WDIR_B/anvio_refine/${line}/input/${line}.bam

  rm $WDIR_B/anvio_refine/${line}/input/contigs*
  rm $WDIR_B/anvio_refine/${line}/input/*.sam
 
  rm $WDIR_B/anvio_refine/${line}/input/Read_1.fastq
  rm $WDIR_B/anvio_refine/${line}/input/Read_2.fastq
done < $WDIR_B/sample_names.txt


PRJNA721298
PRJNA704804
PRJNA541421
PRJNA531756

## Create a contigs DB
# While creating a contigs database k-mer frequencies are caclulated (default 4, change by using --kmer-size),
# contigs longer than 20,000 bp are split into smaller contigs (change size by using --split-length), 
# ORFs are called by using Prodigal (skip by using --skip-gene-calling)and contig splitting will cut mindfully when
# using Prodigal, so that no geens are cut in the middle.
while read line
do
  mkdir -p $WDIR_B/anvio_refine/${line}/output
  anvi-gen-contigs-database -f $WDIR_B/anvio_refine/${line}/input/${line}_concat_assembly.fasta -o $WDIR_B/anvio_refine/${line}/output/CONTIGS.db -n 'Contigs database ${line}' -T 80

## Run anvio hmm profiles
# With anvi'o there are already some hmm models supplied, constituting of multiple published bacterial single-copy gene collections
# Further you can run your own HMM profile (--hmm-profile-dir) or only a specific installed HMM profile (--installed-hmm-profile)
  anvi-run-hmms -c $WDIR_B/anvio_refine/${line}/output/CONTIGS.db --also-scan-trnas -T 80

##Run taxonomy
  anvi-run-ncbi-cogs -c $WDIR_B/anvio_refine/${line}/output/CONTIGS.db --sensitive -T 80 

# Get taxonomy
  anvi-run-scg-taxonomy -c $WDIR_B/anvio_refine/${line}/output/CONTIGS.db -T 80

## Profiling BAM files
# The anvi'o profile db contains specific information about contigs.
# Profiling a bam file creates a single profile reporting properties for each contig in a sample based on mapping.
# Each profile.db links to a contig.db

  anvi-profile -i $WDIR_B/anvio_refine/${line}/input/${line}.bam -c $WDIR_B/anvio_refine/${line}/output/CONTIGS.db --output-dir $WDIR_B/anvio_refine/${line}/profile_db --sample-name ${line} --min-contig-length 2500 --profile-SCVs -T 80

# Import own binning
# binning_results.txt should contain a TAB-delimited file containing infos about which contig belongs to what bin.
  anvi-import-collection $WDIR_B/anvio_refine/${line}/input/${line}_binning_results.txt -p $WDIR_B/anvio_refine/${line}/profile_db/PROFILE.db -c $WDIR_B/anvio_refine/${line}/output/CONTIGS.db -C Metawrap --contigs-mode

# Estimate completeness of bins
  anvi-estimate-genome-completeness -p $WDIR_B/anvio_refine/${line}/profile_db/PROFILE.db -c $WDIR_B/anvio_refine/${line}/output/CONTIGS.db -C Metawrap 

# Estimate taxonomy of bins
  anvi-estimate-scg-taxonomy -p $WDIR_B/anvio_refine/${line}/profile_db/PROFILE.db -c $WDIR_B/anvio_refine/${line}/output/CONTIGS.db -C Metawrap --compute-scg-coverages

# Estimate metabolism (2 h / 1,5 h)
  anvi-run-kegg-kofams -c $WDIR_B/anvio_refine/${line}/output/CONTIGS.db -T 80  
  anvi-estimate-metabolism -p $WDIR_B/anvio_refine/${line}/profile_db/PROFILE.db -c $WDIR_B/anvio_refine/${line}/output/CONTIGS.db -C Metawrap -O ${line} --include-metadata --metagenome-mode


 done < $WDIR_B/sample_names.txt
 
 
 # Refine MAGs
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/Thermoplasmata_MGs/metaG_MAGs"
conda activate anvio-7.1
cd $WDIR
#PRJNA713414 - bin.118, bin 178, bin.187 (not refined)
anvi-show-collections-and-bins -p $WDIR/anvio_refine/PRJNA713414/profile_db/PROFILE.db
anvi-refine -p $WDIR/anvio_refine/PRJNA713414/profile_db/PROFILE.db -c $WDIR/anvio_refine/PRJNA713414/output/CONTIGS.db -C Metawrap --server-only -P 8080 -b bin_187_orig

#PRJNA368391 - bin_183_orig (not refined), bin_90_orig
anvi-show-collections-and-bins -p $WDIR/anvio_refine/PRJNA368391/profile_db/PROFILE.db
anvi-refine -p $WDIR/anvio_refine/PRJNA368391/profile_db/PROFILE.db -c $WDIR/anvio_refine/PRJNA368391/output/CONTIGS.db -C Metawrap --server-only -P 8080 -b bin_90_orig

#PRJNA889212 - bin_2_orig (672 splits, redundancy <40), bin_5_orig (classified by anvio as Bacteria), bin_7_strict and bin_1_orig have been classified as EX4484-6, both were refined
anvi-show-collections-and-bins -p $WDIR/anvio_refine/PRJNA889212/profile_db/PROFILE.db
anvi-refine -p $WDIR/anvio_refine/PRJNA889212/profile_db/PROFILE.db -c $WDIR/anvio_refine/PRJNA889212/output/CONTIGS.db -C Metawrap --server-only -P 8080 -b bin_5_orig

#PRJNA704804 - bin_147_orig
anvi-show-collections-and-bins -p $WDIR/anvio_refine/PRJNA704804/profile_db/PROFILE.db
anvi-refine -p $WDIR/anvio_refine/PRJNA704804/profile_db/PROFILE.db -c $WDIR/anvio_refine/PRJNA704804/output/CONTIGS.db -C Metawrap --server-only -P 8080 -b bin_147_orig

#PRJNA541421 - bin_70_strict (not refined,redundancy =0), bin_72_orig, bin_98_orig
anvi-show-collections-and-bins -p $WDIR/anvio_refine/PRJNA541421/profile_db/PROFILE.db
anvi-refine -p $WDIR/anvio_refine/PRJNA541421/profile_db/PROFILE.db -c $WDIR/anvio_refine/PRJNA541421/output/CONTIGS.db -C Metawrap --server-only -P 8080 -b bin_98_orig

#PRJNA531756 - bin_230_orig, bin_57_strict (completeness <60%), bin_72_orig, bin_93_strict
anvi-show-collections-and-bins -p $WDIR/anvio_refine/PRJNA531756/profile_db/PROFILE.db
anvi-refine -p $WDIR/anvio_refine/PRJNA531756/profile_db/PROFILE.db -c $WDIR/anvio_refine/PRJNA531756/output/CONTIGS.db -C Metawrap --server-only -P 8080 -b bin_93_strict

#PRJNA721298 - bin_46_orig
anvi-show-collections-and-bins -p $WDIR/anvio_refine/PRJNA721298/profile_db/PROFILE.db
anvi-refine -p $WDIR/anvio_refine/PRJNA721298/profile_db/PROFILE.db -c $WDIR/anvio_refine/PRJNA721298/output/CONTIGS.db -C Metawrap --server-only -P 8080 -b bin_46_orig

while read line
do
   anvi-summarize -p $WDIR_B/anvio_refine/${line}/profile_db/PROFILE.db -c $WDIR_B/anvio_refine/${line}/output/CONTIGS.db -o $WDIR_B/anvio_refine/${line}/summary -C Metawrap
done < $WDIR_B/sample_names.txt
# check all MAGs before dereplication
