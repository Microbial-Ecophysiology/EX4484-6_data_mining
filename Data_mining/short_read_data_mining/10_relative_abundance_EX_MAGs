# Run dereplication of all EX4484-6 MAGs for final mapping 
## Running dRep for bin dereplication
## create links from all EX MAGs to $WDIR/dREP/input/MAGs)


WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"
mkdir -p $WDIR/dREP/input/MAGs
cd $WDIR/dREP/input/MAGs

ln -s  $WDIR/ETH_MAGs/ETH_selection_final/MAGs/*.fa ./
ln -s  $WDIR/MAGs/*.fa ./
ln -s  $WDIR/IOW_MAGs/*.fa ./

#create list of bins
ls -d "$PWD"/* > ../bin_list.txt 

# create input file for quality data
cut -d$'\t' -f 1-3 $WDIR/checkM2_out/quality_report.tsv | sed -e '1s/Name/genome/' -e '1s/Completeness/completeness/' -e '1s/Contamination/contamination/' | awk  'NR>=2 {$1=$1".fa"}1' >> $WDIR/dREP/input/checkm2_quality2.csv

# add IOW MAG stats
cut -d$'\t' -f 1-3 $WDIR/IOW_MAGs/checkM2_out/quality_report.tsv | awk  '{$1=$1".fa"}1' | sed '1d' >> $WDIR/dREP/input/checkm2_quality2.csv

# add ETH MAGs
ls -1 $WDIR/ETH_MAGs/ETH_selection_final/MAGs/*.fa | xargs -n1 basename | sed 's/\.fa//' | while read line
do
  grep ${line} $WDIR/ETH_MAGs/checkM2_out/quality_report.tsv | cut -d$'\t' -f 1-3 | awk  '{$1=$1".fa"}1'
done >> $WDIR/dREP/input/checkm2_quality2.csv

cat $WDIR/dREP/input/checkm2_quality2.csv  |  sed -e 's/\t/,/g' -e 's/ /,/g' > $WDIR/dREP/input/checkm2_quality_corrected.csv


# For relative abundance data MAGs
THREADS=60
MINLENGTH=50000
COMPLETENESS=80
CONTAMINATION=5
PANI=0.90
SANI=0.95
MINOVERLAP=0.5

conda activate drep-3.0.0
dRep dereplicate -p $THREADS -g $WDIR/dREP/input/bin_list.txt --length $MINLENGTH --genomeInfo $WDIR/dREP/input/checkm2_quality_corrected.csv --S_algorithm ANImf -pa $PANI -sa $SANI -nc $MINOVERLAP $WDIR/dREP/output > $WDIR/logfiles/EX_MAGs_dRep_dereplication.log 2>&1


## index for mapping against 8287
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"
ls -1 $WDIR/dREP/output/dereplicated_genomes/*.fa | xargs -n 1 basename | sed 's/\.fa//' | while read line
do
 sed "/^>/s/^>/>${line}___/" $WDIR/dREP/output/dereplicated_genomes/${line}".fa"
done >> $WDIR/dREP/index/concat_MAGs_index.fa

mkdir -p $WDIR/dREP/index/EX_MAGs_mapping
# create directory with all MAGs for index
ls -1 $WDIR/dREP/output/dereplicated_genomes/*.fa | xargs -n 1 basename | sed 's/\.fa//' | while read line
do
 sed "/^>/s/^>/>${line}___/" $WDIR/dREP/output/dereplicated_genomes/${line}".fa" > $WDIR/dREP/index/EX_MAGs_mapping/${line}".fa"
done 

conda activate-anvio7.1
bowtie2-build $WDIR/dREP/index/concat_MAGs_index.fa $WDIR/dREP/index/EX4484-6 >> $WDIR/logfiles/"EX4484-6_final_indexing.log" 2>&1
# bowtie 2.3.5
