## Extract only good quality MAGs from ENA assemblies

## To extract only good MAGs from ENA assemblies CheckM2 is run over all MAGs.
## From checkM2 only bins with completeness > 80% and contamination < 10% are chosen.
## A list from these bins will be created and only these MAGs will then be further used for a fastANI search against both Thermoplasmata MAGs I obtained

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/Thermoplasmata/input"
cd $WDIR/input_all

# some sequences have been downloaded multiple times, replicated were removed
rm *fna.gz.1 
rm *fna.gz.2
rm *fna.gz.3

mkdir -p $WDIR/checkM2_out $WDIR/Thermoplasmata_checkm2_filtered

# unzip all MAGs
cd $WDIR/input_all
gunzip *.fna.gz

## run checkM2

module load checkm2
conda activate checkm2

THREADS=60

cd $WDIR

checkm2 predict --threads $THREADS -x fa --input $WDIR/input_all --output-directory $WDIR/checkM2_out

## Extract good Assemblies
## as good assemblies all assemblies with COmpleteness >80 and Contamination <5 were chosen
## all filtered assemblies were copied into new directory input_filtered

awk -v FS="\t" -v OFS="\t" '$2 >= 80 && $3 <= 5' $WDIR/checkM2_out/quality_report.tsv > $WDIR/Thermoplasmata_checkm2_filtered/filtered_MAGs.txt

cut -f1 $WDIR/Thermoplasmata_checkm2_filtered/filtered_MAGs.txt | while read line
do
 cp $WDIR/input_all/${line}".fna"  $WDIR/Thermoplasmata_checkm2_filtered/
done

############################################
## Repeat also for all downloaded Thermoplasmatota MAGs from ENA
mkdir -p $WDIR/input_Thermoplasmatota/checkM2_out 

# unzip all MAGs
cd $WDIR/input_Thermoplasmatota/
gunzip *.fna.gz

## run checkM2
module load checkm2
conda activate checkm2
THREADS=60
checkm2 predict --threads $THREADS -x fna --input $WDIR/input_Thermoplasmatota/ --output-directory $WDIR/input_Thermoplasmatota/checkM2_out

## Extract good Assemblies Thermoplasmatota
## as good assemblies all assemblies with Completeness >80 and Contamination <5 were chosen
## all filtered assemblies were copied into new directory input_filtered

awk -v FS="\t" -v OFS="\t" '$2 >= 80 && $3 <= 5' $WDIR/input_Thermoplasmatota/checkM2_out/quality_report.tsv > $WDIR/input_Thermoplasmatota/filtered_MAGs.txt

cut -f1 $WDIR/input_Thermoplasmatota/filtered_MAGs.txt | while read line
do
 cp $WDIR/input_Thermoplasmatota/${line}".fna"  $WDIR/input_Thermoplasmatota/Thermoplasmatota_checkm2_filtered/
done

############################################
## also run checkm2 for MAGs from literature (Sysuiplasmatales, Gimiplasmatales, Lutacidiplasmatales) and select only good MAGs
checkm2 predict --threads 60 -x fa --input ./ --output-directory ./checkm_out

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/Thermoplasmata"

awk -v FS="\t" -v OFS="\t" '$2 >= 80 && $3 <= 5' /storage/hdd2/mmaeke/Metagenomics/PhD/Thermoplasmata_MAGs/Sysuiplasmatales_MAGs/ncbi-genomes-2022-07-25/checkm_out/quality_report.tsv > $WDIR/literature_MAGs/Sysuiplasmatales_MAGs.txt
awk -v FS="\t" -v OFS="\t" '$2 >= 80 && $3 <= 5' /storage/hdd2/mmaeke/Metagenomics/PhD/Thermoplasmata_MAGs/Lutacidiplasmatales_MAGs/ncbi-genomes-2022-07-25/checkm_out/quality_report.tsv > $WDIR/literature_MAGs/Lutacidiplasmatales_MAGs.txt
awk -v FS="\t" -v OFS="\t" '$2 >= 80 && $3 <= 5' /storage/hdd2/mmaeke/Metagenomics/PhD/Thermoplasmata_MAGs/Gimiplasmatales_MAGs_Hu_et_al_2021/checkm_out/quality_report.tsv > $WDIR/literature_MAGs/Gimiplasmatales_MAGs.txt

cut -f1 $WDIR/literature_MAGs/Sysuiplasmatales_MAGs.txt | while read line
do
 cp /storage/hdd2/mmaeke/Metagenomics/PhD/Thermoplasmata_MAGs/Sysuiplasmatales_MAGs/ncbi-genomes-2022-07-25/${line}".fna"  $WDIR/literature_MAGs/
done

cut -f1 $WDIR/literature_MAGs/Lutacidiplasmatales_MAGs.txt | while read line
do
 cp /storage/hdd2/mmaeke/Metagenomics/PhD/Thermoplasmata_MAGs/Lutacidiplasmatales_MAGs/ncbi-genomes-2022-07-25/${line}".fna"  $WDIR/literature_MAGs/
done

cut -f1 $WDIR/literature_MAGs/Gimiplasmatales_MAGs.txt | while read line
do
 cp /storage/hdd2/mmaeke/Metagenomics/PhD/Thermoplasmata_MAGs/Gimiplasmatales_MAGs_Hu_et_al_2021/${line}".fna"  $WDIR/literature_MAGs/
done
