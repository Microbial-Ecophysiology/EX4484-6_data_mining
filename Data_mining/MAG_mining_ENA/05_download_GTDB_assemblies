## Check ENA assemblies against GTDB genomes
# Download Thermoplasmata metadata file from GTDB
wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/ar53_metadata_r207.tar.gz

#Use checkM output as quality filter
awk -v FS="\t" -v OFS="\t" '$3 >= 80 && $4 <= 5' $WDIR/input/gtdb_genomes/ar53_metadata_r207.tsv > $WDIR/input/gtdb_genomes/arc_gtdb_filtered.txt
#Filter out all Thermoplasmatota MAGs
awk -v FS="\t" -v OFS="\t" '$17 ~ /Thermoplasmatota/' $WDIR/input/gtdb_genomes/arc_gtdb_filtered.txt > $WDIR/input/gtdb_genomes/arc_gtdb_Thermo.txt

#Search for all GTDB MAGs within ENA Thermoplasmata MAGs
cut -f1 $WDIR/input/gtdb_genomes/arc_gtdb_Thermo.txt | cut -d '_' -f3 | sed -r 's/\.[0-9]+/./' | sed 's/^/_/' | grep -c -F -f - <(cut -f1 $WDIR/input/all_Thermo_ENA_MAGs/all_Thermo_ENA.txt)
# 348 / 535 found

# check how many of these sequences are available in refseq
cut -f1 $WDIR/input/all_Thermo_ENA_MAGs/all_Thermo_ENA.txt | cut -d '_' -f2 | sed -r 's/\.[0-9]+/./' | sed 's/^/_/' | grep -v -F -f - <(cut -f1 $WDIR/input/gtdb_genomes/arc_gtdb_Thermo.txt) | cut -d'_' -f3 | sed -r 's/\.[0-9]+/./' | sed 's/^/_/' | grep -F -f - $WDIR/assembly_summary_refseq.txt | wc -l  #191

# check how many of these sequences are available in genbank
cut -f1 $WDIR/input/all_Thermo_ENA_MAGs/all_Thermo_ENA.txt | cut -d '_' -f2 | sed -r 's/\.[0-9]+/./' | sed 's/^/_/' | grep -v -F -f - <(cut -f1 $WDIR/input/gtdb_genomes/arc_gtdb_Thermo.txt) | cut -d'_' -f3 | sed -r 's/\.[0-9]+/./' | sed 's/^/_/' | grep -F -f - $WDIR/assembly_summary_genbank.txt | wc -l # 191

# get a list of all refseq MAGs
cut -f1 $WDIR/input/all_Thermo_ENA_MAGs/all_Thermo_ENA.txt | cut -d '_' -f2 | sed -r 's/\.[0-9]+/./' | sed 's/^/_/' | grep -v -F -f - <(cut -f1 $WDIR/input/gtdb_genomes/arc_gtdb_Thermo.txt) | cut -d'_' -f3 | sed -r 's/\.[0-9]+/./' | sed 's/^/_/' | grep -F -f - $WDIR/assembly_summary_refseq.txt | cut -f20 > $WDIR/input/gtdb_genomes/Links_refseq.txt

# Get 12 refseq assemblies
while read line
do
  ACC=$(basename ${line})
  LINK=$(echo $line | sed "s/$/\/${ACC}_genomic.fna.gz/")
  wget $LINK
done < Links_refseq.txt


## reverse grep to find missing 290 MAGs in genbank file

# create a list of all found refseq MAGs
cut -f1 $WDIR/input/all_Thermo_ENA_MAGs/all_Thermo_ENA.txt | cut -d '_' -f2 | sed -r 's/\.[0-9]+/./' | sed 's/^/_/' | grep -v -F -f - <(cut -f1 $WDIR/input/gtdb_genomes/arc_gtdb_Thermo.txt) | cut -d'_' -f3 | sed -r 's/\.[0-9]+/./' | sed 's/^/_/' | grep -F -f - $WDIR/assembly_summary_refseq.txt > $WDIR/input/gtdb_genomes/refseq_MAGs_all.txt

#create a list of all found genbank MAGs
cut -f1 $WDIR/input/all_Thermo_ENA_MAGs/all_Thermo_ENA.txt | cut -d '_' -f2 | sed -r 's/\.[0-9]+/./' | sed 's/^/_/' | grep -v -F -f - <(cut -f1 $WDIR/input/gtdb_genomes/arc_gtdb_Thermo.txt) | cut -d'_' -f3 | sed -r 's/\.[0-9]+/./' | sed 's/^/_/' | grep -F -f - $WDIR/assembly_summary_genbank.txt > $WDIR/input/gtdb_genomes/genbank_MAGs_all.txt

# reverse search to find all remaining genbank links
cut -f1 $WDIR/input/gtdb_genomes/refseq_MAGs_all.txt | cut -d '_' -f2 | sed -r 's/\.[0-9]+/./' | sed 's/^/_/' | grep -v -F -f - $WDIR/input/gtdb_genomes/genbank_MAGs_all.txt | cut -f20 > Links_genbank.txt


# Get 191 genbank assemblies
while read line
do
  ACC=$(basename ${line})
  LINK=$(echo $line | sed "s/$/\/${ACC}_genomic.fna.gz/")
  wget $LINK
done < Links_genbank.txt


## run checkM2 on GTDB MAGs
module load checkm2
conda activate checkm2
THREADS=60
cd $WDIR

checkm2 predict --threads $THREADS -x fna --input $WDIR/input/gtdb_genomes/ --output-directory $WDIR/input/gtdb_genomes/checkM2_out

awk -v FS="\t" -v OFS="\t" '$2 >= 80 && $3 <= 5' $WDIR/input/gtdb_genomes/checkM2_out/quality_report.tsv  > $WDIR/input/gtdb_genomes/filtered_checkm2_Thermo.txt

mkdir -p $WDIR/input/gtdb_genomes/filtered
cut -f1 $WDIR/input/gtdb_genomes/filtered_checkm2_Thermo.txt | while read line
do
 cp $WDIR/input/gtdb_genomes/${line}".fna"  $WDIR/input/gtdb_genomes/filtered
done
# 132 genomes remain
