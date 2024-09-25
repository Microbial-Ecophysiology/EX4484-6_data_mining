## Download MAGs from NCBI for further analyses

# Create .txt files with all links of found MAGs

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/Thermoplasmata/"
cd $WDIR/input

cat $WDIR/subset_assemblies_scientific_name.txt | sed 1d | cut -f1 | sed s/\"//g | while read SID
do
  grep $SID $WDIR/assembly_summary_refseq.txt | cut -f20 >> input_all/Links.txt
done 

cat $WDIR/subset_assemblies_scientific_name.txt | sed 1d | cut -f1 | sed s/\"//g | while read SID
do
  grep $SID $WDIR/assembly_summary_genbank.txt | cut -f20 >> input_all/Links_genbank.txt
done 

wc -l Links_genbank.txt # 12088 links // orig file 12091, all Refseq genes were also present in Genbank, therefore only Genbank list was used

# Download all MAGs

while read line
do
  ACC=$(basename ${line})
  LINK=$(echo $line | sed "s/$/\/${ACC}_genomic.fna.gz/")
  wget $LINK
done < Links_genbank.txt

############################################
## Download also all Thermoplasmatota MAGs

cat $WDIR/input/input_Thermoplasmatota/results_assembly_Thermoplasmatota_tsv.txt | sed 1d | cut -f1 | sed s/\"//g | while read SID
do
  grep $SID $WDIR/assembly_summary_refseq.txt | cut -f20 >> $WDIR/input/input_Thermoplasmatota/Links_refseq.txt
done 

# create a list of all found refseq MAGs
cut -f1 $WDIR/input/input_Thermoplasmatota/results_assembly_Thermoplasmatota_tsv.txt | cut -d '_' -f2 | sed -r 's/\.[0-9]+/./' | sed 's/^/_/' | grep -F -f - $WDIR/assembly_summary_refseq.txt > $WDIR/input/input_Thermoplasmatota/refseq_MAGs_all.txt

#create a list of all found genbank MAGs
cut -f1 $WDIR/input/input_Thermoplasmatota/results_assembly_Thermoplasmatota_tsv.txt | cut -d '_' -f2 | sed -r 's/\.[0-9]+/./' | sed 's/^/_/' | grep -F -f - $WDIR/assembly_summary_genbank.txt > $WDIR/input/input_Thermoplasmatota/genbank_MAGs_all.txt

# reverse search to find all remaining genbank links
cut -f1 $WDIR/input/input_Thermoplasmatota/refseq_MAGs_all.txt | cut -d '_' -f2 | sed -r 's/\.[0-9]+/./' | sed 's/^/_/' | grep -v -F -f - $WDIR/input/input_Thermoplasmatota/genbank_MAGs_all.txt | cut -f20 > Links_genbank.txt

while read line
do
  ACC=$(basename ${line})
  LINK=$(echo $line | sed "s/$/\/${ACC}_genomic.fna.gz/")
  wget $LINK
done < Links_refseq.txt

while read line
do
  ACC=$(basename ${line})
  LINK=$(echo $line | sed "s/$/\/${ACC}_genomic.fna.gz/")
  wget $LINK
done < Links_genbank.txt
