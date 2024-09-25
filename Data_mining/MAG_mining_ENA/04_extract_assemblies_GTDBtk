## Run GTDB
# Thermoplasmatota assemblies were taxonomically checked using gtdbtk-2.1.0
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/Thermoplasmata"
conda activate gtdbtk-2.1.0
CPUS=60

cd $WDIR
gtdbtk classify_wf --genome_dir $WDIR/input/Thermoplasmata_checkm2_filtered/ --out_dir $WDIR/input/gtdbtk_out --cpus $CPUS -x fna  >> $WDIR/logfiles/Thermo_MAGs"_classify_gtdbtk.log" 2>&1
gtdbtk classify_wf --genome_dir $WDIR/input/literature_MAGs --out_dir $WDIR/input/literature_MAGs/gtdbtk_out --cpus $CPUS -x fna  >> $WDIR/logfiles/literature_MAGs"_classify_gtdbtk.log" 2>&1
gtdbtk classify_wf --genome_dir $WDIR/input/input_Thermoplasmatota/Thermoplasmatota_checkm2_filtered/ --out_dir $WDIR/input/input_Thermoplasmatota/gtdbtk_out --cpus $CPUS -x fna  >> $WDIR/logfiles/Thermoplasmatota_MAGs"_classify_gtdbtk.log" 2>&1


## Extract Thermoplasmatota assemblies
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/Thermoplasmata"
mkdir -p $WDIR/input/Thermoplasmata_checkm2_gtdbtk_filtered
awk -v FS="\t" -v OFS="\t" '$2 ~ /Thermoplasmatota/' $WDIR/input/gtdbtk_out/gtdbtk.ar53.summary.tsv > $WDIR/input/Thermoplasmata_checkm2_gtdbtk_filtered/Thermoplasmata_MAGs.txt
cut -f1 $WDIR/input/Thermoplasmata_checkm2_gtdbtk_filtered/Thermoplasmata_MAGs.txt | while read line
do
 cp $WDIR/input/Thermoplasmata_checkm2_filtered/${line}".fna"  $WDIR/input/Thermoplasmata_checkm2_gtdbtk_filtered/
done

mkdir -p $WDIR/input/input_Thermoplasmatota/Thermoplasmatota_gtdbtk_checkm2_filtered
awk -v FS="\t" -v OFS="\t" '$2 ~ /Thermoplasmatota/' $WDIR/input/input_Thermoplasmatota/gtdbtk_out/gtdbtk.ar53.summary.tsv > $WDIR/input/input_Thermoplasmatota/Thermoplasmatota_MAGs.txt
cut -f1 $WDIR/input/input_Thermoplasmatota/Thermoplasmatota_MAGs.txt | while read line
do
 cp $WDIR/input/input_Thermoplasmatota/Thermoplasmatota_checkm2_filtered/${line}".fna"  $WDIR/input/input_Thermoplasmatota/Thermoplasmatota_gtdbtk_checkm2_filtered
done
# 696 MAGs remain

rm $WDIR/input/input_Thermoplasmatota/*.fna
rm -r $WDIR/input/input_all #to remove all initially downlaoded MAGs
mkdir -p $WDIR/input/all_Thermo_ENA_MAGs

# move all Thermoplasmatota and Thermoplasmata MAGs in one directory via links to make sure there are no further duplicates.
# 714 unique genomes were found
cd $WDIR/input/input_Thermoplasmatota/Thermoplasmatota_gtdbtk_checkm2_filtered/
ls -1 *.fna > ../../input_Thermoplasmatota.txt
cd $WDIR/input/input_Thermoplasmata
ls -1 *.fna > ../input_Thermoplasmata.txt
cd $WDIR/input/
sort input_Thermoplasmata.txt input_Thermoplasmatota.txt | uniq > $WDIR/input/all_Thermo_ENA_MAGs/all_Thermo_ENA.txt
