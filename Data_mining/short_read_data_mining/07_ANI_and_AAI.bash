# FastANI on good quality EX MAGs

# on all MAGs including MAGs of co-authors
mkdir -p $WDIR/fastANI_v2/input
mkdir -p $WDIR/fastANI_v2/output

cd $WDIR/MAGs
ls -d "$PWD"/*.fa >> $WDIR/fastANI_v2/input/reference_list.txt
ls -d "$PWD"/*.fa >> $WDIR/fastANI_v2/input/query_list.txt

cd $WDIR/IOW_MAGs
ls -d "$PWD"/*.fasta >> $WDIR/fastANI_v2/input/reference_list.txt
ls -d "$PWD"/*.fasta >> $WDIR/fastANI_v2/input/query_list.txt

cd $WDIR/ETH_MAGs/ETH_selection_final
ls -d "$PWD"/*.fa >> $WDIR/fastANI_v2/input/reference_list.txt
ls -d "$PWD"/*.fa >> $WDIR/fastANI_v2/input/query_list.txt

conda activate gtdbtk-2.1.0
THREADS=100
fastANI --ql $WDIR/fastANI_v2/input/query_list.txt --rl $WDIR/fastANI_v2/input/reference_list.txt -o $WDIR/fastANI_v2/output/EX_fastANI_out -t $THREADS > $WDIR/logfiles/EX_MAGs"_fastANI_v2.log" 2>&1


# Run FastAAI for all MAGs
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"
mkdir -p $WDIR/fastAAI_v2
conda activate fastAAI-0.1.18
fastaai build_db -g $WDIR/fastAAI_v2/input/ -d EX4484-6_AAI.db -o $WDIR/fastAAI_v2 --threads 60
fastaai db_query -q $WDIR/fastAAI_v2/database/Thermo_AAI.db -t $WDIR/fastAAI_v2/database/Thermo_AAI.db -o $WDIR/fastAAI_v2 --output_style tsv --threads 60 

