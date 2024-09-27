# Check abundance of bins using coverM

conda activate coverm-0.6.1

#EX4484-6 
mkdir -p /storage/hdd1/chh/Mara_metaG_data_mining/v1/coverm/EX/bin_count/ /storage/hdd1/chh/Mara_metaG_data_mining/v1/coverm/EX/mean_cov/ /storage/hdd1/chh/Mara_metaG_data_mining/v1/coverm/EX/mincovfrac10/
ls -1 /storage/hdd1/chh/Mara_metaG_data_mining/v1/mapping_EX/*filt.bam | xargs -n1 basename | sed "s/\..*//" | while read line
do
  coverm genome -b /storage/hdd1/chh/Mara_metaG_data_mining/v1/mapping_EX/${line}".bam" -m count --min-read-percent-identity 50 --min-read-aligned-length 30 --exclude-supplementary --min-covered-fraction 0 --genome-fasta-files /storage/hdd1/chh/Mara_metaG_data_mining/v1/ref_EX_genomes/*.fa -t 50 > /storage/hdd1/chh/Mara_metaG_data_mining/v1/coverm/EX/bin_count/EX4484-6_bin_count_coverm_${line}.txt

 coverm genome -b /storage/hdd1/chh/Mara_metaG_data_mining/v1/mapping_EX/${line}".bam" -m mean --min-read-percent-identity 50 --min-read-aligned-length 30 --exclude-supplementary --min-covered-fraction 0 --genome-fasta-files /storage/hdd1/chh/Mara_metaG_data_mining/v1/ref_EX_genomes/*.fa -t 50 > /storage/hdd1/chh/Mara_metaG_data_mining/v1/coverm/EX/mean_cov/EX4484-6_bin_mean_cov_coverm_${line}.txt

 coverm genome -b /storage/hdd1/chh/Mara_metaG_data_mining/v1/mapping_EX/${line}".bam" -m mean --min-read-percent-identity 50 --min-read-aligned-length 30 --exclude-supplementary --genome-fasta-files /storage/hdd1/chh/Mara_metaG_data_mining/v1/ref_EX_genomes/*.fa -t 50 > /storage/hdd1/chh/Mara_metaG_data_mining/v1/coverm/EX/mincovfrac10/EX4484-6_bin_mean_cov_coverm_mincovfrac10_${line}.txt
done 


# the flag --min-read-aligned-percent from 30 upwards was checked to determine for our best cutoff. We chose 50 as cutoff and will now run this to select all Metagenomes that show presence of our MAGs
mkdir -p /storage/hdd1/chh/Mara_metaG_data_mining/v1/coverm/EX/bin_count_aln_perc_50/ 
ls -1 /storage/hdd1/chh/Mara_metaG_data_mining/v1/mapping_EX/*filt.bam | xargs -n1 basename | sed "s/\..*//" | while read line
do
  coverm genome -b /storage/hdd1/chh/Mara_metaG_data_mining/v1/mapping_EX/${line}".bam" -m count --min-read-percent-identity 50 --min-read-aligned-percent 50 --min-read-aligned-length 30 --exclude-supplementary --min-covered-fraction 0 --genome-fasta-files /storage/hdd1/chh/Mara_metaG_data_mining/v1/ref_EX_genomes/*.fa -t 100 > /storage/hdd1/chh/Mara_metaG_data_mining/v1/coverm/EX/bin_count_aln_perc_50/EX4484-6_bin_count_coverm_${line}.txt
done 

# same for 95% similarity
mkdir -p /storage/hdd1/chh/Mara_metaG_data_mining/v1/coverm/EX/bin_count_sim_95_aln_perc_50/ 
ls -1 /storage/hdd1/chh/Mara_metaG_data_mining/v1/mapping_EX/*filt.bam | xargs -n1 basename | sed "s/\..*//" | while read line
do
  coverm genome -b /storage/hdd1/chh/Mara_metaG_data_mining/v1/mapping_EX/${line}".bam" -m count --min-read-percent-identity 95 --min-read-aligned-percent 50 --min-read-aligned-length 30 --exclude-supplementary --min-covered-fraction 0 --genome-fasta-files /storage/hdd1/chh/Mara_metaG_data_mining/v1/ref_EX_genomes/*.fa -t 100 > /storage/hdd1/chh/Mara_metaG_data_mining/v1/coverm/EX/bin_count_sim_95_aln_perc_50/EX4484-6_bin_count_coverm_${line}.txt
done 
 
