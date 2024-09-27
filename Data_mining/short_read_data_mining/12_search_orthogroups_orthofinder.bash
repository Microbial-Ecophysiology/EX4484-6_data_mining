# Run orthofinder
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"

mkdir -p $WDIR/orthofinder_v1/input
# copy .faa sequences of all MAGs to input directory

ls -1 $WDIR/mmseq2/prodigal_out/*.faa | xargs -n 1 basename | sed 's/\.faa//' | while read line
do
 sed "/^>/s/^>/>${line}___/" $WDIR/mmseq2/prodigal_out/${line}".faa" > $WDIR/orthofinder_v1/input/${line}".faa"
done 

ls -1 $WDIR/ETH_MAGs/ETH_selection_final/prodigal/*.faa | xargs -n 1 basename | sed 's/\.faa//' | while read line
do
 sed "/^>/s/^>/>${line}___/" $WDIR/ETH_MAGs/ETH_selection_final/prodigal/${line}".faa" > $WDIR/orthofinder_v1/input/${line}".faa"
done 

conda activate orthofinder-2.5.5

cd $WDIR/orthofinder_v1
orthofinder -f $WDIR/orthofinder_v1/input/ -t 100 -a 12 -M msa

# Orthofinder analyze results
# in R group specific orthogroups for all found families were selected and saved
# now select all fasta header from group specific orthogroups, to search within annotation file for possible annotation

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"

while read SID 
do
  cut -f1 $WDIR/orthofinder_v1/${SID}"_orthogroups.txt" | sed 1d | while read line
  do
    grep -e ">" $WDIR/orthofinder_v1/input/OrthoFinder/Results_Nov10/Orthogroup_Sequences/${line}".fa"
  done >> $WDIR/orthofinder_v1/"fasta_header_orthogroups_"${SID}.txt
done < $WDIR/orthofinder_v1/Families.txt


# rerun for all Thermoplasmatota

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"

mkdir -p $WDIR/orthofinder_allThermo_v1/input
# copy .faa sequences of all MAGs to input directory

ls -1 $WDIR/mmseq2/prodigal_out/*.faa | xargs -n 1 basename | sed 's/\.faa//' | while read line
do
 sed "/^>/s/^>/>${line}___/" $WDIR/mmseq2/prodigal_out/${line}".faa" > $WDIR/orthofinder_allThermo_v1/input/${line}".faa"
done 

ls -1 $WDIR/ETH_MAGs/ETH_selection_final/prodigal/*.faa | xargs -n 1 basename | sed 's/\.faa//' | while read line
do
 sed "/^>/s/^>/>${line}___/" $WDIR/ETH_MAGs/ETH_selection_final/prodigal/${line}".faa" > $WDIR/orthofinder_allThermo_v1/input/${line}".faa"
done 

# add all MAGs from the phylogenomic tree
ls -1 $WDIR/mmseq2/prodigal_out_marker_gene_tree/*.faa | xargs -n 1 basename | sed 's/\.faa//' | while read line
do
 sed "/^>/s/^>/>${line}___/" $WDIR/mmseq2/prodigal_out_marker_gene_tree/${line}".faa" > $WDIR/orthofinder_allThermo_v1/input/${line}".faa"
done 
conda activate orthofinder-2.5.5

cd $WDIR/orthofinder_allThermo_v1
orthofinder -f $WDIR/orthofinder_allThermo_v1/input/ -t 100 -a 12 -s $WDIR/orthofinder_allThermo_v1/XaxY-UqpaqsQj0g8TZACjQ_newick.txt #rooted phylogenomic tree

