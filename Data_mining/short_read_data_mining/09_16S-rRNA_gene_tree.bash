# Prepare ARB tree for 16S rRNA gene phylogenetic tree
#ARB
# download 11764 sequences of Thermoplasmatota from SILVA
# sequence length > 1300
# pintail quality >= 30
# sequence quality >50
# alignment quality >50
# allow merging duplicates
# 10769 new species transferred and added to tree via quick add marked

# For new 16S sequences alignment via SINA (SINA (v1.2.11))
# reject sequences below identity of 70%
# Importing sequences into ARB using fasta w gaps.
# Go to align and check on all sequences. Select all full length 16S sequences after manual alignment --> 271 sequences remain
# Quick add marked to existing tree: pos_var_ssuref:archaea + termini , only columns 123456.0 (valid columns 1370)
# export with vertical gaps, filters as above, 351 sequences, incl 10 outgroup were selected, position 1109 - 41790


## Run Barrnap on all ENA data mined Thermoplasmatota assemblies to retrieve 16S rRNA genes
conda activate barrnap-0.9
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/Thermoplasmata"

THREADS=60

#run before adding Thermoplasmatota MAGs
ls -1 $WDIR/dRep/input/MAGs/*.fna | sed 's/\.fna//' | while read line
do
 barrnap --threads $THREADS --kingdom arc --outseq ${line}.gff --quiet < ${line}".fna"
done

#run after also adding Thermpoplasmatota MAGs
ls -1 $WDIR/dRep/input/MAGs/*.fa | sed 's/\.fa//' | while read line
do
 barrnap --threads $THREADS --kingdom arc --outseq ${line}.gff --quiet < ${line}".fa"
done

# copy files that don't yet exist from first run (total 168 files)
ls -1 $WDIR/16S_genes/*.gff | xargs -n1 basename | grep -v -F -f - <(ls -1 $WDIR/dRep/input/MAGs/*.gff) | xargs -n1 basename | while read line
do
 cp $WDIR/dRep/input/MAGs/$line $WDIR/16S_genes/new_16S/$line
done

# this gives me all present rRNA genes including 5S and 23S
# now I need to get only 16S genes

mv dRep/input/MAGs/*.gff 16S_genes/
ls -1 $WDIR/16S_genes/new_16S/*.gff | xargs -n1 basename | sed  's/\.gff//' | while read line
do
 grep -A1 '^>16S' $WDIR/16S_genes/new_16S/${line}.gff | sed "/^>/s/^>/>${line}_/"   >> $WDIR/16S_genes/new_16S/concat_16S_new_MAGs.gff
done
# from all dereplicated genomes 431 16S rRNA genes were retrieved for a 16S rRNA gene phylogenetic tree


# extract 16S sequences from new MAGs retrieved through short read data mining
conda activate barrnap-0.9
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"

THREADS=60
ls -1 $WDIR/IOW_MAGs/*fasta | sed 's/\.fasta//' | while read line
do
 barrnap --threads $THREADS --kingdom arc --outseq ${line}.gff --quiet < ${line}".fasta"
done

# this gives me all present rRNA genes including 5S and 23S
# now I need to get only 16S genes

mv $WDIR/IOW_MAGs/*.gff $WDIR/tree/16S/16S_genes/
ls -1 $WDIR/tree/16S/16S_genes/*.gff | xargs -n1 basename | sed  's/\.gff//' | while read line
do
 grep -A1 '^>16S' $WDIR/tree/16S/16S_genes/${line}.gff | sed "/^>/s/^>/>${line}_/"   >> $WDIR/tree/16S/16S_genes/concat_16S_MAGs.gff
done

#IOW MAGS will be added with shorter MAGs and ASVs

# in total 10 16S sequences were found
# proceed by aligning sequences with SINA, then add sequences to the existing ARB tree containing Thermoplasmatota and 16S sequences from MAG mining. 
# Remove sequences which are not full length (add later) - PRJNA541421 bin 70, PRJNA531756 bin 93, PRJNA541421 bin 72
# alignment (E.coli base 52 (1109) - 1391 (41790). In total 354 sequences remain (file: MSA_seqs_plus_metaG_mining)

#run modeltest with 16S alignment
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"

conda activate modeltest-ng-0.1.7
mkdir -p $WDIR/tree/16S/modeltest
THREADS=100
modeltest-ng -i $WDIR/tree/16S/alignment/MSA_seqs_plus_metaG_mining2.fasta -o $WDIR/tree/16S/modeltest/16S_Thermo_1109_41790 -d nt -p $THREADS -t ml
# chosen model GTR+I+G4

# raxml tree
conda activate raxml-ng-1.1.0
mkdir -p $WDIR/tree/16S/raxml
cd $WDIR/tree/16S/raxml

# check alignment
raxml-ng --check --msa $WDIR/tree/16S/alignment/MSA_seqs_plus_metaG_mining2.fasta --model GTR+I+G4 --prefix T1
#check failed therefore changing header 
cat $WDIR/tree/16S/alignment/MSA_seqs_plus_metaG_mining2.fasta | sed "/^>/s/[^[:alnum:]>-]/_/g" > $WDIR/tree/16S/alignment/MSA_seqs_plus_metaG_mining2_fixed_header.fasta
raxml-ng --check  --msa $WDIR/tree/16S/alignment/MSA_seqs_plus_metaG_mining2_fixed_header.fasta --model GTR+I+G4 --prefix T2
#new test gives a passed check

#infer default starting trees
raxml-ng --msa $WDIR/tree/16S/alignment/MSA_seqs_plus_metaG_mining2_fixed_header.fasta --msa-format fasta --model GTR+I+G4 --prefix T3 --threads 8 --seed 2 
# infer 50 starting trees to check if results become better
raxml-ng --msa $WDIR/tree/16S/alignment/MSA_seqs_plus_metaG_mining2_fixed_header.fasta --msa-format fasta --model GTR+I+G4 --tree pars{25},rand{25} --prefix T4 --threads 12 --seed 2 
# infer 70 starting trees 
raxml-ng --msa $WDIR/tree/16S/alignment/MSA_seqs_plus_metaG_mining2_fixed_header.fasta --msa-format fasta --model GTR+I+G4 --tree pars{35},rand{35} --prefix T5 --threads 12 --seed 2 

grep "Final LogLikelihood:" T{3,4,5}.raxml.log
# T4 and T5 are the same likelihoods, therefore going with T4 (50 starting trees), however rf distances are different and there are 50 distinct topologies

# Run bootstraps
raxml-ng --bootstrap --msa $WDIR/tree/16S/alignment/MSA_seqs_plus_metaG_mining2_fixed_header.fasta --msa-format fasta --model GTR+I+G4 --prefix T6 --seed 2 --threads 12 --bs-trees 1500

raxml-ng --bsconverge --bs-trees T6.raxml.bootstraps --prefix T7 --seed 2 --threads 12 --bs-cutoff 0.03
# convergence after 1300 trees

#Calculating branch support using the best Raxml tree from T4
raxml-ng --support --tree T4.raxml.bestTree --bs-trees T6.raxml.bootstraps --prefix T8 --threads 12
raxml-ng --support --tree T4.raxml.bestTree --bs-trees T6.raxml.bootstraps --prefix T9 --threads 12 --bs-metric tbe

################# Adding more sequences to tree ###################
# After having the 16S tree more shorter sequences will be added to the phylogenetic tree
# Sequences were obtained through 16S rRNA Amplicon sequencing and SIP on protein samples, additionally shorter 16S sequences will be added

# Using tree T8.raxml.support (/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs/tree/16S/raxml/T8.raxml.support)
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"

# aligning ASV sequences with mothur
module load mothur/1.45.3
#copy 16S MSalignment file to directory with new species to be aligned
cp $WDIR/tree/16S/alignment/MSA_seqs_plus_metaG_mining2_fixed_header.fasta $WDIR/tree/16S/ASVs/

#filter alignment before so that alignment lengths are equal, successful worked
cd $WDIR/tree/16S/ASVs
mothur > filter.seqs(fasta=MSA_seqs_plus_metaG_mining2_fixed_header.fasta)

#Length of filtered alignment: 1303
#Number of columns removed: 0
#Length of the original alignment: 1303
#Number of sequences used to construct filter: 354

mothur
align.seqs(candidate=ASVs_MAG16S_short.fasta, template=MSA_seqs_plus_metaG_mining2_fixed_header.fasta)

#[WARNING]: 4 of your sequences generated alignments that eliminated too many bases, a list is provided in ASVs_MAG16S_short.flip.accnos.
#OTU9919790217396_methanomass_protein
#OTU22040881119034_methanomass_protein
#OTU84494960656293_methanomass_protein
#OTU83444820243530_uncultured_Thermo

#add species to existing tree
conda activate Pbdas-v0.3.8
mkdir -p $WDIR/tree/16S/epa_placements/
epa-ng --ref-msa $WDIR/tree/16S/ASVs/MSA_seqs_plus_metaG_mining2_fixed_header.filter.fasta --tree $WDIR/tree/16S/raxml/T8.raxml.support --query $WDIR/tree/16S/ASVs/ASVs_MAG16S_short.align --model GTR+I+G4 --outdir $WDIR/tree/16S/epa_placements

# put new branches in existing tree using gappa
mkdir -p $WDIR/tree/16S/epa_placements/tree_epa_concat
gappa examine graft --jplace-path $WDIR/tree/16S/epa_placements/epa_result.jplace --out-dir $WDIR/tree/16S/epa_placements/tree_epa_concat --file-prefix 16S_tree_plusASVs_ --threads 60 --log-file $WDIR/tree/16S/epa_placements/tree_epa_concat/logfile


### add all ETH 16S genes to existing 16S tree

# start with extracting 16S sequences from new MAGs
conda activate barrnap-0.9
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"

THREADS=60

ls -1 $WDIR/ETH_MAGs/ETH_selection_final/MAGs/*fa | sed 's/\.fa//' | while read line
do
 barrnap --threads $THREADS --kingdom arc --outseq ${line}.gff --quiet < ${line}".fa"
done


ls -1 $WDIR/ETH_MAGs/ETH_selection_final/16S_genes/*.gff | xargs -n1 basename | sed  's/\.gff//' | while read line
do
 grep -A1 '^>16S' $WDIR/ETH_MAGs/ETH_selection_final/16S_genes/${line}.gff | sed "/^>/s/^>/>${line}_/"   >> $WDIR/ETH_MAGs/ETH_selection_final/16S_genes/concat_16S_MAGs.gff
done

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"

# aligning ASV sequences with mothur
module load mothur/1.45.3
mothur
align.seqs(candidate=concat_16S_ETH_MAGs.gff, template=MSA_seqs_plus_metaG_mining2_fixed_header.fasta)

# combine ASVs_MAG16S_short.align and concat_16S_ETH_MAGs.align
cat ASVs_MAG16S_short.align concat_16S_ETH_MAGs.align > 16S_genes_all.align

#add species to existing tree
conda activate Pbdas-v0.3.8
mkdir -p $WDIR/tree/16S/epa_placements_v2/
epa-ng --ref-msa $WDIR/tree/16S/ASVs/MSA_seqs_plus_metaG_mining2_fixed_header.filter.fasta --tree $WDIR/tree/16S/raxml/T8.raxml.support --query $WDIR/tree/16S/ASVs/16S_genes_all.align --model GTR+I+G4 --outdir $WDIR/tree/16S/epa_placements_v2

# put new branches in existing tree using gappa
mkdir -p $WDIR/tree/16S/epa_placements_v2/tree_epa_concat
gappa examine graft --jplace-path $WDIR/tree/16S/epa_placements_v2/epa_result.jplace --out-dir $WDIR/tree/16S/epa_placements_v2/tree_epa_concat --file-prefix 16S_tree_plusASVs_ETH_ --threads 60 --log-file $WDIR/tree/16S/epa_placements_v2/tree_epa_concat/logfile
