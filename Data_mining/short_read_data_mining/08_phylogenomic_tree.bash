# Run marker gene tree, also including all MAGs from co-authors
# cp own EX MAG, EX MAGs from MAG mining, IOW MAGs, ETH MAGs plus Outgroup MAGs to /storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs/tree/marker_gene_v3/input
# in total 432 MAGs incl. outgroup, novel EX MAGs with checkm2 completeness <80 were removed from MAG list, also two MAGs classified as bacteria
# ones we already have from ETH: SPET22-1_SAMN18353074_MAG_00000107.fa, SPET22-1_SAMN18353075_MAG_00000011.fa, SPET22-1_SAMN18353078_MAG_00000111.fa, SPET22-1_SAMN18353079_MAG_00000036.fa, DOMB MAGs all Guaymas basin; ZHAN22-2_SAMN18200488_MAG_00000104.fa (cold seep), were removed from selection, 13 MAGs remain

## run gtdbtk de novo workflow for a multiple sequence alignment
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"
conda activate gtdbtk-2.1.0
THREADS=30
gtdbtk de_novo_wf --genome_dir $WDIR/tree/marker_gene_v3/input --archaea --outgroup_taxon p__Halobacteriota --out_dir $WDIR/tree/marker_gene_v3/gtdbtk_de_novo_out --cpus $THREADS --prot_model LG --extension fa 

## run modeltest
conda activate modeltest-ng-0.1.7
mkdir -p $WDIR/tree/marker_gene_v3/modeltest
THREADS=100
modeltest-ng -i $WDIR/tree/marker_gene_v3/gtdbtk_de_novo_out/align/gtdbtk.ar53.user_msa.fasta -o $WDIR/tree/marker_gene_v3/modeltest/model_test -d aa -p $THREADS -t ml
# always results in error
# LH mismatch: 0.000000000000  != -2997523.715124555863
#modeltest-ng: /opt/conda/conda-bld/modeltest-ng_1646252428554/work/libs/pll-modules/src/algorithm/algo_search.c:1282: #pllmod_algo_spr_round: Assertion `fabs(loglh - best_lh) < 1e-6' failed.
#Aborted
# try running iqtree modeltest with program ModelFinder

conda activate iqtree2-2.2.2.7
iqtree2 -s $WDIR/tree/marker_gene_v3/gtdbtk_de_novo_out/align/gtdbtk.ar53.user_msa.fasta -m MF -T 60 
#415 sequences failed composition chi2 test (p-value<5%; df=19)

# best model LG+F+I+R10

## run raxml-ng
#check if MSA can be read, prefix for avoiding oerwriting of other output files
#MSA will be stored in binary format
#memory requirements will be estimated

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"
conda activate raxml-ng-1.1.0
mkdir -p $WDIR/tree/marker_gene_v3/raxml
cd $WDIR/tree/marker_gene_v3/raxml
raxml-ng --check --msa $WDIR/tree/marker_gene_v3/gtdbtk_de_novo_out/align/gtdbtk.ar53.user_msa.fasta --msa-format fasta --model LG+F+I+R10 --prefix T1

raxml-ng --parse --msa $WDIR/tree/marker_gene_v3/gtdbtk_de_novo_out/align/gtdbtk.ar53.user_msa.fasta --msa-format fasta --model LG+F+I+R10 --prefix T2
#25 threads recommended

# infer tree with model from modeltest
raxml-ng --msa $WDIR/tree/marker_gene_v3/gtdbtk_de_novo_out/align/gtdbtk.ar53.user_msa.fasta --msa-format fasta --model LG+F+I+R10 --prefix T3 --threads 80 --seed 2 
# infer more starting trees to check if this gives better results
raxml-ng --msa $WDIR/tree/marker_gene_v3/gtdbtk_de_novo_out/align/gtdbtk.ar53.user_msa.fasta --msa-format fasta --model LG+F+I+R10 --prefix T4 --tree pars{15},rand{15} --threads 80 --seed 2 
# run one quick and dirty tree
raxml-ng --search1 --msa $WDIR/tree/marker_gene_v3/gtdbtk_de_novo_out/align/gtdbtk.ar53.user_msa.fasta --msa-format fasta --model LG+F+I+R10 --prefix T5 --threads 80 --seed 2 

grep "Final LogLikelihood:" T{3,4,5,6}.raxml.log
# T3.raxml.log:Final LogLikelihood: -2432462.353340
# T4.raxml.log:Final LogLikelihood: -2432462.364990
# T5.raxml.log:Final LogLikelihood: -2432462.388813
# T6.raxml.log:Final LogLikelihood: -2432462.364990

# looks like T4 and T6 have same likelihood, but using only 20 starting trees gave best results. Checking for distinct topologies

# check for distinct topologies of tree T3
raxml-ng --rfdist --tree T3.raxml.mlTrees --prefix RF3
# Average absolute RF distance in this tree set: 7.568421
# Average relative RF distance in this tree set: 0.008821
# Number of unique topologies in this tree set: 11

# check for distinct topologies of tree T4
raxml-ng --rfdist --tree T4.raxml.mlTrees --prefix RF4
# Average absolute RF distance in this tree set: 7.668966
# Average relative RF distance in this tree set: 0.008938
# Number of unique topologies in this tree set: 15

# check for distinct topologies of tree T6
raxml-ng --rfdist --tree T6.raxml.mlTrees --prefix RF6
# Average absolute RF distance in this tree set: 7.620513
# Average relative RF distance in this tree set: 0.008882
# Number of unique topologies in this tree set: 14

# looks like T3 has least number of distinct topologies and best likelihood. Going with T3

# Run bootstraps
raxml-ng --bootstrap --msa $WDIR/tree/marker_gene_v3/gtdbtk_de_novo_out/align/gtdbtk.ar53.user_msa.fasta --msa-format fasta --model LG+F+I+R10 --prefix T7 --seed 2 --threads 100 --bs-trees 100

#check if bootstrapping reached convergence
raxml-ng --bsconverge --bs-trees T7.raxml.bootstraps --prefix T8 --seed 2 --threads 20 --bs-cutoff 0.03
#reached after 50
raxml-ng --bsconverge --bs-trees T7.raxml.bootstraps --prefix T9 --seed 2 --threads 14 --bs-cutoff 0.02
#reached after 100

#Calculating branch support using the best Raxml tree from T3 (Felsenstein bootstrap, FBP)
raxml-ng --support --tree T3.raxml.bestTree --bs-trees T7.raxml.bootstraps --prefix T10 --threads 14 
#also calculate transfer bootstrap expectation support metric, which is more appropirate for very large trees Transfer bootstrap TBE) # using this tree for phylogenomic tree
raxml-ng --support --tree T4.raxml.bestTree --bs-trees T7.raxml.bootstraps --prefix T11 --threads 14 --bs-metric tbe
