# Annotation of all EX MAGs

### Run diamond scan with KEGG against MAGs of interest

# In /storage/hdd6/DB/KEGG/release_20221204 already present are scan db and pepunit db (both of which are required)

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"
DB="/storage/hdd6/DB/KEGG/release_20221204"
conda activate checkm2

mkdir -p $WDIR/annotation/kegg/tmp
mkdir -p $WDIR/annotation/kegg/kegg_annotation_dmnd/
mkdir -p $WDIR/annotation/kegg/kegg_annotation_selfblast/
mkdir -p $WDIR/annotation/kegg/kegg_annotation_parsed/
THREADS=100

# scan DB
# use input from prodigal

ls -d $WDIR/mmseq2/prodigal_out/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  diamond blastp --db $DB/diamond/kegg_genes.dmnd --query $WDIR/mmseq2/prodigal_out/${line}.faa --out $WDIR/annotation/kegg/${line}".kegg.annotation" --evalue 1e-10 --very-sensitive --top 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads $THREADS --tmpdir $WDIR/annotation/kegg/tmp
  diamond makedb --in $WDIR/mmseq2/prodigal_out/${line}.faa -d $WDIR/annotation/kegg/kegg_annotation_dmnd/cds_db_${line}
  diamond blastp -d $WDIR/annotation/kegg/kegg_annotation_dmnd/cds_db_${line}.dmnd -q $WDIR/mmseq2/prodigal_out/${line}.faa -k 1 -o $WDIR/annotation/kegg/kegg_annotation_selfblast/${line}".selfblast" -f 6
done
 
 
 
module load R/4.2.2
ls -d $WDIR/annotation/kegg/*.kegg.annotation | sed "s/\.kegg.annotation//" | xargs -n 1 basename |  while read line
do
  $WDIR/parse_kegg.R -b $WDIR/annotation/kegg/${line}".kegg.annotation" -s $WDIR/annotation/kegg/kegg_annotation_selfblast/${line}".selfblast" -m $DB/mapping_info -c 0.4 -t 120 -o $WDIR/annotation/kegg/kegg_annotation_parsed_new/${line}".kegg.parsed"
done



# Prediction of peptidases 
# new MEROPS was released, therefore run new MEROPS and SignalP test, diamond 2.0.15
conda activate checkm2

mkdir -p $WDIR/annotation/MEROPS/scan $WDIR/annotation/MEROPS/pepunit $WDIR/annotation/MEROPS/tmp $WDIR/annotation/MEROPS/peptide_blast $WDIR/annotation/MEROPS/out_merops_parsed

DB="/storage/hdd6/DB/MEROPS/v12.4/diamond"
THREADS=100

# scan DB
# use input from prodigal

ls -d $WDIR/mmseq2/prodigal_out/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
 do
  diamond blastp --db $DB/merops_scan.dmnd --query $WDIR/mmseq2/prodigal_out/${line}.faa --out $WDIR/annotation/MEROPS/scan/${line}".scan.blastout" --very-sensitive --evalue 1e-10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads $THREADS
 done


ls -d $WDIR/mmseq2/prodigal_out/*.faa| sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  diamond blastp --db $DB/merops_pepunit.dmnd --query $WDIR/mmseq2/prodigal_out/${line}.faa --out  $WDIR/annotation/MEROPS/pepunit/${line}".pepunit.blastout" --very-sensitive --evalue 1e-10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads $THREADS
done


# check for false positives by only keeping those hits from pepunit that were also present in scan.db
ls -1 $WDIR/annotation/MEROPS/pepunit | sed 's/\.pepunit.blastout//' | while read line
do
 cat $WDIR/annotation/MEROPS/pepunit/${line}".pepunit.blastout" | cut -f1 | grep -F -f - <(cut -f1 $WDIR/annotation/MEROPS/scan/${line}".scan.blastout") | uniq  > $WDIR/annotation/MEROPS/tmp/${line}"_hits"
done


ls -1 $WDIR/annotation/MEROPS/tmp | sed 's/\_hits//' | while read line
do
 cat $WDIR/annotation/MEROPS/tmp/${line}"_hits" | while read IID
 do
  grep ${IID} $WDIR/annotation/MEROPS/pepunit/${line}.pepunit.blastout >> $WDIR/annotation/MEROPS/peptide_blast/${line}".peptide.blastout"
 done
done


# parse files
module load R/4.2.2
ls -d $WDIR/mmseq2/prodigal_out/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  ./parse_merops.R -f $WDIR/annotation/MEROPS/scan/${line}".scan.blastout" -p $WDIR/annotation/MEROPS/pepunit/${line}".pepunit.blastout" -m /storage/hdd6/DB/MEROPS/v12.4/metadata_pepunit.txt -s $WDIR/annotation/kegg/kegg_annotation_selfblast/${line}".selfblast" -c 0.4 -o $WDIR/annotation/MEROPS/out_merops_parsed/${line}".merops.parsed"
done


# Prediction of extracellular peptidases
# run signalp on samples
conda activate signalp-6.0g_fast
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"

mkdir -p $WDIR/annotation/signalp
cd $WDIR/annotation/signalp
ls -d $WDIR/mmseq2/prodigal_out/*.faa | sed "s/\.faa//"  | xargs -n 1 basename  | while read line
do
  signalp6 -fasta $WDIR/mmseq2/prodigal_out/${line}.faa -format none -org other -od $WDIR/annotation/signalp/${line}_signalp
done



# dbCAN annotation
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"
DB="/storage/hdd6/DB/dbCAN/v3.0.7"

conda activate dbcan-3.0.7

mkdir -p $WDIR/annotation/dbCAN
cd $WDIR/annotation/dbCAN
ls -d $WDIR/mmseq2/prodigal_out/*.faa | sed "s/\.faa//"  | xargs -n 1 basename |  while read line
do
   #sed '/^\\#\\#FASTA$/,$d'  $WDIR/Annotation/MetaErg/output/${line}/data/all.gff > $WDIR/dbCAN/tmp_${line}.gff
   run_dbcan $WDIR/mmseq2/prodigal_out/${line}.faa protein --dbCANFile $DB/dbCAN.txt --out_dir $WDIR/annotation/dbCAN/${line} --db_dir $DB >> $WDIR/logfiles/dbCAN_${line}.log 2>&1
done
rm $WDIR/annotation/dbCAN/tmp*



# NR annotation
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"
DB="/storage/hdd6/DB/Diamond/NR_130523"
conda activate checkm2 # diamond 2.0.15

mkdir -p $WDIR/annotation/NR/tmp
mkdir -p $WDIR/annotation/NR/NR_annotation_dmnd/
mkdir -p $WDIR/annotation/NR/NR_annotation_selfblast/
mkdir -p $WDIR/annotation/NR/NR_annotation_parsed/
THREADS=100

ls -d $WDIR/mmseq2/prodigal_out/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  diamond blastp --db $DB/diamond_nr.dmnd --query $WDIR/mmseq2/prodigal_out/${line}.faa --out $WDIR/annotation/NR/${line}".NR.annotation" --evalue 1e-10 --sensitive --top 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxids sscinames salltitles -b 10 --threads $THREADS --tmpdir $WDIR/annotation/NR/tmp
  diamond makedb --in $WDIR/mmseq2/prodigal_out/${line}.faa -d $WDIR/annotation/NR/NR_annotation_dmnd/cds_db_${line}
  diamond blastp -d $WDIR/annotation/NR/NR_annotation_dmnd/cds_db_${line}.dmnd -q $WDIR/mmseq2/prodigal_out/${line}.faa -k 1 -o $WDIR/annotation/NR/NR_annotation_selfblast/${line}".selfblast" -f 6
done

module load R/4.2.2
ls -d $WDIR/annotation/NR/*.NR.annotation | sed "s/\.NR.annotation//" | xargs -n 1 basename |  while read line
do
  $WDIR/parse_NR.R -b $WDIR/annotation/NR/${line}".NR.annotation" -s $WDIR/annotation/NR/NR_annotation_selfblast/${line}".selfblast"  -c 0.4 -o $WDIR/annotation/NR/NR_annotation_parsed_new/${line}".NR.parsed"
done


# Annotation of EX MAGs using the Ocean Genome Atlas

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"
DB="/storage/hdd6/DB/Diamond/OM-RGC_v2"

conda activate checkm2

mkdir -p $WDIR/annotation/OM-RGCv2/tmp
mkdir -p $WDIR/annotation/OM-RGCv2/OM-RGCv2_annotation_dmnd/
mkdir -p $WDIR/annotation/OM-RGCv2/OM-RGCv2_annotation_selfblast/
mkdir -p $WDIR/annotation/OM-RGCv2/OM-RGCv2_annotation_parsed/
THREADS=100

# annotation for EX MAGs
# use input from prodigal

ls -d $WDIR/mmseq2/prodigal_out/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  diamond blastp --db $DB/OM-RGC_v2_ref.dmnd --query $WDIR/mmseq2/prodigal_out/${line}.faa --out $WDIR/annotation/OM-RGCv2/${line}".OM-RGCv2.annotation" --evalue 1e-10 --very-sensitive --top 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads $THREADS --tmpdir $WDIR/annotation/OM-RGCv2/tmp
done


# parse annotation
ls -d $WDIR/annotation/OM-RGCv2/*.OM-RGCv2.annotation | sed "s/\.OM-RGCv2.annotation//" | xargs -n 1 basename |  while read line
do
  $WDIR/parse_OM-RGCv2.R -b $WDIR/annotation/OM-RGCv2/OM-RGCv2/${line}".OM-RGCv2.annotation" -s $WDIR/annotation/NR/NR_annotation_selfblast/${line}".selfblast" -m $DB/OM-RGC_v2.tsv -c 0.4 -t 60 -o $WDIR/annotation/OM-RGCv2/OM-RGCv2_annotation_parsed/${line}".OM-RGCv2.parsed"
done


# Transporter using TCDB database
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"
DB="/storage/hdd6/DB/TCDB/v122023/diamond"
conda activate checkm2

mkdir -p $WDIR/annotation/TCDB/tmp
mkdir -p $WDIR/annotation/TCDB/TCDB_annotation_parsed/

THREADS=100

# EX MAGs
ls -d $WDIR/mmseq2/prodigal_out/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  diamond blastp --db $DB/tcdb.dmnd --query $WDIR/mmseq2/prodigal_out/${line}.faa --out $WDIR/annotation/TCDB/${line}".TCDB.annotation" --evalue 1e-10 --very-sensitive --top 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads $THREADS --tmpdir $WDIR/annotation/TCDB/tmp
done

# parse
module load R/4.2.2
ls -d $WDIR/annotation/TCDB/*.TCDB.annotation | sed "s/\.TCDB.annotation//" | xargs -n 1 basename |  while read line
do
  $WDIR/parse_TCDB.R -b $WDIR/annotation/TCDB/${line}".TCDB.annotation" -s $WDIR/annotation/kegg/kegg_annotation_selfblast/${line}".selfblast" -m /storage/hdd6/DB/TCDB/v122023 -c 0.4 -t 120 -o $WDIR/annotation/TCDB/TCDB_annotation_parsed/${line}".TCDB.parsed"
done


##### Repeat annotation also for all ETH_MAGs

### Run diamond scan with KEGG against ETH MAGs of interest

# In /storage/hdd6/DB/KEGG/release_20221204 already present are scan db and pepunit db (both of which are required)

WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"
DB="/storage/hdd6/DB/KEGG/release_20221204"
conda activate checkm2

mkdir -p $WDIR/annotation_ETH/kegg/tmp
mkdir -p $WDIR/annotation_ETH/kegg/kegg_annotation_dmnd/
mkdir -p $WDIR/annotation_ETH/kegg/kegg_annotation_selfblast/
mkdir -p $WDIR/annotation_ETH/kegg/kegg_annotation_parsed/
THREADS=100

# scan DB
# use input from prodigal

ls -d $WDIR/ETH_MAGs/ETH_selection_final/prodigal/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  diamond blastp --db $DB/diamond/kegg_genes.dmnd --query $WDIR/ETH_MAGs/ETH_selection_final/prodigal/${line}.faa --out $WDIR/annotation_ETH/kegg/${line}".kegg.annotation" --evalue 1e-10 --very-sensitive --top 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads $THREADS --tmpdir $WDIR/annotation_ETH/kegg/tmp
  diamond makedb --in $WDIR/ETH_MAGs/ETH_selection_final/prodigal/${line}.faa -d $WDIR/annotation_ETH/kegg/kegg_annotation_dmnd/cds_db_${line}
  diamond blastp -d $WDIR/annotation_ETH/kegg/kegg_annotation_dmnd/cds_db_${line}.dmnd -q $WDIR/ETH_MAGs/ETH_selection_final/prodigal/${line}.faa -k 1 -o $WDIR/annotation_ETH/kegg/kegg_annotation_selfblast/${line}".selfblast" -f 6
done
 
module load R/4.2.2
ls -d $WDIR/ETH_MAGs/ETH_selection_final/prodigal/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  $WDIR/parse_kegg.R -b $WDIR/annotation_ETH/kegg/${line}".kegg.annotation" -s $WDIR/annotation_ETH/kegg/kegg_annotation_selfblast/${line}".selfblast" -m $DB/mapping_info -c 0.4 -t 120 -o $WDIR/annotation_ETH/kegg/kegg_annotation_parsed/${line}".kegg.parsed"
done


# Prediction of peptidases 
# new MEROPS was released, therefore run new MEROPS and SignalP test, diamond 2.0.15
conda activate checkm2

mkdir -p $WDIR/annotation_ETH/MEROPS/scan $WDIR/annotation_ETH/MEROPS/pepunit $WDIR/annotation_ETH/MEROPS/tmp $WDIR/annotation_ETH/MEROPS/peptide_blast $WDIR/annotation_ETH/MEROPS/out_merops_parsed

DB="/storage/hdd6/DB/MEROPS/v12.4/diamond"
THREADS=100

# scan DB
# use input from prodigal

ls -1 $WDIR/ETH_MAGs/ETH_selection_final/prodigal/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
 do
  diamond blastp --db $DB/merops_scan.dmnd --query $WDIR/ETH_MAGs/ETH_selection_final/prodigal/${line}.faa --out $WDIR/annotation_ETH/MEROPS/scan/${line}".scan.blastout" --very-sensitive --evalue 1e-10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads $THREADS
 done

ls -1 $WDIR/ETH_MAGs/ETH_selection_final/prodigal/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  diamond blastp --db $DB/merops_pepunit.dmnd --query $WDIR/ETH_MAGs/ETH_selection_final/prodigal/${line}.faa --out  $WDIR/annotation_ETH/MEROPS/pepunit/${line}".pepunit.blastout" --very-sensitive --evalue 1e-10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads $THREADS
done


# check for false positives by only keeping those hits from pepunit that were also present in scan.db
ls -1 $WDIR/annotation_ETH/MEROPS/pepunit | sed 's/\.pepunit.blastout//' | while read line
do
 cat $WDIR/annotation_ETH/MEROPS/pepunit/${line}".pepunit.blastout" | cut -f1 | grep -F -f - <(cut -f1 $WDIR/annotation_ETH/MEROPS/scan/${line}".scan.blastout") | uniq  > $WDIR/annotation_ETH/MEROPS/tmp/${line}"_hits"
done


ls -1 $WDIR/annotation_ETH/MEROPS/tmp | sed 's/\_hits//' | while read line
do
 cat $WDIR/annotation_ETH/MEROPS/tmp/${line}"_hits" | while read IID
 do
  grep ${IID} $WDIR/annotation_ETH/MEROPS/pepunit/${line}.pepunit.blastout >> $WDIR/annotation_ETH/MEROPS/peptide_blast/${line}".peptide.blastout"
 done
done


# parse files
module load R/4.2.2
ls -1 $WDIR/ETH_MAGs/ETH_selection_final/prodigal/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  $WDIR/parse_merops.R -f $WDIR/annotation_ETH/MEROPS/scan/${line}".scan.blastout" -p $WDIR/annotation_ETH/MEROPS/pepunit/${line}".pepunit.blastout" -m /storage/hdd6/DB/MEROPS/v12.4/metadata_pepunit.txt -s $WDIR/annotation_ETH/kegg/kegg_annotation_selfblast/${line}".selfblast" -c 0.4 -o $WDIR/annotation_ETH/MEROPS/out_merops_parsed/${line}".merops.parsed"
done


# Prediction of extracellular peptidases
# run signalp on samples
conda activate signalp-6.0g_fast
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"

mkdir -p $WDIR/annotation_ETH/signalp
cd $WDIR/annotation_ETH/signalp
ls -d $WDIR/ETH_MAGs/ETH_selection_final/prodigal/*.faa | sed "s/\.faa//"  | xargs -n 1 basename  | while read line
do
  signalp6 -fasta $WDIR/ETH_MAGs/ETH_selection_final/prodigal/${line}.faa -format none -org other -od $WDIR/annotation_ETH/signalp/${line}_signalp
done


# dbCAN annotation
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"
DB="/storage/hdd6/DB/dbCAN/v3.0.7"

conda activate dbcan-3.0.7

mkdir -p $WDIR/annotation_ETH/dbCAN
cd $WDIR/annotation_ETH/dbCAN
ls -d $WDIR/ETH_MAGs/ETH_selection_final/prodigal/*.faa | sed "s/\.faa//"  | xargs -n 1 basename |  while read line
do
   #sed '/^\\#\\#FASTA$/,$d'  $WDIR/Annotation/MetaErg/output/${line}/data/all.gff > $WDIR/dbCAN/tmp_${line}.gff
   run_dbcan $WDIR/ETH_MAGs/ETH_selection_final/prodigal/${line}.faa protein --dbCANFile $DB/dbCAN.txt --out_dir $WDIR/annotation_ETH/dbCAN/${line} --db_dir $DB >> $WDIR/logfiles/dbCAN_${line}.log 2>&1
done
rm $WDIR/annotation_ETH/dbCAN/tmp*



# NR annotation
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"
DB="/storage/hdd6/DB/Diamond/NR_130523"
conda activate checkm2 # diamond 2.0.15

mkdir -p $WDIR/annotation_ETH/NR/tmp
mkdir -p $WDIR/annotation_ETH/NR/NR_annotation_dmnd/
mkdir -p $WDIR/annotation_ETH/NR/NR_annotation_selfblast/
mkdir -p $WDIR/annotation_ETH/NR/NR_annotation_parsed/
THREADS=60

ls -d $WDIR/ETH_MAGs/ETH_selection_final/prodigal/*.faa| sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  diamond blastp --db $DB/diamond_nr.dmnd --query $WDIR/ETH_MAGs/ETH_selection_final/prodigal/${line}.faa --out $WDIR/annotation_ETH/NR/${line}".NR.annotation" --evalue 1e-10 --sensitive --top 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle staxids sscinames salltitles -b 10 --threads $THREADS --tmpdir $WDIR/annotation_ETH/NR/tmp
  diamond makedb --in $WDIR/ETH_MAGs/ETH_selection_final/prodigal/${line}.faa  -d $WDIR/annotation_ETH/NR/NR_annotation_dmnd/cds_db_${line}
  diamond blastp -d $WDIR/annotation_ETH/NR/NR_annotation_dmnd/cds_db_${line}.dmnd -q $WDIR/ETH_MAGs/ETH_selection_final/prodigal/${line}.faa -k 1 -o $WDIR/annotation_ETH/NR/NR_annotation_selfblast/${line}".selfblast" -f 6
done

module load R/4.2.2
ls -d $WDIR/ETH_MAGs/ETH_selection_final/prodigal/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  $WDIR/parse_NR.R -b $WDIR/annotation_ETH/NR/${line}".NR.annotation" -s $WDIR/annotation_ETH/NR/NR_annotation_selfblast/${line}".selfblast"  -c 0.4 -o $WDIR/annotation_ETH/NR/NR_annotation_parsed/${line}".NR.parsed"
done

ls -1 $WDIR/annotation/NR/NR_annotation_parsed_new/*best.txt | sed 's/\.NR.parsed_best.txt//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/annotation/NR/NR_annotation_parsed_new/${line}".NR.parsed_best.txt" >> $WDIR/annotation/NR/NR_annotation_parsed_new/concat_test.txt
done


# Annotation of ETH MAGs using the Ocean Genome Atlas

mkdir -p $WDIR/annotation_ETH/OM-RGCv2/tmp
mkdir -p $WDIR/annotation_ETH/OM-RGCv2/OM-RGCv2_annotation_dmnd/
mkdir -p $WDIR/annotation_ETH/OM-RGCv2/OM-RGCv2_annotation_selfblast/
mkdir -p $WDIR/annotation_ETH/OM-RGCv2/OM-RGCv2_annotation_parsed/

ls -d $WDIR/ETH_MAGs/ETH_selection_final/prodigal/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  diamond blastp --db $DB/OM-RGC_v2_ref.dmnd --query $WDIR/ETH_MAGs/ETH_selection_final/prodigal/${line}.faa --out $WDIR/annotation_ETH/OM-RGCv2/${line}".OM-RGCv2.annotation" --evalue 1e-10 --very-sensitive --top 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads $THREADS --tmpdir $WDIR/annotation_ETH/OM-RGCv2/tmp
done

# parse annotation
ls -d $WDIR/annotation_ETH/OM-RGCv2/*.OM-RGCv2.annotation | sed "s/\.OM-RGCv2.annotation//" | xargs -n 1 basename |  while read line
do
  $WDIR/parse_OM-RGCv2.R -b $WDIR/annotation_ETH/OM-RGCv2/OM-RGCv2/${line}".OM-RGCv2.annotation" -s $WDIR/annotation_ETH/NR/NR_annotation_selfblast/${line}".selfblast" -m $DB/OM-RGC_v2.tsv -c 0.4 -t 60 -o $WDIR/annotation_ETH/OM-RGCv2/OM-RGCv2_annotation_parsed/${line}".OM-RGCv2.parsed"
done


# Transporter using TCDB database
WDIR="/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs"
DB="/storage/hdd6/DB/TCDB/v122023/diamond"
conda activate checkm2

mkdir -p $WDIR/annotation_ETH/TCDB/tmp
mkdir -p $WDIR/annotation_ETH/TCDB/TCDB_annotation_parsed/

THREADS=100

# ETH MAGs
ls -d $WDIR/ETH_MAGs/ETH_selection_final/prodigal/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  diamond blastp --db $DB/tcdb.dmnd --query $WDIR/ETH_MAGs/ETH_selection_final/prodigal/${line}.faa --out $WDIR/annotation_ETH/TCDB/${line}".TCDB.annotation" --evalue 1e-10 --very-sensitive --top 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads $THREADS --tmpdir $WDIR/annotation_ETH/TCDB/tmp
done

# parse
module load R/4.2.2
ls -d $WDIR/ETH_MAGs/ETH_selection_final/prodigal/*.faa | sed "s/\.faa//" | xargs -n 1 basename |  while read line
do
  $WDIR/parse_TCDB.R -b $WDIR/annotation_ETH/TCDB/${line}".TCDB.annotation" -s $WDIR/annotation_ETH/kegg/kegg_annotation_selfblast/${line}".selfblast" -m /storage/hdd6/DB/TCDB/v122023 -c 0.4 -t 120 -o $WDIR/annotation_ETH/TCDB/TCDB_annotation_parsed/${line}".TCDB.parsed"
done



# parse results from annotations

# kegg annotation
ls -1 $WDIR/annotation/kegg/kegg_annotation_parsed_new/*best.txt | sed 's/\.kegg.parsed_best.txt//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/annotation/kegg/kegg_annotation_parsed_new/${line}".kegg.parsed_best.txt" >> $WDIR/annotation_ETH/concat_kegg_annotation_best.txt
done

ls -1 $WDIR/annotation_ETH/kegg/kegg_annotation_parsed2/*best.txt | sed 's/\.kegg.parsed_best.txt//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/annotation_ETH/kegg/kegg_annotation_parsed2/${line}".kegg.parsed_best.txt" >> $WDIR/annotation_ETH/concat_kegg_annotation_best.txt
done

# merops
ls -1 $WDIR/annotation/MEROPS/out_merops_parsed/*best.txt | sed 's/\.merops.parsed_best.txt//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/annotation/MEROPS/out_merops_parsed/${line}".merops.parsed_best.txt" >> $WDIR/annotation_ETH/concat_merops_annotation_best.txt
done

ls -1 $WDIR/annotation_ETH/MEROPS/out_merops_parsed/*best.txt | sed 's/\.merops.parsed_best.txt//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/annotation_ETH/MEROPS/out_merops_parsed/${line}".merops.parsed_best.txt" >> $WDIR/annotation_ETH/concat_merops_annotation_best.txt
done

# signalp
ls -d $WDIR/mmseq2/prodigal_out/*.faa | sed "s/\.faa//" | xargs -n 1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/annotation/signalp/${line}"_signalp"/prediction_results.txt >> $WDIR/annotation_ETH/concat_signalp_annotation.txt
done

ls -1 $WDIR/ETH_MAGs/ETH_selection_final/prodigal/*.faa | sed "s/\.faa//" | xargs -n 1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/annotation_ETH/signalp/${line}"_signalp"/prediction_results.txt >> $WDIR/annotation_ETH/concat_signalp_annotation.txt
done

# dbCAN
ls -d $WDIR/mmseq2/prodigal_out/*.faa | sed "s/\.faa//" | xargs -n 1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/annotation/dbCAN/${line}/overview.txt >> $WDIR/annotation_ETH/concat_dbCAN_annotation.txt
done

ls -1 $WDIR/ETH_MAGs/ETH_selection_final/prodigal/*.faa | sed "s/\.faa//" | xargs -n 1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/annotation_ETH/dbCAN/${line}/overview.txt >> $WDIR/annotation_ETH/concat_dbCAN_annotation.txt
done

# NR
ls -1 $WDIR/annotation/NR/NR_annotation_parsed_new/*best.txt | sed 's/\.NR.parsed_best.txt//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/annotation/NR/NR_annotation_parsed_new/${line}".NR.parsed_best.txt" >> $WDIR/annotation_ETH/concat_NR_annotation_best.txt
done

ls -1 $WDIR/annotation_ETH/NR/NR_annotation_parsed2/*best.txt | sed 's/\.NR.parsed_best.txt//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/annotation_ETH/NR/NR_annotation_parsed2/${line}".NR.parsed_best.txt" >> $WDIR/annotation_ETH/concat_NR_annotation_best.txt
done

# OM-RGC-v2
ls -1 $WDIR/annotation/OM-RGCv2/OM-RGCv2_annotation_parsed/*best.txt | sed 's/\.OM-RGCv2.parsed_best.txt//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/annotation/OM-RGCv2/OM-RGCv2_annotation_parsed/${line}".OM-RGCv2.parsed_best.txt" >> $WDIR/annotation_ETH/concat_OM-RGCv2_annotation_best.txt
done

ls -1 $WDIR/annotation_ETH/OM-RGCv2/OM-RGCv2_annotation_parsed/*best.txt | sed 's/\.OM-RGCv2.parsed_best.txt//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/annotation_ETH/OM-RGCv2/OM-RGCv2_annotation_parsed/${line}".OM-RGCv2.parsed_best.txt" >> $WDIR/annotation_ETH/concat_OM-RGCv2_annotation_best.txt
done


# TCDB
ls -1 $WDIR/annotation/TCDB/TCDB_annotation_parsed/*best.txt | sed 's/\.TCDB.parsed_best.txt//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/annotation/TCDB/TCDB_annotation_parsed/${line}".TCDB.parsed_best.txt" >> $WDIR/annotation_ETH/concat_TCDB_annotation_best.txt
done

ls -1 $WDIR/annotation_ETH/TCDB/TCDB_annotation_parsed/*best.txt | sed 's/\.TCDB.parsed_best.txt//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" $WDIR/annotation_ETH/TCDB/TCDB_annotation_parsed/${line}".TCDB.parsed_best.txt" >> $WDIR/annotation_ETH/concat_TCDB_annotation_best.txt
done


# concatenate InterProscan results
ls -1 /storage/hdd1/chh/Mara_metaG_data_mining/EX_annotation/interpro_out/*tsv | sed 's/\.iprout.tsv//' | xargs -n1 basename | while read line
do
	sed "s/^/${line}___/" /storage/hdd1/chh/Mara_metaG_data_mining/EX_annotation/interpro_out/${line}".iprout.tsv" >> /storage/hdd1/chh/Mara_metaG_data_mining/EX_annotation/interpro_out/concat_interproscan_annotation.txt
done
