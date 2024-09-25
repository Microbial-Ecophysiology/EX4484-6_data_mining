#!/bin/bash

SID=$1
CPU=$2
PPLACER_CPU=$3
complete=$4
contam=$5

# running last steps in bin refinement manually
mkdir $SID/tmp
checkm lineage_wf -x fa $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_bins $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_bins.checkm -t $CPU --tmpdir $SID/tmp --pplacer_threads $PPLACER_CPU
# this duplicates the above error
# shortening path to tmpdir --> this seems to be very annoying behavior that we just have to live with
/opt/moep/metaWRAP/bin/metawrap-scripts/summarize_checkm.py $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_bins.checkm/storage/bin_stats_ext.tsv | (read -r; printf "%s\n" "$REPLY"; sort) > $SID/metaG_M>rm -r $SID/tmp

# collect best bins
mkdir -p $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_best_bins
for i in $(/opt/moep/metaWRAP/bin/metawrap-scripts/choose_best_bin.py $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_bins.stats $complete $contam)
do
  echo "Copying best bin: $i"
  cp $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_bins/${i}.fa $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_best_bins
done
ls -l $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_best_bins | grep orig | wc -l
ls -l $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_best_bins | grep strict | wc -l
ls -l $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_best_bins | grep permissive | wc -l

mkdir $SID/metaG_MAGs/reassembled_bins/metawrap/work_files
mv $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_bins $SID/metaG_MAGs/reassembled_bins/metawrap/work_files/
mv $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_bins.checkm $SID/metaG_MAGs/reassembled_bins/metawrap/work_files/
mv $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_bins.stats $SID/metaG_MAGs/reassembled_bins/metawrap/work_files/
mv $SID/metaG_MAGs/reassembled_bins/metawrap/reads_for_reassembly $SID/metaG_MAGs/reassembled_bins/metawrap/work_files/
mv $SID/metaG_MAGs/reassembled_bins/metawrap/binned_assembly $SID/metaG_MAGs/reassembled_bins/metawrap/work_files/
mv $SID/metaG_MAGs/reassembled_bins/metawrap/reassemblies $SID/metaG_MAGs/reassembled_bins/metawrap/work_files/
mv $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_best_bins $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_bins

# re-running checkm on final bin set
mkdir $SID/tmp
checkm lineage_wf -x fa $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_bins $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_bins.checkm -t $CPU --tmpdir $SID/tmp --pplacer_threads $PPLACER_CPU
/opt/moep/metaWRAP/bin/metawrap-scripts/summarize_checkm.py $SID/metaG_MAGs/reassembled_bins/metawrap/reassembled_bins.checkm/storage/bin_stats_ext.tsv | (read -r; printf "%s\n" "$REPLY"; sort) >  $SID/metaG_>rm -r $SID/tmp
