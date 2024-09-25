# Script to combine annotation results with ETH MAGs

library(tidyverse)
library(reshape2)

setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Data_mining/Annotation_v3")

# load annotation ####
load("annotation_parsed.RData")

## kegg ----
kegg <- read.csv("full/concat_kegg_annotation_best.txt", 
                 sep="\t", 
                 header=TRUE, 
                 fill=TRUE, 
                 quote="",
                 col.names=c("qseqid", "kegg_sseqid","kegg_pident","kegg_length","kegg_mismatch","kegg_gapopen","kegg_qstart",
                             "kegg_qend","kegg_sstart","kegg_send","kegg_evalue","kegg_bitscore","kegg_bsr","kegg_KO","kegg_KO_name","kegg_taxon","kegg_gene"))
kegg <- subset(kegg, kegg_sseqid != "sseqid")

## NR ----
NR <- read.csv("full/concat_NR_annotation_best.txt", 
               sep="\t", 
               header=TRUE, 
               fill=TRUE, 
               quote="",
               col.names=c("qseqid", "NR_sseqid","NR_pident","NR_length","NR_mismatch","NR_gapopen","NR_qstart",
                           "NR_qend","NR_sstart","NR_send","NR_evalue","NR_bitscore","NR_qlen","NR_slen","NR_stitle","NR_staxids","NR_sscinames","NR_salltitles","NR_bsr"))
NR <- subset(NR, NR_sseqid != "sseqid")

## OM-RGCv2 ----
OM_RGC <- read.csv("full/concat_OM-RGCv2_annotation_best.txt", 
               sep="\t", 
               header=TRUE, 
               fill=TRUE, 
               quote="",
               col.names=c("qseqid", "OM_RGC_sseqid","OM_RGC_pident","OM_RGC_length","OM_RGC_mismatch","OM_RGC_gapopen","OM_RGC_qstart",
                           "OM_RGC_qend","OM_RGC_sstart","OM_RGC_send","OM_RGC_evalue","OM_RGC_bitscore","OM_RGC_bsr","OM_RGC_KO","OM_RGC_COG","OM_RGC_taxpath"))
OM_RGC <- subset(OM_RGC, OM_RGC_sseqid != "sseqid")


## dbCAN ----
dbCAN <- read.csv("CAZy/concat_dbCAN_annotation.txt", 
                  sep="\t", 
                  header=TRUE, 
                  fill=TRUE, 
                  quote="",
                  col.names=c("Gene.ID","dbCAN_EC","dbCAN_HMMER","dbCAN_eCAMI","dbCAN_DIAMOND","dbCAN_no.Tools"))
dbCAN$Gene.ID <- gsub(".faa","",as.character(dbCAN$Gene.ID))

## merops ----
merops <- read.csv("Peptidases/concat_merops_annotation_best.txt", 
                   sep="\t", 
                   header=TRUE, 
                   fill=TRUE, 
                   quote="",
                   col.names=c("qseqid", "merops_sseqid","merops_pident","merops_length","merops_mismatch","merops_gapopen","merops_qstart",
                               "merops_qend","merops_sstart","merops_send","merops_evalue","merops_bitscore","merops_bsr","merops_name","merops_id","merops_subfamily","merops_unit","merops_source"))

merops <- subset(merops, merops_sseqid != "sseqid")
merops$qseqid <- gsub(".faa","",as.character(merops$qseqid))
merops$qseqid <- gsub(".genes","",as.character(merops$qseqid))

## signalp ----
#first line of output and # in second line need to be removed before loading into R
signalp <- read.csv("Peptidases/concat_signalp_annotation.txt", 
                    sep="\t", 
                    header=TRUE, 
                    fill=TRUE, 
                    quote="")
signalp$ID <- gsub(".faa","",as.character(signalp$ID))
signalp$ID <- gsub(".genes","",as.character(signalp$ID))


## TCDB
TCDB <- read.csv("TCDB/concat_TCDB_annotation_best.txt", 
                   sep="\t", 
                   header=TRUE, 
                   fill=TRUE, 
                   quote="",
                   col.names=c("qseqid", "TCDB_sseqid","TCDB_pident","TCDB_length","TCDB_mismatch","TCDB_gapopen","TCDB_qstart",
                               "TCDB_qend","TCDB_sstart","TCDB_send","TCDB_evalue","TCDB_bitscore","TCDB_bsr","TCDB_tc_id","TCDB_tc_fam","TCDB_family","TCDB_go","TCDB_pfam","TCDB_chebi"))
TCDB <- subset(TCDB, TCDB_sseqid != "sseqid")


# Parse annotation ----
# Combine annotation data of kegg, NR, dbCAN and merops
annotation <- full_join(kegg, NR, by = c("qseqid")) %>%
  full_join(., OM_RGC, by = c("qseqid")) %>%
  full_join(., dbCAN, by = c("qseqid" = "Gene.ID")) %>%
  full_join(., merops, by = "qseqid", suffix = c("_kegg","_merops")) %>%
  full_join(.,TCDB, by ="qseqid") %>%
  separate(col = qseqid, into = c("genome","cluster"), sep = "___")

annotation$genome <- gsub(".genes","",as.character(annotation$genome))


write.table(annotation, "annotation_parsed.txt", 
            quote = FALSE,
            sep="\t",
            col.names = TRUE,
            row.names = FALSE)
