# Script to check Transporter in MAGs

library(tidyverse)
library(reshape2)

setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Data_mining/Annotation_v3")

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

## TCDB
TCDB <- read.csv("TCDB/concat_TCDB_annotation_best.txt", 
                   sep="\t", 
                   header=TRUE, 
                   fill=TRUE, 
                   quote="",
                   col.names=c("qseqid", "TCDB_sseqid","TCDB_pident","TCDB_length","TCDB_mismatch","TCDB_gapopen","TCDB_qstart",
                               "TCDB_qend","TCDB_sstart","TCDB_send","TCDB_evalue","TCDB_bitscore","TCDB_bsr","TCDB_tc_id","TCDB_tc_fam","TCDB_family","TCDB_go","TCDB_pfam","TCDB_chebi"))
TCDB <- subset(TCDB, TCDB_sseqid != "sseqid")


# heatmap of Transporters
# combine table with Transporters in first few columns
Transporters_table <- full_join(TCDB, kegg, by = c("qseqid")) %>%
  full_join(., NR, by = c("qseqid")) %>%
  separate(col = qseqid, into = c("genome","cluster"), sep = "___")

# remove all rows from dataframe which don't have any annotation through NR (or kegg)as NR has annotation for all
Transporters_short <- Transporters_table[complete.cases(Transporters_table[, 3:17]) & complete.cases(Transporters_table[, 37]), ]
# 1375 Transporters remain

write.table(Transporters_short, "TCDB/Transporters_short_parsed.txt", 
            quote = FALSE,
            sep="\t",
            col.names = TRUE,
            row.names = FALSE)

## Create a matrix of TCDB families against genomes ####
matrix_TCDB <- table(Transporters_short$genome, Transporters_short$TCDB_tc_fam, Transporters_short$TCDB_tc_id, Transporters_short$TCDB_family)
matrix_TCDB_short <- table(Transporters_short$genome, Transporters_short$TCDB_tc_fam, Transporters_short$TCDB_family)
# melt matrix into tabular view
matrix_TCDB_melt <- melt(matrix_TCDB)
matrix_TCDB_melt_short <- melt(matrix_TCDB_short)
colnames(matrix_TCDB_melt) <- c("genome","tc_fam","tc_id","TCDB_family","value")
colnames(matrix_TCDB_melt_short) <- c("genome","tc_fam","TCDB_family","value")
matrix_TCDB_melt <- matrix_TCDB_melt[- grep("PRJNA713414", matrix_TCDB_melt$genome),]
matrix_TCDB_melt_short <- matrix_TCDB_melt_short[- grep("PRJNA713414", matrix_TCDB_melt_short$genome),]
matrix_TCDB_melt$genome <- gsub(".genes","",as.character(matrix_TCDB_melt$genome))
matrix_TCDB_melt_short$genome <- gsub(".genes","",as.character(matrix_TCDB_melt_short$genome))

matrix_TCDB_melt$genome <- factor(matrix_TCDB_melt$genome, levels = c("SUTE22-1_SAMN10231914_MAG_00000186","SUTE22-1_SAMN10231894_MAG_00000173","SUTE22-1_SAMN10231913_MAG_00000219","SUTE22-1_SAMN10231893_MAG_00000067",
                                                                    "SUTE22-1_SAMN10231914_MAG_00000292","SUTE22-1_SAMN10231913_MAG_00000177",
                                                                    "SUTE22-1_SAMN10231904_MAG_00000012","SCRA20-1_SAMN11854494_MAG_00000079","PRJNA704804_bin_147_orig_refined-contigs","TARA_SAMEA2620113_MAG_00000097",
                                                                    "PRJNA541421_bin_98_orig_refined-contigs","PRJNA531756_bin_72_orig_refined-contigs",
                                                                    "GCA_021160765.1_ASM2116076v1_genomic","GCA_021158585.1_ASM2115858v1_genomic",
                                                                    "GCA_003663625.1_ASM366362v1_genomic","PRJNA368391_bin_90_orig-contigs","GCA_003650025.1_ASM365002v1_genomic","GCA_011039765.1_ASM1103976v1_genomic",
                                                                    "ZHEN22-1_SAMN22703512_MAG_00000312","ZHEN22-1_SAMN22703514_MAG_00000188","ZHEN22-1_SAMN22703513_MAG_00000265",
                                                                    "ZORZ22-1_SAMN30647027_MAG_00000060","GCA_023145185.1_ASM2314518v1_genomic", 
                                                                    "PRJNA531756_bin_93_strict_refined-contigs","GCA_016928095.1_ASM1692809v1_genomic","PRJNA889212_bin_1_orig_refined-contigs", "PRJNA541421_bin_72_orig_refined-contigs","PRJNA721298_bin_46_orig_refined-contigs",
                                                                    "EMB267_Co_bin_434_strict_1","PRJNA368391_bin_183_orig-contigs","GCA_003663595.1_ASM366359v1_genomic",
                                                                    "MSM105_N25025F_bin_277_ori_permissive_1","PRJNA541421_bin_70_strict-contigs","PRJNA889212_bin_7_strict_refined-contigs","E3_1_d157_spades_bin_16_orig_1_refined-contigs"),
                                 labels = c("SAMN10231914_1", "SAMN10231894_1","SAMN10231913_1", "SAMN10231893_1",
                                            "SAMN10231914_2","SAMN10231913_2", 
                                            "SAMN10231904_1","SAMN11854494_1","PRJNA704804_1","SAMEA2620113_1",
                                            "PRJNA541421_3","PRJNA531756_1",
                                            "GCA_021160765.1","GCA_021158585.1",
                                            "GCA_003663625.1","PRJNA368391_2","GCA_003650025.1","GCA_011039765.1",
                                            "SAMN22703512_1","SAMN22703514_1","SAMN22703513_1",
                                            "SAMN30647027_1","GCA_023145185.1",
                                            "PRJNA531756_2","GCA_016928095.1","PRJNA889212_1","PRJNA541421_2","PRJNA721298_1",
                                            "EMB267","PRJNA368391_1","GCA_003663595.1","MSM105","PRJNA541421_1","PRJNA889212_2","E3_1_d157"),ordered =T)

matrix_TCDB_melt_short$genome <- factor(matrix_TCDB_melt_short$genome, levels = c("SUTE22-1_SAMN10231914_MAG_00000186","SUTE22-1_SAMN10231894_MAG_00000173","SUTE22-1_SAMN10231913_MAG_00000219","SUTE22-1_SAMN10231893_MAG_00000067",
                                                                      "SUTE22-1_SAMN10231914_MAG_00000292","SUTE22-1_SAMN10231913_MAG_00000177",
                                                                      "SUTE22-1_SAMN10231904_MAG_00000012","SCRA20-1_SAMN11854494_MAG_00000079","PRJNA704804_bin_147_orig_refined-contigs","TARA_SAMEA2620113_MAG_00000097",
                                                                      "PRJNA541421_bin_98_orig_refined-contigs","PRJNA531756_bin_72_orig_refined-contigs",
                                                                      "GCA_021160765.1_ASM2116076v1_genomic","GCA_021158585.1_ASM2115858v1_genomic",
                                                                      "GCA_003663625.1_ASM366362v1_genomic","PRJNA368391_bin_90_orig-contigs","GCA_003650025.1_ASM365002v1_genomic","GCA_011039765.1_ASM1103976v1_genomic",
                                                                      "ZHEN22-1_SAMN22703512_MAG_00000312","ZHEN22-1_SAMN22703514_MAG_00000188","ZHEN22-1_SAMN22703513_MAG_00000265",
                                                                      "ZORZ22-1_SAMN30647027_MAG_00000060","GCA_023145185.1_ASM2314518v1_genomic", 
                                                                      "PRJNA531756_bin_93_strict_refined-contigs","GCA_016928095.1_ASM1692809v1_genomic","PRJNA889212_bin_1_orig_refined-contigs", "PRJNA541421_bin_72_orig_refined-contigs","PRJNA721298_bin_46_orig_refined-contigs",
                                                                      "EMB267_Co_bin_434_strict_1","PRJNA368391_bin_183_orig-contigs","GCA_003663595.1_ASM366359v1_genomic",
                                                                      "MSM105_N25025F_bin_277_ori_permissive_1","PRJNA541421_bin_70_strict-contigs","PRJNA889212_bin_7_strict_refined-contigs","E3_1_d157_spades_bin_16_orig_1_refined-contigs"),
                                  labels = c("SAMN10231914_1", "SAMN10231894_1","SAMN10231913_1", "SAMN10231893_1",
                                             "SAMN10231914_2","SAMN10231913_2", 
                                             "SAMN10231904_1","SAMN11854494_1","PRJNA704804_1","SAMEA2620113_1",
                                             "PRJNA541421_3","PRJNA531756_1",
                                             "GCA_021160765.1","GCA_021158585.1",
                                             "GCA_003663625.1","PRJNA368391_2","GCA_003650025.1","GCA_011039765.1",
                                             "SAMN22703512_1","SAMN22703514_1","SAMN22703513_1",
                                             "SAMN30647027_1","GCA_023145185.1",
                                             "PRJNA531756_2","GCA_016928095.1","PRJNA889212_1","PRJNA541421_2","PRJNA721298_1",
                                             "EMB267","PRJNA368391_1","GCA_003663595.1","MSM105","PRJNA541421_1","PRJNA889212_2","E3_1_d157"),ordered =T)


## plot heatmap for all TCDB families ####
ggplot(matrix_TCDB_melt_short, aes(x = genome, y = TCDB_family)) +
  geom_point(aes(size=value, fill=value, col=value), shape = 21) +
  theme_bw() +
  coord_flip() +
  scale_y_discrete(position="right")+
  ggtitle("TCDB family overview")+
  xlab("") + ylab("")+
  theme(plot.title = element_text(size=7, face = "bold"),
        axis.text.x = element_text(angle=-90, hjust=1, vjust=-0.8, color='grey40', size=7),
        axis.text.y =element_text(color='grey40',size=7),
        legend.position = 'right',
        legend.title = element_text(size=7),
        legend.text = element_text(size=7),
        axis.title=element_text(face="bold", size=7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_color_viridis_c("no. of genes", guide = "legend", option="mako", direction=-1,limits = c(0.1,20), breaks = c(0,2,4,6,8,10,12,14,16))+
  scale_fill_viridis_c("no. of genes",guide = "legend", option="mako", direction=-1,limits = c(0.1,20), breaks = c(0,2,4,6,8,10,12,14,16))+
  scale_size_area("no. of genes",max_size = 8,limits = c(0.1,20), breaks = c(0,2,4,6,8,10,12,14,16))

ggsave("TCDB/TCDB_families_overview.png", plot = last_plot(), width = 25, height = 30, units = "cm", device = "png" )
