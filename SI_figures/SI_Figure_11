#Script to combine annotation for peptidase and extracellular peptidase heatmaps

library(tidyverse)
library(reshape2)

setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Data_mining/Annotation_v3")

# this file was created during Data mining annotation of EX MAGs
annotation <- read.table("annotation_parsed.txt", 
                         sep="\t", 
                         header=TRUE, 
                         fill=TRUE, 
                         quote="")

#count different peptidase families within MAGs
Peptidases <- aggregate(annotation$merops_sseqid ~ annotation$genome, data = annotation, FUN = function(x) sum(!is.na(x)))
print(Peptidases)

# Peptidase heatmap ----
## kegg ----
kegg <- read.csv("full/concat_kegg_annotation_best.txt", 
                 sep="\t", 
                 header=TRUE, 
                 fill=TRUE, 
                 quote="",
                 col.names=c("qseqid", "kegg_sseqid","kegg_pident","kegg_length","kegg_mismatch","kegg_gapopen","kegg_qstart",
                             "kegg_qend","kegg_sstart","kegg_send","kegg_evalue","kegg_bitscore","kegg_bsr","kegg_KO","kegg_KO_name","kegg_taxon","kegg_gene"))
kegg <- subset(kegg, kegg_sseqid != "sseqid")
kegg$qseqid <- gsub(".genes","",as.character(kegg$qseqid))

## NR ----
NR <- read.csv("full/concat_NR_annotation_best.txt", 
               sep="\t", 
               header=TRUE, 
               fill=TRUE, 
               quote="",
               col.names=c("qseqid", "NR_sseqid","NR_pident","NR_length","NR_mismatch","NR_gapopen","NR_qstart",
                           "NR_qend","NR_sstart","NR_send","NR_evalue","NR_bitscore","NR_qlen","NR_slen","NR_stitle","NR_staxids","NR_sscinames","NR_salltitles","NR_bsr"))
NR <- subset(NR, NR_sseqid != "sseqid")
NR$qseqid <- gsub(".genes","",as.character(NR$qseqid))


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

# combine table with Peptidases in first few columns
Peptidases_table <- full_join(merops, kegg, by = c("qseqid")) %>%
  full_join(., NR, by = c("qseqid")) %>%
  separate(col = qseqid, into = c("genome","cluster"), sep = "___")

write.table(Peptidases_table, "Peptidases/Peptidases_table_parsed.txt", 
            quote = FALSE,
            sep="\t",
            col.names = TRUE,
            row.names = FALSE)

# remove all rows from dataframe which don't have any annotation through NR or kegg
Peptidases_short <- Peptidases_table[complete.cases(Peptidases_table[, 3:19]) & complete.cases(Peptidases_table[, 36:53]), ]
Peptidases_short$genome <- gsub(".genes","",as.character(Peptidases_short$genome))
# 634 Peptidases remain

write.table(Peptidases_short, "Peptidases/Peptidases_short_parsed.txt", 
            quote = FALSE,
            sep="\t",
            col.names = TRUE,
            row.names = FALSE)

# Create a matrix of Merops Peptidase families against genomes
matrix_Pep <- table(Peptidases_short$genome, Peptidases_short$merops_subfamily)
# melt matrix into tabular view
matrix_Pep_melt <- melt(matrix_Pep)
colnames(matrix_Pep_melt) <- c("genome","merops_subfamily","value")

matrix_Pep_melt$genome <- factor(matrix_Pep_melt$genome, levels = c("SUTE22-1_SAMN10231914_MAG_00000186","SUTE22-1_SAMN10231894_MAG_00000173","SUTE22-1_SAMN10231913_MAG_00000219","SUTE22-1_SAMN10231893_MAG_00000067",
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

# plot heatmap
Pep_plot <- ggplot(matrix_Pep_melt[], aes(x = genome, y = merops_subfamily)) +
  geom_point(aes(size=value, fill=value, col=value), shape = 21)+
  theme_bw() +
  coord_flip() +
  scale_y_discrete(position="right")+
  xlab("genome") + ylab("merops subfamily")+
  theme(axis.text.x = element_text(angle=-90, hjust=1, vjust=-0.8, color='grey40', size=7),
        axis.text.y =element_blank(),
        legend.position = 'right',
        legend.title = element_text(size=7, face = "bold"),
        legend.text = element_text(size=7),
        axis.title=element_text(face="bold", size=7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_color_viridis_c("peptidase counts", guide = "legend", option="mako", direction=-1)+
  scale_fill_viridis_c("peptidase counts", guide = "legend", option="mako", direction=-1)+
  scale_size_area("peptidase counts",max_size = 5)

ggsave("Peptidases/Peptidase_families_overview.png", plot = last_plot(), width = 18, height = 21, units = "cm", device = "png" )
ggsave("Peptidases/Peptidase_families_overview.pdf", plot = last_plot(), width = 18, height = 21, units = "cm", device = "pdf" )


# Signal P Heatmap ----
# Search different signal peptides 
# (signal peptides are short amino acid sequences, which control protein secretion and translocation)

#rename ID in signalp table
signalp$ID <- gsub(" #.*","", signalp$ID ) 

#combine signalp and merops tables
SPs <- full_join(signalp, merops, by = c("ID"="qseqid" ))

# filter only those rows in which there is data for signalp and merops predicition
SPs_prediction <- SPs[complete.cases(SPs$merops_id), ]
unique(SPs_prediction$Prediction)
# "OTHER" "SP" , so there is either no signal peptide or signal peptides classified as SP (Sec/SPI)

# reduce data set
SPs_prediction <- separate(SPs_prediction, ID, into = c("genome","cluster"), sep = "___")



# Summarize the data to get the count of SP predictions
summarized_data <- SPs_prediction %>%
  filter(Prediction == "SP") %>%
  group_by(genome, merops_subfamily) %>%
  summarise(count = n()) %>%
  ungroup()

# Identify merops_subfamilies with SP counts
relevant_merops_subfamilies <- summarized_data %>%
  filter(count > 0) %>%
  pull(merops_subfamily) %>%
  unique()

# Get a complete list of genomes and relevant merops_subfamilies
all_genomes <- unique(SPs_prediction$genome)
complete_data <- expand.grid(genome = all_genomes, merops_subfamily = relevant_merops_subfamilies)

# Merge the complete data with the summarized data
merged_data <- complete_data %>%
  left_join(summarized_data, by = c("genome", "merops_subfamily")) %>%
  replace_na(list(count = 0))

merged_data$genome <- factor(merged_data$genome, levels = c("SUTE22-1_SAMN10231914_MAG_00000186","SUTE22-1_SAMN10231894_MAG_00000173","SUTE22-1_SAMN10231913_MAG_00000219","SUTE22-1_SAMN10231893_MAG_00000067",
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



SP_plot <- ggplot(merged_data, aes(x = genome, y = merops_subfamily)) +
  geom_point(aes(size=count, fill=count, col=count), shape = 21)+
  theme_bw() +
  coord_flip() +
  scale_y_discrete(position="right")+
  xlab("genome") + ylab("merops subfamily")+
  theme(axis.text.x = element_text(angle=-90, hjust=1, vjust=-0.8, color='grey40', size=7),
        axis.text.y =element_text(color='grey40',size=7),
        legend.position = 'right',
        legend.title = element_text(size=7, face = "bold"),
        legend.text = element_text(size=7),
        axis.title=element_text(face="bold", size=7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_color_viridis_c("signal peptidase \ncounts",guide = "legend", option="rocket", direction=-1,breaks = c(0, 1, 2))+
  scale_fill_viridis_c("signal peptidase \ncounts",guide = "legend", option="rocket", direction=-1, breaks = c(0, 1, 2))+
  scale_size_area("signal peptidase \ncounts",max_size = 3,
                  breaks=c(0,1,2),
                  labels = c("0","1","2"))

ggsave("Peptidases/Peptidase_SPs_families_overview.png", plot = last_plot(), width = 12, height = 21, units = "cm", device = "png" )
ggsave("Peptidases/Peptidase_SPs_families_overview.pdf", plot = last_plot(), width = 12, height = 21, units = "cm", device = "pdf" )


library(patchwork)

(SP_plot | Pep_plot ) + 
  plot_layout(guides="collect", axis_titles = "collect", widths = c(3,30)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 7, face = "bold"))

ggsave("Peptidases/Pep_all.png", plot = last_plot(), width = 18, height = 21, units = "cm", device = "png" )
ggsave("Peptidases/Pep_all.pdf", plot = last_plot(), width = 18, height = 21, units = "cm", device = "pdf" )


