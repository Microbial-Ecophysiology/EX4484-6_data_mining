
# parse annotation

library(tidyverse)
library(reshape)
library(patchwork)
library(RColorBrewer)
library(vegan)
library(ImpactEffectsize)

setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Data_mining/Annotation_v3/hypothetical_allThermo")


load("hypothetical_unknown_OG_analysis.RData")

# load annotation ####

## kegg ----
kegg <- read.csv("concat_kegg_annotation_best.txt", 
                 sep="\t", 
                 header=TRUE, 
                 fill=TRUE, 
                 quote="",
                 col.names=c("qseqid", "kegg_sseqid","kegg_pident","kegg_length","kegg_mismatch","kegg_gapopen","kegg_qstart",
                             "kegg_qend","kegg_sstart","kegg_send","kegg_evalue","kegg_bitscore","kegg_bsr","kegg_KO","kegg_KO_name","kegg_taxon","kegg_gene"))
kegg <- subset(kegg, kegg_sseqid != "sseqid")

kegg_EX <- read.csv("../full/concat_kegg_annotation_best.txt", 
                 sep="\t", 
                 header=TRUE, 
                 fill=TRUE, 
                 quote="",
                 col.names=c("qseqid", "kegg_sseqid","kegg_pident","kegg_length","kegg_mismatch","kegg_gapopen","kegg_qstart",
                             "kegg_qend","kegg_sstart","kegg_send","kegg_evalue","kegg_bitscore","kegg_bsr","kegg_KO","kegg_KO_name","kegg_taxon","kegg_gene"))
kegg_EX <- subset(kegg_EX, kegg_sseqid != "sseqid")
# kegg_EX$qseqid <- gsub(".genes", "", kegg_EX$qseqid)
# kegg_EX$qseqid <- gsub("refined-contigs", "refined.contigs", kegg_EX$qseqid)


## NR ----
NR <- read.csv("concat_NR_annotation_best.txt", 
               sep="\t", 
               header=TRUE, 
               fill=TRUE, 
               quote="",
               col.names=c("qseqid", "NR_sseqid","NR_pident","NR_length","NR_mismatch","NR_gapopen","NR_qstart",
                           "NR_qend","NR_sstart","NR_send","NR_evalue","NR_bitscore","NR_qlen","NR_slen","NR_stitle","NR_staxids","NR_sscinames","NR_salltitles","NR_bsr"))
NR <- subset(NR, NR_sseqid != "sseqid")

NR_EX <- read.csv("../full/concat_NR_annotation_best.txt", 
               sep="\t", 
               header=TRUE, 
               fill=TRUE, 
               quote="",
               col.names=c("qseqid", "NR_sseqid","NR_pident","NR_length","NR_mismatch","NR_gapopen","NR_qstart",
                           "NR_qend","NR_sstart","NR_send","NR_evalue","NR_bitscore","NR_qlen","NR_slen","NR_stitle","NR_staxids","NR_sscinames","NR_salltitles","NR_bsr"))
NR_EX <- subset(NR_EX, NR_sseqid != "sseqid")
# NR_EX$qseqid <- gsub(".genes", "", NR_EX$qseqid)
# NR_EX$qseqid <- gsub("refined-contigs", "refined.contigs", NR_EX$qseqid)


# Combine annotation data of kegg and NR
annotation <- full_join(kegg, NR, by = c("qseqid")) 
annotation_EX <- full_join(kegg_EX, NR_EX, by = c("qseqid")) 
all.equal(colnames(annotation), colnames(annotation_EX))
annotation_all <- rbind(annotation, annotation_EX)

# remove qseqid's of outgroup taxa
# Korarchaeia, Methanomethylicia, Methanosarcinia, Nitrososphaeria, Bathy and 2 EX MAGs, which are duplicated in the data set

# GCA_001918745.1_ASM191874v1_genomic, GCA_001593935.1_ASM159393v1_genomic, GCA_018396415.1_ASM1839641v1_genomic, GCA_003601775.1_ASM360177v1_genomic
# GCF_000019605.1_ASM1960v1_genomic, GCA_001717035.1_ASM171703v1_genomic, GCA_000200715.1_ASM20071v1_genomic,
# GCF_000007345.1_ASM734v1_genomic, GCF_000970085.1_ASM97008v1_genomic, GCF_000970205.1_ASM97020v1_genomic,GCF_000970265.1_ASM97026v1_genomic,
# GCF_000970285.1_ASM97028v1_genomic, GCF_001304615.2_ASM130461v2_genomic, GCF_002287235.1_ASM228723v1_genomic, 
# PRJNA713414_bin_118_orig_1_refined-contigs,PRJNA713414_bin_178_orig_1_refined-contigs

accession_numbers_to_exclude <- c(
  "GCA_001918745.1_ASM191874v1_genomic", "GCA_001593935.1_ASM159393v1_genomic", 
  "GCA_018396415.1_ASM1839641v1_genomic", "GCA_003601775.1_ASM360177v1_genomic",
  "GCF_000019605.1_ASM1960v1_genomic", "GCA_001717035.1_ASM171703v1_genomic", 
  "GCA_000200715.1_ASM20071v1_genomic", "GCF_000007345.1_ASM734v1_genomic", 
  "GCF_000970085.1_ASM97008v1_genomic", "GCF_000970205.1_ASM97020v1_genomic", 
  "GCF_000970265.1_ASM97026v1_genomic", "GCF_000970285.1_ASM97028v1_genomic", 
  "GCF_001304615.2_ASM130461v2_genomic", "GCF_002287235.1_ASM228723v1_genomic",
  "PRJNA713414_bin_118_orig_1_refined-contigs","PRJNA713414_bin_178_orig_1_refined-contigs"
)

# Create a pattern that matches any qseqid starting with any of the accession numbers
pattern_exclude <- paste0("^(", paste(accession_numbers_to_exclude, collapse = "|"), ")")

# Filter out rows with qseqid matching the pattern
annotation_filtered <- annotation_all %>%
  filter(!grepl(pattern_exclude, qseqid))

annotation_thermo <- annotation
annotation <- annotation_all

annotation$kegg_KO[is.na(annotation$kegg_KO)] <- ""
annotation$kegg_gene[is.na(annotation$kegg_gene)] <- ""
annotation$NR_stitle[is.na(annotation$NR_stitle)] <- ""

uniq_genes <- unique(annotation$kegg_gene[!is.na(annotation$kegg_KO)])
uniq_stitle <- unique(annotation$NR_stitle)
grep("hypothetical|[pP]rotein of unknown function", uniq_genes, value = T)
head(grep("hypothetical|[pP]rotein of unknown function", uniq_stitle, value = T))


tmp <- data.frame(
  tmp_KO = ifelse(annotation$kegg_KO != "", "annotated", ""),
  tmp_gene = ifelse(!grepl("hypothetical|[pP]rotein of unknown function", annotation$kegg_gene) & annotation$kegg_gene != "", "annotated", ""),
  tmp_stitle = ifelse(!grepl("hypothetical|[pP]rotein of unknown function", annotation$NR_stitle) & annotation$NR_stitle != "", "annotated", "")
)
annotation$annotation_status <- ifelse(apply(tmp, 1, function(x) sum(x == "") == 3), "hypothetical", "annotated")
tmp <- annotation[, c("qseqid", "annotation_status")]

# load Orthogroups.txt
orthogroups <- scan("../../Orthofinder_allThermo/Orthogroups.txt", what = "character", sep = "\n")
length(orthogroups)
tmp_df <- data.frame(
  do.call(
    "rbind",
    strsplit(orthogroups, ": ")
  )
)
tmp_list <- strsplit(tmp_df$X2, " ", fixed = T)
tmp_df_filt <- tmp_df[sapply(tmp_list, length) > 1, ]
tmp_list_filt <- tmp_list[sapply(tmp_list, length) > 1]
orthogroups_map <- data.frame(
  OG = rep(tmp_df_filt$X1, sapply(tmp_list_filt, length)),
  qseqid = unlist(tmp_list_filt)
)
table(table(orthogroups_map$OG))
orthogroups_annotation_status <- left_join(orthogroups_map, tmp, by = c("qseqid")) 
orthogroups_annotation_status$annotation_status[is.na(orthogroups_annotation_status$annotation_status)] <- "unknown"
orthogroups_annotation_status$annotation_status <- as.factor(orthogroups_annotation_status$annotation_status)
table(orthogroups_annotation_status$annotation_status)
OG_annotation_status <- data.frame(do.call("rbind", by(orthogroups_annotation_status$annotation_status, orthogroups_annotation_status$OG, table)))
dim(OG_annotation_status)
OG_annotation_status$OG_status <- ifelse(OG_annotation_status$annotated > 0, "annotated", ifelse(OG_annotation_status$hypothetical > 0, "hypothetical", "unknown"))
table(OG_annotation_status$OG_status)

OG_counts <- read.table("../../Orthofinder_allThermo/Orthogroups.GeneCount.tsv", h = T, sep = "\t", check.names = F)
OG_counts$OG_status <- OG_annotation_status[OG_counts$Orthogroup, "OG_status"]

# remove outgroup MAGs as above 
OG_counts_filtered <- OG_counts %>%
  select(-all_of(accession_numbers_to_exclude)) # excluding outgroup MAGs with accession_numbers_to_exclude

rm(OG_counts)

# filter rows only present in EX genomes
# vegan decostand pa

OG_counts_pa <- decostand(OG_counts_filtered[, 2:419], "pa")

# rowSums will give you number of genomes present
OG_counts_filtered$n_genomes_shared <- rowSums(OG_counts_pa)

# boxplot over OG annotation status
plot_all <- boxplot(OG_counts_filtered$n_genomes_shared ~OG_counts_filtered$OG_status)

# filter rows only present in EX genomes
EX_MAGs_unique_tmp <- annotation_EX %>% 
  separate(col = qseqid, into = c("genome","cluster"), sep = "___")
EX_MAGs_unique <- unique(EX_MAGs_unique_tmp$genome)
# remove PRJNA713414 from list
remove <- c("PRJNA713414_bin_118_orig_1_refined-contigs","PRJNA713414_bin_178_orig_1_refined-contigs")
EX_MAGs_unique <- setdiff(EX_MAGs_unique, remove)
EX_MAGs_unique <- gsub(".genes$", "", EX_MAGs_unique)

# Create a logical vector indicating rows with counts in EX_MAGs_unique
has_counts_in_EX_MAGs <- rowSums(OG_counts_filtered[, EX_MAGs_unique]) > 0 # 6231 OG in total

# Create a logical vector indicating rows with counts in other MAGs
has_counts_in_other_MAGs <- rowSums(OG_counts_filtered[, c(2:419)[!(colnames(OG_counts_filtered)[2:419] %in% EX_MAGs_unique)]]) > 0

## Orthogroups of EX ####
subset_og_counts_EX_all <- OG_counts_filtered[has_counts_in_EX_MAGs, ]
subset_og_counts_EX_all$n_shared_EX_only <- rowSums(subset_og_counts_EX_all[, EX_MAGs_unique] > 0)
subset_og_counts_EX_all$n_shared_EX_only_perc <- subset_og_counts_EX_all$n_shared_EX_only / length(EX_MAGs_unique) * 100
boxplot(subset_og_counts_EX_all$n_shared_EX_only ~ subset_og_counts_EX_all$OG_status)
boxplot(subset_og_counts_EX_all$n_shared_EX_only_perc ~ subset_og_counts_EX_all$OG_status)

## presence of EX orthogroups in non-EX MAGs ####
subset_og_counts_EX_all$n_shared_nonEX_only <- rowSums(subset_og_counts_EX_all[, c(2:419)[!(colnames(OG_counts_filtered)[2:435] %in% EX_MAGs_unique)]] > 0)
subset_og_counts_EX_all$n_shared_nonEX_only_perc <- subset_og_counts_EX_all$n_shared_nonEX_only / (418 - length(EX_MAGs_unique)) * 100
boxplot(subset_og_counts_EX_all$n_shared_nonEX_only ~ subset_og_counts_EX_all$OG_status)
boxplot(subset_og_counts_EX_all$n_shared_nonEX_only_perc ~ subset_og_counts_EX_all$OG_status)
sum(subset_og_counts_EX_all$n_shared_nonEX_only > 0 & subset_og_counts_EX_all$OG_status == "unknown") #88
sum(subset_og_counts_EX_all$n_shared_nonEX_only > 0 & subset_og_counts_EX_all$OG_status == "hypothetical") #603


## boxplots ####

df_boxplot <- subset_og_counts_EX_all %>% 
  select(Orthogroup, n_shared_EX_only_perc, n_shared_nonEX_only_perc, OG_status )

df_boxplot_melt <- melt(df_boxplot)

install.packages("ggpubr")
library(ggpubr)


ggplot()+
  geom_boxplot(data = df_boxplot_melt, 
               aes(x = OG_status, y =value, fill = variable))+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size= 7),
    axis.title = element_text(size= 7, face = "bold"),
    legend.title=element_text(size= 7, face = "bold"),
    legend.text = element_text(size= 7))+
  xlab("")+
  ylab("Thermoplasmatota containing OGs found in EX4484-6 MAGs [%]")+ 
  scale_fill_manual("shared orthogroups", values = c("#ef8a62", "#2166ac"),
                    labels = c("EX4484-6", "other Thermoplasmatota"))+ 
  stat_compare_means(data =df_boxplot_melt, aes(x = OG_status, y =value, group=variable),
                     label = "p.signif", method = "wilcox.test")+
  scale_y_continuous(limits = c(0,105))

ggsave("boxplot_shared_OGs.png", width = 18, height = 12, units = "cm")
ggsave("boxplot_shared_OGs.pdf", width = 18, height = 12, units = "cm")
ggsave("boxplot_shared_OGs.emf", width = 18, height = 12, units = "cm")


# test for differences between orthogroups
results_test <- lapply(
  unique(df_boxplot_melt$OG_status),
  function(x) {
    wilcox.test(value ~ variable, data = df_boxplot_melt[df_boxplot_melt$OG_status == x, ], conf.int = TRUE)
  }
)
names(results_test) <- unique(df_boxplot_melt$OG_status)

# add impact effectsize manually to plots
# Impact Effectsize to calculate the Effectsize for each of the OG status ####

# Impact = the main effect size measure
# MorphDiff = the extend of the group-difference in the shapes of the pdf
# CTDiff =  the extend of the difference in the group medians

## annotated ####
Impact(Data = (df_boxplot_melt$value[df_boxplot_melt$OG_status == "annotated"]), 
       Cls = (df_boxplot_melt$variable[df_boxplot_melt$OG_status == "annotated"]))
# $Impact
# [1] -0.3495771
# 
# $MorphDiff
# [1] 0.3526224
# 
# $CTDiff
# [1] 0.3344095

## hypothetical ####
Impact(Data = (df_boxplot_melt$value[df_boxplot_melt$OG_status == "hypothetical"]), 
       Cls = (df_boxplot_melt$variable[df_boxplot_melt$OG_status == "hypothetical"]))
# $Impact
# [1] -1.195459
# 
# $MorphDiff
# [1] 1.444463
# 
# $CTDiff
# [1] 0.8758858

## unknown ####
Impact(Data = (df_boxplot_melt$value[df_boxplot_melt$OG_status == "unknown"]), 
       Cls = (df_boxplot_melt$variable[df_boxplot_melt$OG_status == "unknown"]))
# $Impact
# [1] -2.26678
# 
# $MorphDiff
# [1] 1.822755
# 
# $CTDiff
# [1] 2.26678



# check OG abundance for each class to see if some classes stand out ####
# subset OG
tmp_og <- data.frame(OG_counts_filtered$Orthogroup, OG_counts_pa, check.names = FALSE)[has_counts_in_EX_MAGs, ]

tmp_og_melt <- melt(tmp_og) %>% 
  filter(value > 0)

# add taxonomy from gtdb
gtdb <- read.table("GTDB_classification_mod.txt", h = T, sep = "\t")

#check, if all accession numbers are present
sum(unique(tmp_og_melt$variable) %in% rownames(gtdb))

# melt dataframe
tmp_og_melt_2 <- left_join(tmp_og_melt, 
          gtdb,
          by = join_by("variable" == "user_genome")) %>% 
  select(-domain, -phylum, -order)

# calculate in how many MAGs of class each OG is present
tmp_og_cast <- cast(tmp_og_melt_2, "OG_counts_filtered$Orthogroup ~ class", value = "value", fun.aggregate = sum)

# filter gtdb table and remove outgroup lineages:
gtdb_filtered <- gtdb %>%
  filter(!user_genome %in% accession_numbers_to_exclude) # excluding outgroup MAGs with accession_numbers_to_exclude

# calculate percentage of presence for each of the classes
tmp_og_cast_perc <- apply(
  tmp_og_cast %>% column_to_rownames("OG_counts_filtered$Orthogroup"),
  1,
  function(x) {
    x / table(gtdb_filtered$class) * 100
  }
) %>% t()


# add annotation status
tmp_og_cast_perc <- merge(tmp_og_cast_perc,
                  OG_annotation_status,
                  by = 0) %>%
  select(-annotated, -hypothetical, -unknown)

rownames(tmp_og_cast_perc) <- tmp_og_cast_perc[,1]
tmp_og_cast_perc <- tmp_og_cast_perc[,-1]

tmp_og_cast_perc_melt <- melt(tmp_og_cast_perc)


# plot for all classes
ggplot()+
  geom_boxplot(data = tmp_og_cast_perc_melt, 
               aes(x = OG_status, y =value, fill = variable))+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size= 7),
    axis.title = element_text(size= 7, face = "bold"),
    legend.title=element_text(size= 7, face = "bold"),
    legend.text = element_text(size= 7))+
  xlab("")+
  ylab("Thermoplasmatota containing OGs found in EX4484-6 MAGs [%]")+
  scale_fill_manual("class", values = c("#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf", "#2166ac"))

ggsave("boxplot_shared_OGs_per_class.png", width = 18, height = 12, units = "cm")

# based on these results all classes except EX4484-6 show similar results, 
# therefore going with boxplot combining all other Thermoplasmatota classes as nonEX


save.image("hypothetical_unknown_OG_analysis.RData")



# hexplots ####

## for orthogroups with annotation "hypothetical" ####
ggplot()+
  geom_hex(data = subset_og_counts_EX_all[subset_og_counts_EX_all$OG_status == "hypothetical",], 
           aes(x = n_shared_EX_only_perc, y =n_shared_nonEX_only_perc),color = "white")+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size= 7),
    axis.title = element_text(size= 7, face = "bold"),
    legend.title = element_text(size= 7, face = "bold"),
    legend.text = element_text(size= 7))+
  xlab("shared orthogroups within EX4484-6 [%]")+
  ylab("shared orthogroups within other Thermoplasmatota [%]")+
  scale_fill_gradientn("hypothetical \northogroups" ,colours=rev(brewer.pal(10,"RdGy")),
                       limits=c(0,165),
                       breaks = c(0,40,80,120,160))+
  guides(fill = guide_colourbar(barwidth = 1,
                                barheight = 10))

ggsave("hexplot_hypothetical_OG.png", width = 14, height = 12, units = "cm")

## for orthogroups with annotation "unknown" ####
ggplot()+
  geom_hex(data = subset_og_counts_EX_all[subset_og_counts_EX_all$OG_status == "unknown",], 
           aes(x = n_shared_EX_only_perc, y =n_shared_nonEX_only_perc),color = "white",bins = 20)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size= 7),
    axis.title = element_text(size= 7, face = "bold"),
    legend.title = element_text(size= 7, face = "bold"),
    legend.text = element_text(size= 7))+
  xlab("shared orthogroups within EX4484-6 [%]")+
  ylab("shared orthogroups within other Thermoplasmatota [%]")+
  scale_fill_gradientn("unknown \northogroups", colours=rev(brewer.pal(10,"RdGy")),
                       limits=c(0,1200),
                       breaks = c(0,250,500,750,1100))+
  guides(fill = guide_colourbar(barwidth = 1,
                                barheight = 10))

ggsave("hexplot_unknown_OG.png", width = 14, height = 12, units = "cm")


## for orthogroups with annotation "unknown" ####
ggplot()+
  geom_hex(data = subset_og_counts_EX_all[subset_og_counts_EX_all$OG_status == "annotated",], 
           aes(x = n_shared_EX_only_perc, y =n_shared_nonEX_only_perc),color = "white")+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size= 7),
    axis.title = element_text(size= 7, face = "bold"),
    legend.title = element_text(size= 7, face = "bold"),
    legend.text = element_text(size= 7))+
  xlab("shared orthogroups within EX4484-6 [%]")+
  ylab("shared orthogroups within other Thermoplasmatota [%]")+
  scale_fill_gradientn("annotated \northogroups", colours=rev(brewer.pal(10,"RdGy")),
                       limits=c(0,130),
                       breaks = c(0,25,50,75,100,125))+
  guides(fill = guide_colourbar(barwidth = 1,
                                barheight = 10))

ggsave("hexplot_annotated_OG.png", width = 14, height = 12, units = "cm")
