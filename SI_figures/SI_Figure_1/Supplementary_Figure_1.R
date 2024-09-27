## script for subsetting a taxmap object to a final long format table with taxa above a certain threshold

# packages ####
library(tidyverse)
library(taxa)
library(metacoder)
library(phyloseq)
library(hues)
library(cowplot)
library(devEMF)


# import data ####
setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Nadja_BSc/16S_Amplicon/dada2_analysis_new_14102023/Plot_rel_abundance_Archaea/")

## ASV table raw counts
asv <- read_tsv("../Archaea_all/results/asv_table_tax_Lib28_arc_GTDB_mod.txt")

## metadata, need to have same sample names as in taxmap
mdata <- read_tsv("../meta_data_mod.txt")

# create taxmap object and do some quality control ####

## get sample names
sample.names <- mdata$Sample

## make taxmap object ####
obj <- parse_tax_data(
  tax_data = asv,            # which table to use
  class_cols = c("Kingdom", "Phylum", "Class", "Order", "Family","Genus"),   # which columns have the taxonomy
  named_by_rank = T          # keep the column names as rank names in command 'taxon_ranks'
)
obj
# the object has some additional columns now with the ASV ID and the taxonomy

### remove taxonomy columns from ASV count table, only keep ASV_ID and taxon_ID column
obj$data$tax_data <- obj$data$tax_data[c("taxon_id", "ASV_ID", sample.names)]

### rename table to asv_table
names(obj$data) <- "asv_table"

obj
# 367 taxa (including subtaxa)
# 1,142 ASVs (=number of rows of tax_data table)
ini_reads <- sum(obj$data$asv_table[, sample.names])
ini_reads
# 376132 reads in original table

## modify taxmap object ####
### remove any ASVs which are not Archaea
obj <- metacoder::filter_taxa(obj,
                              taxon_names == "Archaea",
                              subtaxa = T)
sum(obj$data$asv_table[, sample.names])/ini_reads
#  1

### remove doubletons and singletons
obj$data$asv_table <- metacoder::zero_low_counts(obj, "asv_table",
                                                 min_count = 3,
                                                 use_total = T,
                                                 other_cols = T)
#Zeroing 147 of 1142 rows with total counts less than 3

### check for empty ASVs
no_reads <- rowSums(obj$data$asv_table[, sample.names]) == 0
sum(no_reads)  # 450 empty ASVs

# remove empty ASVs
obj <- metacoder::filter_obs(obj, "asv_table",
                             ! no_reads,
                             drop_taxa = T)
obj
# taxa reduced to 268
# ASVs reduced to 692
sum(obj$data$asv_table[, sample.names])/ini_reads
# 99.94% of the reads kept

# check taxonomy table if everything looks ok
print(taxonomy_table(obj), n = 300)

## calculate further tables ####
## Calculate relative abundance
# relative abundance per ASV
obj$data$rel_abd <- calc_obs_props(obj, "asv_table", other_cols = T)
# relative abundance per taxon
obj$data$rel_tax_abd <- calc_taxon_abund(obj, "rel_abd")
print(obj)

## save or load created taxmap object ####
save(obj, file = "taxmap_MAC_arc_GTDB_mod.RData")

load("taxmap_MAC_arc_GTDB_mod.RData")
# here named 'obj'


# create phyloseq object ####
# maybe need to change names of otu_table, otu_id_col and sample names (= column name of column with sample names)
objphy <- as_phyloseq(obj, otu_table = "asv_table", otu_id_col = "ASV_ID", sample_data = mdata, sample_id_col = "Sample")

ntaxa(objphy)                  # how many taxa (including subtaxa) - 692
nsamples(objphy)               # how many samples - 8
sample_names(objphy)           # check sample names
sample_variables(objphy)       # sample variables from mdata
otu_table(objphy)              # ASV count table
tax_table(objphy)              # taxonomy for ASVs
taxa_names(objphy)             # ASV IDs

rank_names(objphy)             # check rank names
## if not already named like this, change rank names
colnames(tax_table(objphy)) <- c("Kingdom", "Phylum", "Class", "Order", "Family","Genus")


## calculate relative abundance
opr <- transform_sample_counts(objphy, function(x) x/sum(x))

## melt phyloseq object into dataframe
phylo_melt <- psmelt(opr)

### save object phylo_melt used to perform all following operations and remove not needed objects
save(phylo_melt, file = "arc_mod_GTDB_MAC_seq.RData")
rm(obj, objphy, opr, asv)
load("arc_mod_GTDB_MAC_seq.RData")



# rename taxa depending on max abundance ####
## use function 'sort_abundant_taxa'
source("C:/Users/admin/OneDrive/PhD/Nadja_BSc/16S_Amplicon/dada2_analysis_new_14102023/filter_taxa_above_threshold_corrected.R")



## taxa above 1% in at least one samples
# enter threshold in 'abundance_threshold' in percent
# this function creates a list as output
# for each rank (Phylum, Class, etc.) it contains a vector with all taxa which are in at least one sample
#above the set threshold
# it contains an abundance table similar to the input, with additional modified taxon names
#in which, for each rank, taxa which were below the threshold were renamed as e.g. 'other_o_Aminicenantales_<1%'
# but be careful, it does this via making it NA, so any before unclassified ASV on that rank, will be renamed
#by its supertaxon! but you can always recheck in the not modified taxon name

ta1 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 1, Abundance = "Abundance")
str(ta1)
# 14 genera, 13 families, 11 orders, 10 classes, 5 phyla

## taxa above 2% in at least one samples
ta2 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 2, Abundance = "Abundance")
str(ta2)
# 8 genera, 11 families, 8 orders, 8 classes, 5 phyla

# use 1% threshold
rm(ta2)


# order additional metadata columns for later plotting ####
ta1_original <- ta1

ta1$ASV_table_taxa_abv_1_perc$C_source <- factor(ta1$ASV_table_taxa_abv_1_perc$C_source, levels=c("Control",
                                                                                                  "Protein"),ordered=T)

ta1$ASV_table_taxa_abv_1_perc$Incubation_time <- factor(ta1$ASV_table_taxa_abv_1_perc$Incubation_time, levels = c("98", "157"), ordered = T)

save(ta1, file = "arc_mod_GTDB_long_rel_abd_table_bel1perc.RData")
load("arc_mod_GTDB_long_rel_abd_table_bel1perc.RData")

# plots ####
## plot on class level, not sorted ####

### which classes will be plotted
unique(ta1$ASV_table_taxa_abv_1_perc$Class_mod)

### subset to only column with class_mod as single taxonomy and Phylum_mod to group them by
ta1_class <- select(ta1$ASV_table_taxa_abv_1_perc, OTU, Sample, Class_mod, Phylum_mod, Abundance, Incubation_time, Sediment_depth, AB,  BES, C_source) %>% 
  mutate(Class_mod=factor(Class_mod, levels=c("Heimdallarchaeia","Lokiarchaeia","Thorarchaeia","other_p_Asgardarchaeota_<1%", "Methanosarcinia", "Syntropharchaeia", "other_p_Halobacteriota_<1%",
                                              "Nanoarchaeia","E2","EX4484-6","Thermoplasmata","other_p_Thermoplasmatota_<1%","Bathyarchaeia","other_p_Thermoproteota_<1%" , "other_Archaea_<1%"),
                          labels = c("Heimdallarchaeia","Lokiarchaeia","Thorarchaeia","other Asgardarchaeota", "Methanosarcinia", "Syntropharchaeia", "other Halobacteriota",
                                     "Nanoarchaeia","E2","EX4484-6","Thermoplasmata","other Thermoplasmatota","Bathyarchaeia","other Thermoproteota%" , "other Archaea"),
                          ordered = T), 
         Phylum_mod=factor(Phylum_mod)) %>% 
  group_by(Class_mod)



### sum ASVs of same taxonomy
ta1_class_s <- ta1_class %>% 
  group_by(Sample, Class_mod, Phylum_mod, Incubation_time, Sediment_depth, AB,  BES, C_source) %>% 
  summarise(Abundance = sum(Abundance), .groups = "keep") %>% 
  ungroup()


# Relabel the samples based on their day and replicate *(replicate 2)
ta1_class_s$Sample <- factor(ta1_class_s$Sample, levels = c("A3.1.day98_1","A3.2.day98_1","A3.1.day157_1","A3.2.day157_1",
                                                            "E3.1.day98","E3.2.day98","E3.1.day157","E3.2.day157"),
                             labels = c("98", "98*", "157", "157*",
                                        "98", "98*", "157", "157*"), ordered = T)



# color vector, manually selected
cols_all = c("#00441b","#a6dba0","#1b7837","#e0e0e0","#cab2d6","#40004b","#e0e0e0","#fde0dd","#abd9e9","#4393c3","#081d58","#e0e0e0","#ffffb3","#e0e0e0","#969696")


# plot only replicate one of Control and Protein amended samples for manuscript
ggplot(ta1_class_s[ta1_class_s$C_source %in% c("Control","Protein") & ta1_class_s$BES == "noBES" & ta1_class_s$AB == "AB" & ta1_class_s$Sediment_depth =="upper" & ta1_class_s$Sample %in% c("98","157"),], aes(x = Sample, y = Abundance, fill = Class_mod)) + 
  geom_bar(position = "fill", stat = "identity", col = "black", linewidth = 0.1) + 
  scale_fill_manual("Taxa", values = cols_all) + theme_bw() +
  theme(plot.title = element_text(size = 7, face = "bold",hjust=0),
        axis.text.x=element_text(size=7, hjust = 1, vjust = 0.5, angle = 90, color = "black"),
        axis.text.y=element_text(size=7, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face="bold", size=7, color = "black"),
        axis.title.x = element_text(face="bold", size=7, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=6),
        legend.title=element_text(size=7, face="bold"),
        legend.key.size = unit(0.35, "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size = 7, face="bold")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  ggtitle("Relative 16S rRNA gene abundance of archaea ") + 
  labs(x="Day", y = "relative abundance [%]")+
  guides(fill=guide_legend(ncol=1)) +
  facet_grid(~ C_source,scales = "free")

ggsave("arc_all_protein_AB_EX4484-6.png", width = 11, height = 10, units = "cm")



# Plots for Supplemental ####

## plot all Protein + AB amended samples for Archaea (Sl) ####
all_arc <- ggplot(ta1_class_s[ta1_class_s$C_source %in% c("Control","Protein") & ta1_class_s$BES == "noBES" & ta1_class_s$AB == "AB" & ta1_class_s$Sediment_depth =="upper",], aes(x = Sample, y = Abundance, fill = Class_mod)) + 
  geom_bar(position = "fill", stat = "identity", col = "black", linewidth = 0.1) + 
  scale_fill_manual("Taxa", values = cols_all) + theme_bw() +
  theme(plot.title = element_text(size = 7, face = "bold",hjust=0),
        axis.text.x=element_text(size=7, hjust = 1, vjust = 0.5, angle = 0, color = "black"),
        axis.text.y=element_text(size=7, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face="bold", size=7, color = "black"),
        axis.title.x = element_text(face="bold", size=7, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=6),
        legend.title=element_text(size=7, face="bold"),
        legend.key.size = unit(0.35, "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size = 7, face="bold")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  ggtitle("Relative 16S rRNA gene abundance of archaea ") + 
  labs(x="Day", y = "relative abundance [%]")+
  guides(fill=guide_legend(ncol=1)) +
  facet_grid(~ C_source,scales = "free")



## select only EX4484-6, color all ASV's with relative abundance >1 % and set all other taxa as other Archaea ####
## find higher abundance ASVs
list_taxa_above_threshold <- function(input_table, rank, threshold, abundance_column) {
  df1 <- input_table %>% 
    group_by(across({{rank}})) %>% 
    summarise(max_abundance = max({{abundance_column}})) %>% 
    filter(max_abundance > threshold/100 & is.na({{rank}}) == FALSE) %>% 
    arrange(desc(max_abundance))
  out_vector <- pull(df1, {{rank}})
  return(out_vector)
}

EX4484_6_test <- ta1_class %>% 
  filter(Class_mod == "EX4484-6")

abundance_threshold <- 1
ASV_ab_threshold_EX <- list_taxa_above_threshold(EX4484_6_test, OTU, abundance_threshold, Abundance)

EX4484_6_test_filt <- EX4484_6_test %>% 
  mutate(ASV_mod = if_else(OTU %in% ASV_ab_threshold_EX, OTU, paste0("other ASVs < ", abundance_threshold, "%"))) %>% 
  group_by(ASV_mod, Sample, Incubation_time, Sediment_depth, AB,  BES, C_source) %>% 
  summarise(Abundance = sum(Abundance), .groups = "keep") %>% 
  ungroup() %>% 
  left_join(unique(ta1_class[,c("OTU", "Class_mod",  "Phylum_mod")]), 
            by = c("ASV_mod" = "OTU"))

unique((
  EX4484_6_test_filt %>% 
    subset(ASV_mod != "other ASVs < 1%") %>% 
    mutate(ASV_taxa = paste(Class_mod, ASV_mod)) %>% 
    arrange(ASV_taxa)
)$ASV_taxa)


# Relabel the samples based on their day and replicate *(replicate 2)
EX4484_6_test_filt$Sample <- factor(EX4484_6_test_filt$Sample, levels = c("A3.1.day98_1","A3.2.day98_1","A3.1.day157_1","A3.2.day157_1",
                                                                          "E3.1.day98","E3.2.day98","E3.1.day157","E3.2.day157"),
                                    labels = c("98", "98*", "157", "157*",
                                               "98", "98*", "157", "157*"), ordered = T)


# Reorder Taxa in legend
EX4484_6_test_filt$ASV_mod <- factor(EX4484_6_test_filt$ASV_mod, levels = c("sq2","other ASVs < 1%"), ordered = T)


# set colors
cols_EX = c("#4393c3","#969696")

#plot
EX_ASV <- ggplot(EX4484_6_test_filt[EX4484_6_test_filt$C_source %in% c("Control","Protein") & EX4484_6_test_filt$BES == "noBES" & EX4484_6_test_filt$AB == "AB" & EX4484_6_test_filt$Sediment_depth == "upper",], aes(x = Sample, y = Abundance, fill = ASV_mod)) + 
  geom_bar(position = "stack", stat = "identity", col = "black", linewidth = 0.1) + 
  scale_fill_manual("Taxa", values = cols_EX) + theme_bw() +
  theme(plot.title = element_text(size = 7, face = "bold",hjust=0),
        axis.text.x=element_text(size=7, hjust = 1, vjust = 0.5, angle = 0, color = "black"),
        axis.text.y=element_text(size=7, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face="bold", size=7, color = "black"),
        axis.title.x = element_text(face="bold", size=7, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=6),
        legend.title=element_text(size=7, face="bold"),
        legend.key.size = unit(0.35, "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size = 7, face="bold")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), labels = scales::percent) +
  ggtitle("EX4484-6 ASVs ") + 
  labs(x="Day", y = "relative abundance [%]")+
  guides(fill=guide_legend(ncol=1)) +
  facet_grid(~ C_source,scales = "free")



## select only Lokiarchaeia, color all ASV's with relative abundance >1 % and set all other taxa as other Archaea ####
## find higher abundance ASVs

Loki_test <- ta1_class %>% 
  filter(Class_mod == "Lokiarchaeia") 

abundance_threshold <- 1
ASV_ab_threshold_Loki <- list_taxa_above_threshold(Loki_test, OTU, abundance_threshold, Abundance)

Loki_test_filt <- Loki_test %>% 
  mutate(ASV_mod = if_else(OTU %in% ASV_ab_threshold_Loki, OTU, paste0("other ASVs < ", abundance_threshold, "%"))) %>% 
  group_by(ASV_mod, Sample, Incubation_time, Sediment_depth, AB,  BES, C_source) %>% 
  summarise(Abundance = sum(Abundance), .groups = "keep") %>% 
  ungroup() %>% 
  left_join(unique(ta1_class[,c("OTU", "Class_mod",  "Phylum_mod")]), 
            by = c("ASV_mod" = "OTU"))

Loki_test_filt <- Loki_test_filt %>% 
  mutate(Class_mod = if_else(is.na(Class_mod), "others", Class_mod))

unique((
  Loki_test_filt %>% 
    subset(ASV_mod != "other ASVs < 1%") %>% 
    mutate(ASV_taxa = paste(Class_mod, ASV_mod)) %>% 
    arrange(ASV_taxa)
)$ASV_taxa)


# Reorder and relabel the samples based on their day and replicate *(replicate 2)
Loki_test_filt$Sample <- factor(Loki_test_filt$Sample, levels = c("A3.1.day98_1","A3.2.day98_1","A3.1.day157_1","A3.2.day157_1",
                                                                  "E3.1.day98","E3.2.day98","E3.1.day157","E3.2.day157"),
                                labels = c("98", "98*", "157", "157*",
                                           "98", "98*", "157", "157*"), ordered = T)


# Reorder Taxa in legend
Loki_test_filt$ASV_mod <- factor(Loki_test_filt$ASV_mod, levels = c("sq1","sq3","sq6","sq10","sq11","sq12","sq15",
                                                                    "sq18","sq29","sq32","sq41","sq48","sq55","sq57","sq131","other ASVs < 1%"), ordered = T)


# set colors
cols_Loki = c("#a6dba0","#1b7837","#a6cee3","#1f78b4","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#8dd3c7","#ffffb3","#80b1d3","#fccde5","#bc80bd","#969696")

#plot
Loki_ASV <- ggplot(Loki_test_filt[Loki_test_filt$C_source %in% c("Control","Protein") & Loki_test_filt$BES == "noBES" & Loki_test_filt$AB == "AB" & Loki_test_filt$Sediment_depth == "upper",], aes(x = Sample, y = Abundance, fill = ASV_mod)) + 
  geom_bar(position = "stack", stat = "identity", col = "black", linewidth = 0.1) + 
  scale_fill_manual("Taxa", values = cols_Loki) + theme_bw() +
  theme(plot.title = element_text(size = 7, face = "bold",hjust=0),
        axis.text.x=element_text(size=7, hjust = 1, vjust = 0.5, angle = 0, color = "black"),
        axis.text.y=element_text(size=7, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face="bold", size=7, color = "black"),
        axis.title.x = element_text(face="bold", size=7, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=6),
        legend.title=element_text(size=7, face="bold"),
        legend.key.size = unit(0.35, "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size = 7, face="bold")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1), labels = scales::percent) +
  ggtitle("Lokiarchaeia ASVs ") + 
  labs(x="Day", y = "relative abundance [%]")+
  guides(fill=guide_legend(ncol=1)) +
  facet_grid(~ C_source,scales = "free")

#Combine Plots for all samples
library(patchwork)

layout <- "
AA
BB
CC"

all_arc + EX_ASV + Loki_ASV + plot_layout(design = layout) + plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 7),
        legend.justification = "left")

ggsave("arc_all_protein_AB_EX4484-6_Sl.png", width = 16, height = 26, units = "cm")
ggsave("arc_all_protein_AB_EX4484-6_Sl.pdf", width = 16, height = 26, units = "cm")
