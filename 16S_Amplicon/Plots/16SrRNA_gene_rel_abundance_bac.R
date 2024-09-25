## script for subsetting a taxmap object to a final long format table with taxa above a certain threshold

# input:
# taxmap object (created in first steps here)
# undesired taxa sorted out
# asv counts, no relative abundance calculated yet
# mdata associated with asv table

# output:
# 1. list of taxa above certain threshold (can be set) in at least one sample for each rank (genus, family, etc.) 
# 2. long format table with samples and individual ASVs with additional columns for each rank renamed taxa name if ASV was below threshold
# 3. long format table with samples and taxa summed up per taxa from table 2


# packages ####
library(tidyverse)
library(taxa)
library(metacoder)
library(phyloseq)
library(hues)
library(cowplot)
library(devEMF)


# import data ####
setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Nadja_BSc/16S_Amplicon/dada2_analysis_new_14102023/Plot_rel_abundance_Bacteria/")

## ASV table raw counts
asv <- read_tsv("../Bacteria_all/results/asv_table_tax_Lib28_bac_GTDB_mod.txt")

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
# 2206 taxa (including subtaxa)
# 5256 ASVs (=number of rows of tax_data table)
ini_reads <- sum(obj$data$asv_table[, sample.names])
ini_reads
# 460332 reads in original table

## modify taxmap object ####
### remove any ASVs which are not Bacteria
obj <- metacoder::filter_taxa(obj,
                              taxon_names == "Bacteria",
                              subtaxa = T)
sum(obj$data$asv_table[, sample.names])/ini_reads
#  1

### remove doubletons and singletons
obj$data$asv_table <- metacoder::zero_low_counts(obj, "asv_table",
                                                 min_count = 3,
                                                 use_total = T,
                                                 other_cols = T)
#Zeroing 781 of 5256 rows with total counts less than 3

### check for empty ASVs
no_reads <- rowSums(obj$data$asv_table[, sample.names]) == 0
sum(no_reads)  # 3007 empty ASVs

# remove empty ASVs
obj <- metacoder::filter_obs(obj, "asv_table",
                             ! no_reads,
                             drop_taxa = T)
obj
# taxa reduced to 1195
# ASVs reduced to 2249
sum(obj$data$asv_table[, sample.names])/ini_reads
# 99.75% of the reads kept

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
save(obj, file = "taxmap_MAC_bac_GTDB_mod.RData")

load("taxmap_MAC_bac_GTDB_mod.RData")
# here named 'obj'


# create phyloseq object ####
# maybe need to change names of otu_table, otu_id_col and sample names (= column name of column with sample names)
objphy <- as_phyloseq(obj, otu_table = "asv_table", otu_id_col = "ASV_ID", sample_data = mdata, sample_id_col = "Sample")

ntaxa(objphy)                  # how many taxa (including subtaxa) - 2249
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
save(phylo_melt, file = "bac_GTDB_MAC_seq_mod.RData")
rm(obj, objphy, opr, asv)
load("bac_GTDB_MAC_seq_mod.RData")



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
# 10 genera, 10 families, 13 orders, 13 classes, 10 phyla

## taxa above 2% in at least one samples
ta2 <- sort_abundant_taxa(table = phylo_melt, abundance_threshold = 2, Abundance = "Abundance")
str(ta2)
# 7 genera, 8 families, 9 orders, 11 classes, 8 phyla

# use 2% threshold
rm(ta1)



# order additional metadata columns for later plotting ####
ta2_original <- ta2

ta2$ASV_table_taxa_abv_2_perc$C_source <- factor(ta2$ASV_table_taxa_abv_2_perc$C_source, levels=c("Control",
                                                                                                  "Protein"),ordered=T)

ta2$ASV_table_taxa_abv_2_perc$Incubation_time <- factor(ta2$ASV_table_taxa_abv_2_perc$Incubation_time, levels = c("98", "157"), ordered = T)


save(ta2, file = "bac_GTDB_long_rel_abd_table_bel2perc_mod.RData")
load("bac_GTDB_long_rel_abd_table_bel2perc_mod.RData")

# plots ####
## plot on class level, not sorted ####

### which classes will be plotted
unique(ta2$ASV_table_taxa_abv_2_perc$Class_mod)

### subset to only column with class_mod as single taxonomy and Phylum_mod to group them by
ta2_class <- select(ta2$ASV_table_taxa_abv_2_perc, OTU, Sample, Class_mod, Phylum_mod, Abundance, Incubation_time, Sediment_depth, AB,  BES, C_source) %>% 
  mutate(Class_mod=factor(Class_mod, levels=c("Aminicenantia", "other_p_Acidobacteriota_<2%","JS1","Bacilli","Campylobacteria","Anaerolineae","Dehalococcoidia","other_p_Chloroflexota_<2%",
                                              "Desulfobacteria","Desulfovibrionia","other_p_Desulfobacterota_<2%", "Phycisphaerae","Planctomycetia","other_p_Planctomycetota_<2%",
                                              "Gammaproteobacteria","other_p_Pseudomonadota_<2%", "other_Bacteria_<2%"), 
                          ordered = T), 
         Phylum_mod=factor(Phylum_mod)) %>% 
  group_by(Class_mod)



### sum ASVs of same taxonomy
ta2_class_s <- ta2_class %>% 
  group_by(Sample, Class_mod, Phylum_mod, Incubation_time, Sediment_depth, AB,  BES, C_source) %>% 
  summarise(Abundance = sum(Abundance), .groups = "keep") %>% 
  ungroup()


# Relabel the samples based on their day and replicate *(replicate 2)
ta2_class_s$Sample <- factor(ta2_class_s$Sample, levels = c("A3.1.day98_1","A3.2.day98_1","A3.1.day157_1","A3.2.day157_1",  
                                                            "E3.2.day98","E3.1.day98", "E3.1.day157","E3.2.day157"),
                             labels = c("98", "98*", "157", "157*",
                                        "98", "98*", "157", "157*"), ordered = T)



# color vector, manually selected

cols_all = c("#8e0152","#e0e0e0","#f1b6da","#ffffb3","#dfc27d","#abd9e9","#4393c3","#e0e0e0","#238b45","#a6dba0","#e0e0e0","#40004b","#cab2d6","#e0e0e0","#053061","#e0e0e0","#969696")

# plot only replicate one of Control and Protein amended samples for manuscript
ggplot(ta2_class_s[ta2_class_s$C_source %in% c("Control","Protein") & ta2_class_s$BES == "noBES" & ta2_class_s$AB == "AB" & ta2_class_s$Sediment_depth =="upper",], aes(x = Sample, y = Abundance, fill = Class_mod)) + 
  geom_bar(position = "fill", stat = "identity", col = "black", linewidth = 0.1) + 
  scale_fill_manual("Taxa", values = cols_all) + theme_bw() +
  theme(plot.title = element_text(size = 7, face = "bold",hjust=0),
        axis.text.x=element_text(size=7, hjust = 1, vjust = 0.5, angle = 0),
        axis.text.y=element_text(size=7),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face="bold", size=7),
        axis.title.x = element_text(face="bold", size=7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=6),
        legend.title=element_text(size=7, face="bold"),
        legend.key.size = unit(0.35, "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size = 7, face="bold")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  ggtitle("Relative 16S rRNA gene abundance of bacteria ") + 
  labs(x="Day", y = "Relative abundance")+
  guides(fill=guide_legend(ncol=1)) +
  facet_grid(~ C_source,scales = "free")

ggsave("bac_all_protein_AB_Sl.png", width = 15, height = 10, units = "cm")


