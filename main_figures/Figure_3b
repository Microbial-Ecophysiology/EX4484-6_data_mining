# Script to combine results of post analysis coverm over 8573 metagenomic runs

library(reshape2)
library(tidyverse)
library(ggplot2)
library(maps)
library(patchwork)


# import data ####
setwd("C:/Dokumente und Einstellungen/admin/OneDrive/PhD/Data_mining/CoverM/v3/")

load("Thermo_rel_abundance.RData")

# Using 95 % similarity and breadth of coverage with 50%

breadth_95 <- read.table("Thermo_breadth_95.txt", header = T)
rel_95 <- read.table("Thermo_rel_95.txt", header = T)


# Only samples with a breadth of coverage of at least 50% are chosen for further analyses ####
# Create a logical condition to check if any value in a row is greater than 0.5
condition_95 <- apply(breadth_95, 1, function(row) any(row > 0.5))

# for 50 % breadth of coverage
# Subset the data frame based on the condition
subsetted_data_95 <- breadth_95[condition_95, ]

# Extract the row names from the subsetted_data
breadth_50_accessions_95 <- rownames(subsetted_data_95)

# Subset rel_95 based on the interesting row names
rel_95_samples <- rel_95[breadth_50_accessions_95, , drop = FALSE]

# further, rel_95_samples is used to combine this data with initial sample list and create a heatmap based on the sample environments
meta <- read.csv("results_MG_mining.txt", 
                 sep="\t", 
                 header=TRUE, 
                 fill=TRUE, 
                 quote="")

# run accession in column 79
rownames(meta) <- meta[,79]

# Subset meta based on the interesting row names
meta_subset_samples_95 <- meta[breadth_50_accessions_95, , drop = FALSE]

# remove large meta file
rm(meta)

# subset columns to 10: base count, 62: lat, 70: lon, 92: scientific name
meta_subset_samples_95_short <- select(meta_subset_samples_95,"lat", "lon")

sample_merged_95 <- merge(rel_95_samples, meta_subset_samples_95_short, 
                          by = 'row.names', all = TRUE) 

# remove unmapped column
sample_merged_95 <- sample_merged_95[, -2]

save.image("Thermo_rel_abundance.RData")


# add information about Taxonomy ####

sample_merged_95_melt <- melt(sample_merged_95, id.vars = c("Row.names","lat","lon"))

Taxonomy <- read.csv("GTDB_classification.txt", 
                     sep="\t", 
                     header=TRUE, 
                     fill=TRUE, 
                     quote="")

tax_merged_95 <- merge(sample_merged_95_melt, Taxonomy, 
                       by.x = "variable", by.y = "user_genome", all = TRUE) 


# remove rows with NA in Row.names, as these are the rows with redundant EX MAGs
tax_merged_95 <- tax_merged_95[!is.na(tax_merged_95$Row.names),]

unique(tax_merged_95$class) # 11 classes

# "EX4484-6", "Thermoplasmata", "Poseidoniia","E2", "UBA186","SW-10-69-26"
# outgroups: "Korarchaeia", "Methanomethylicia", "Bathyarchaeia","Methanosarcinia", "Nitrososphaeria"

# subset all Thermoplasmatota
tax_Thermo_95 <- tax_merged_95[tax_merged_95$phylum == 'Thermoplasmatota',]
write.table(tax_Thermo_95, "tax_Thermo_all_95.txt", row.names = T, col.names = T, sep = "\t", quote = F)

tax_Thermo_95 <- tax_Thermo_95[tax_Thermo_95$value != 0, ]
write.table(tax_Thermo_95, "tax_Thermo_95.txt", row.names = T, col.names = T, sep = "\t", quote = F)

# From this data, now subset all EX4484-6 MAGs
tax_EX_95 <- tax_merged_95[tax_merged_95$class == 'EX4484-6',]
tax_EX_95 <- tax_EX_95[tax_EX_95$value != 0, ]
write.table(tax_EX_95, "tax_EX_95.txt", row.names = F, col.names = T, sep = "\t", quote = F)


# subset all unique run accessions and create a file with detailed environment for all EX4484-6 MAGs
unique_EX_run_acc <- unique(tax_EX_95$Row.names) # 128 unique run_accessions
subset_EX_samples <- meta_subset_samples_95[meta_subset_samples_95$run_accession %in% unique_EX_run_acc,]
write.table(subset_EX_samples, "subset_EX_samples.txt", row.names = T, col.names = T, sep = "\t", quote = F)

mat_EX <- reshape::cast(tax_EX_95, variable~Row.names, value = "value", fun.aggregate = "mean", fill = 0)

tax_EX_95_melt <- melt(mat_EX)
colnames(tax_EX_95_melt) <- c("MAG","run_accession","rel_abundance")

# reload manually refined data
EX_sample_info <- read.csv("subset_EX_samples_mod.txt", 
                           sep="\t", 
                           header=TRUE, 
                           fill=TRUE, 
                           quote="")


# subset columns 
EX_sample_info_short <- select(EX_sample_info, "X", "lat", "lon", "environment_detailed", "info_sampling", "environment_combined", "water_depth", "filter_category")

# add family 
family <- read.csv("new_family_groups.txt", 
                   sep="\t", 
                   header=TRUE, 
                   fill=TRUE, 
                   quote="")

# add info about environment and family to tax_EX_95

tax_EX_95_melt <- full_join(x = tax_EX_95_melt, y = family, by = c("MAG"="MAG")) %>%
  full_join(., EX_sample_info_short, by = c("run_accession" = "X"), suffix = c("_EX","_run_acc"))

tax_EX_95_melt <- tax_EX_95_melt %>% drop_na(run_accession)

tax_EX_95_melt$MAG <- factor(tax_EX_95_melt$MAG, levels = c("E3_1_d157_spades_bin_16_orig_1_refined.contigs","PRJNA889212_bin_7_strict_refined.contigs","PRJNA541421_bin_70_strict.contigs","MSM105_N25025F_bin_277_ori_permissive_1","PRJNA368391_bin_183_orig.contigs",
                                                            "EMB267_Co_bin_434_strict_1","PRJNA721298_bin_46_orig_refined.contigs","PRJNA889212_bin_1_orig_refined.contigs","GCA_016928095.1_ASM1692809v1_genomic",
                                                            "PRJNA531756_bin_93_strict_refined.contigs","GCA_023145185.1_ASM2314518v1_genomic","ZORZ22.1_SAMN30647027_MAG_00000060","ZHEN22.1_SAMN22703514_MAG_00000188",
                                                            "GCA_003650025.1_ASM365002v1_genomic","GCA_021160765.1_ASM2116076v1_genomic","PRJNA531756_bin_72_orig_refined.contigs","PRJNA541421_bin_98_orig_refined.contigs",
                                                            "TARA_SAMEA2620113_MAG_00000097","PRJNA704804_bin_147_orig_refined.contigs","SUTE22.1_SAMN10231904_MAG_00000012","SUTE22.1_SAMN10231914_MAG_00000292","SUTE22.1_SAMN10231914_MAG_00000186"),
                             labels = c("E3_1_d157","PRJNA889212_2","PRJNA541421_1","MSM105","PRJNA368391_1","EMB267","PRJNA721298_1","PRJNA889212_1","GCA_016928095.1","PRJNA531756_2","GCA_023145185.1",
                                        "SAMN30647027_1","SAMN22703514_1","GCA_003650025.1","GCA_021160765.1","PRJNA531756_1","PRJNA541421_3","SAMEA2620113_1","PRJNA704804_1","SAMN10231904_1","SAMN10231914_2","SAMN10231914_1"),ordered=T)

tax_EX_95_melt$filter_category <- factor(tax_EX_95_melt$filter_category, levels = c("free-living","particle attached", "bulk", "no meta data available"), ordered = T)

tax_EX_95_melt$rel_abundance_nas <- ifelse(tax_EX_95_melt$rel_abundance==0,NA,tax_EX_95_melt$rel_abundance)


# Figure 3b ####
## combined plot for all families  ####

tax_EX_95_map  <- tax_EX_95_melt[tax_EX_95_melt$rel_abundance != 0, ]

vals <- c("cold seep" = "#000000",  "hydrothermal vent" = "#e31a1c","lake sediment"= "#d9d9d9", "marine sediment" = "#74add1", "marine water" = "#07316A")



facet_labels <- c(
  `1A` = "Family 1A",
  `1B` = "Family 1B",
  `2` = "Family 2",
  `3A` = "Family 3A",
  `3B` = "Family 3B",
  `4` = "Family 4"
)


## relative abundance [%] ####
World <- map_data("world")

ggplot() +
  geom_polygon(data = World, aes(x=long, y =lat, group = group), fill="grey")+
  geom_point(data=tax_EX_95_map, 
             aes(x=lon, y=lat,  fill = environment_detailed ,size = rel_abundance),shape = 21, colour = "#252525")+ 
  ggtitle("")+
  labs(x="longitude", y="latitude",size="relative abundance [%]", color = "environment", fill = "environment")+
  theme(plot.title = element_text(size = 7, face = "bold"),
        legend.title = element_text(size=7, face="bold"),
        legend.text = element_text(size=7),
        axis.title = element_text(size = 7, face = "bold"),
        axis.text = element_text(size = 7),
        legend.key=element_rect(fill="white"),
        legend.position="right",
        legend.box="vertical", 
        legend.margin=margin(),
        panel.background = element_rect(fill = "#f0f0f0"),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text.x = element_text(size = 7, face="bold"))+
  scale_color_manual(values = vals, limits = names(vals))+
  scale_fill_manual(values = vals, limits = names(vals))+
  scale_size_continuous(limits=c(0.00000001,0.35), breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3), range=c(0,4))+
  guides(color = guide_legend(order = 2), size =guide_legend(order=1),fill = guide_legend(override.aes = list(size=3)))+
  facet_wrap(~Family, ncol = 3,labeller = as_labeller(facet_labels))

ggsave("Overview_MG_searched_relabd_95similarity.png", width = 19, height = 8, device="png", units = "cm")
ggsave("Overview_MG_searched_relabd_95similarity.pdf", width = 19, height = 8, device="pdf", units = "cm")
ggsave("Overview_MG_searched_relabd_95similarity.emf", width = 19, height = 8, device="emf", units = "cm")
