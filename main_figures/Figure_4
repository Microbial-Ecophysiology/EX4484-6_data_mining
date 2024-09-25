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


# Figure 4 ####
## scatter plot of relative abundance for each redundant sample ####
rel_abundance_MAG <- ggplot() +
  geom_jitter(data=tax_EX_95_melt[tax_EX_95_melt$environment_detailed != "marine water",],  width =0.25,
              aes(x=MAG, y=rel_abundance_nas, colour = environment_detailed), size = 1)+
  geom_jitter(data=tax_EX_95_melt,  width =0.25,
              aes(x=MAG, y=rel_abundance_nas, colour = environment_detailed, shape = filter_category), size = 0.9)+
  ggtitle("")+
  labs(x="", y="relative abundance [%]", color = "environment", shape = "marine water size fraction")+
  theme_bw()+
  theme(plot.title = element_text(size = 7, face = "bold"),
        legend.title = element_text(size=7, face="bold", color="black"),
        legend.text = element_text(size=7, color="black"),
        plot.background = element_rect(fill = 'white', color = NA),
        axis.title = element_text(size = 7, face= "bold", color="black"),
        axis.text.x = element_text(size = 7, angle = 45, hjust =1, color="black"),
        axis.text = element_text(size = 7, color="black"),
        legend.key=element_rect(fill="white"),
        legend.position="right",
        legend.box="vertical", 
        legend.margin=margin(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "#d9d9d9", linewidth = 0.25, linetype = 3))+
  scale_y_continuous(breaks = scales::pretty_breaks(n=7), expand = c(0,0), limits = c(0,0.36))+
  scale_color_manual(values = c("cold seep" = "#000000",  "hydrothermal vent" = "#e31a1c","lake sediment"= "#d9d9d9", "marine sediment" = "#74add1", "marine water" = "#07316A"),
                     labels = c("cold seep fluid", "marine hydrothermal sediment", "lake sediment", "marine sediment","marine water"),
                     drop = F)+
  scale_shape_manual(labels = c("free-living (0.2-5 µm)", "particle attached (>5 µm)", "bulk (>0.2 µm)","no meta data available"),values = c(1, 17, 15, 4), na.translate = F)+
  guides(shape = guide_legend(order = 2, override.aes = list(size=2)),col = guide_legend(order = 1, override.aes = list(size=2)))

ggsave(plot = rel_abundance_MAG, "environment_scatter_plot.png", width = 18, height = 10, units = "cm")


## average abundance per MAG

n_observations_tax_EX_95_melt <- tax_EX_95_melt %>% 
  group_by(Family) %>%
  summarise(count = n())

# 1A: 1280
# 1B: 256
# 2: 256
# 3A: 384
# 3B: 128
# 4: 128

n_observations_tax_EX_95_total <- tax_EX_95_melt %>% 
  summarise(count = n())

# total: 2432

plot_1 <- ggplot(data=subset(tax_EX_95_melt, !is.na(rel_abundance_nas)), aes(x= Family, y=rel_abundance_nas))+
  geom_boxplot(linewidth= 0.5, fatten = 1, outlier.size = 0.7)+
  ggtitle("")+
  labs(x="Family", y="relative abundance [%]", color = "environment")+
  theme_bw()+
  theme(plot.title = element_text(size = 7, face = "bold"),
        legend.title = element_text(size=7, face="bold", color="black"),
        legend.text = element_text(size=7, color="black"),
        axis.title = element_text(size = 7, face= "bold", color="black"),
        axis.text.x = element_text(size = 7, color="black"),
        axis.text = element_text(size = 7, color="black"),
        legend.key=element_rect(fill="white"),
        legend.position="right",
        legend.box="vertical", 
        legend.margin=margin(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))+
  scale_y_continuous(breaks = scales::pretty_breaks(n=7), expand = c(0,0))


plot_2 <- ggplot(data=subset(tax_EX_95_melt, !is.na(rel_abundance_nas)), aes(x= "total" , y=rel_abundance_nas))+
  geom_boxplot(linewidth = 0.5, fatten = 1, width=0.9, outlier.size = 0.7)+
  ggtitle("")+
  labs(x="total", y="", color = "environment")+
  theme_bw()+
  theme(plot.title = element_text(size = 7, face = "bold"),
        legend.title = element_text(size=7, face="bold", color="black"),
        legend.text = element_text(size=7, color="black"),
        axis.title = element_text(size = 7, face = "bold", color="black"),
        axis.text = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key=element_rect(fill="white"),
        legend.position="right",
        legend.box="vertical", 
        legend.margin=margin(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA))+
  scale_y_continuous(breaks = scales::pretty_breaks(n=7), expand = c(0,0))

box_rel_abundance <- wrap_elements( plot_1 + plot_2 + plot_layout(widths = c(6,1))+ 
                                      plot_annotation(title = "b") &
                                      theme(plot.title = element_text(face="bold", size = 7),
                                            panel.background = element_rect(fill = "transparent"),
                                            panel.border = element_rect(fill = "NA"),
                                            plot.background = element_rect(fill = "transparent", color = NA)))

rel_abundance_MAG_plot <- wrap_elements( rel_abundance_MAG +  
                                           plot_annotation(title = "a") &
                                           theme(plot.title = element_text(face="bold", size = 7)))

combined_plot <- rel_abundance_MAG_plot + inset_element(box_rel_abundance, left = 0.09, bottom = 0.55, right = 0.55, top = 0.92)

ggsave(plot =combined_plot ,"environment_box_plot_mod.png", width = 18, height = 17, units = "cm")
ggsave(plot =combined_plot ,"environment_box_plot_mod.pdf", width = 18, height = 17, units = "cm")
ggsave(plot =combined_plot ,"environment_box_plot_mod.emf", width = 18, height = 17, units = "cm")

save.image("Thermo_rel_abundance.RData")
load("Thermo_rel_abundance.RData")
