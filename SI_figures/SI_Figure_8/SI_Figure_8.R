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

# only check relative abundance for detected Thermoplasmatota
tax_Thermo_95 <- tax_Thermo_95[tax_Thermo_95$value != 0, ]
write.table(tax_Thermo_95, "tax_Thermo_95.txt", row.names = T, col.names = T, sep = "\t", quote = F)


#add new order based on more recent publications
new_order <- read.csv("GTDB_new_order.txt", 
                     sep="\t", 
                     header=TRUE, 
                     fill=TRUE, 
                     quote="")

tax_merged_95_new_order <- merge(tax_Thermo_95, new_order, 
                       by.x = "variable", by.y = "user_genome", all = TRUE) 
tax_merged_95_new_order <- tax_merged_95_new_order[!is.na(tax_merged_95_new_order$Row.names),]


n_observations_plot1 <- tax_merged_95_new_order %>% count(new.order)

tax_merged_95_new_order$class <- factor(tax_merged_95_new_order$class, levels = c("E2", "EX4484-6","Poseidoniia","SW-10-69-26","Thermoplasmata","UBA186"), ordered = T)
tax_merged_95_new_order$new.order <- factor(tax_merged_95_new_order$new.order, levels = c("DHVEG-1","UBA202","UBA9212",
                                                                              "EX4484-6_1","EX4484-6_2","EX4484-6_3","EX4484-6_4",
                                                                              "Poseidoniales","MGIII",
                                                                              "JACQPN01","SW-10-69-26",
                                                                              "Aciduliprofundales","Ca. Angelarchaeales","Ca. Gimiplasmatales","Ca. Lutacidiplasmatales","Ca. Sysuiplasmatales","Methanomassiliicoccales",
                                                                              "Thermoplasmatales","ARK-15","B87-G9","PWKY01","SG8-5","Thermoplasmata_unclassified",
                                                                              "DTKX01",
                                                                              "UBA186","UBA287"),
                                            labels = c("DHVEG-1 (N=181)","UBA202 (N=81)","UBA9212 (N=4)",
                                                       "EX4484-6_1 (N=95)","EX4484-6_2 (N=18)","EX4484-6_3 (N=90)","EX4484-6_4 (N=7)",
                                                       "Poseidoniales (N=12537)","MGIII (N=3232)",
                                                       "JACQPN01 (N=0)","SW-10-69-26 (N=0)",
                                                       "Aciduliprofundales (N=18)","Ca. Angelarchaeales (N=60)","Ca. Gimiplasmatales (N=49)","Ca. Lutacidiplasmatales (N=902)","Ca. Sysuiplasmatales (N=113)","Methanomassiliicoccales (N=224)",
                                                       "Thermoplasmatales (N=2291)","ARK-15 (N=9)","B87-G9 (N=9)","PWKY01 (N=90)","SG8-5 (N=17)","Thermoplasmata_unclassified (N=2)",
                                                       "DTKX01 (N=2)",
                                                       "UBA186 (N=1)","UBA287 (N=0)"),ordered = T)

## plot relative abundance [%] all Thermo####
cols <- c("#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf","#2166ac")

plot1 <- ggplot(data=tax_merged_95_new_order,
                             aes(x=new.order, y=value, fill = class)) +
  geom_boxplot(outlier.colour="black", outlier.shape=1,
               outlier.size=0.9, linewidth = 0.2) +
  labs(x="", y="relative abundance [%] (sqrt)",color="class")+
  theme_bw()+
  theme(axis.title.y = element_text(face="bold", size=7),
        axis.title.x = element_text(face="bold", size=7),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7),
        axis.text.y = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7, face = "bold"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill='white', color=NA),
        legend.background = element_rect(fill='white'))+
  scale_y_sqrt(expand = c(0,0),limits= c(0,50),breaks = c(0,1,5,10,20,30,40))+
  scale_fill_manual(values = cols, drop = F)

ggsave(plot=rel_abd_all_Thermo,"rel_abd_all_Thermo.png", width = 18, height = 15, units = "cm")
ggsave(plot=rel_abd_all_Thermo,"rel_abd_all_Thermo.pdf", width = 18, height = 15, units = "cm")
ggsave(plot=rel_abd_all_Thermo,"rel_abd_all_Thermo.emf", width = 18, height = 15, units = "cm")

save.image("Thermo_rel_abundance.RData")




# Combine percentage of unknown genes with relative abundance data

setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Data_mining/Annotation_v3/hypothetical_allThermo")

# load data ####
rel_abundance <- read.csv("../../CoverM/v3/rel_abundance_tax_Thermo_all_95_mod_v2.txt", 
                          sep="\t", 
                          header=TRUE, 
                          fill=TRUE, 
                          quote="",
                          col.names = c("genome","accession","lat","lon","rel_abundance","domain","phylum","class","order"))

EX_data <- read.csv("../hypothetical/unknown_percent_EX4484-6.txt",
                    sep="\t", 
                    header=TRUE, 
                    fill=TRUE, 
                    quote="",
                    col.names = c("genome","annotated","hypothetical","unknown","percentage_annotated","percentage_hypothetical","percentage_unknown"))

Taxonomy_EX <- read.csv("GTDB_classification_mod.txt",
                        sep="\t", 
                        header=TRUE, 
                        fill=TRUE, 
                        quote="")

Thermo_data <- read.csv("unknown_percent_Thermoplasmatota.txt",
                        sep="\t", 
                        header=TRUE, 
                        fill=TRUE, 
                        quote="")

new_order <- read.csv("GTDB_new_order.txt",
                      sep="\t", 
                      header=TRUE, 
                      fill=TRUE, 
                      quote="")


## add information about taxonomy for EX MAGs ####
EX_data_merged <- merge(EX_data, Taxonomy_EX, 
                        by.x = "genome", by.y = "user_genome", all = TRUE) 


# remove all non EX4484-6 MAGs
EX_data_merged <- EX_data_merged[EX_data_merged$class == "EX4484-6",]
EX_data_merged <- subset(EX_data_merged, select = -c(percentage_annotated, percentage_hypothetical) )

# combine Thermoplasmatota with EX table
all_Thermo <- rbind(EX_data_merged, Thermo_data)


# combine rel abundance dataframe with Thermo percentage hypothetical
all_Thermo_subset <- all_Thermo %>%
  select(genome, percentage_unknown)

df_all <- full_join(x = rel_abundance, y = all_Thermo_subset, by = c("genome"))

# remove redundant EX_MAGs and remove all rows, in which rel abundance = 0
df_all <- df_all[!is.na(df_all$domain),]

# plot fraction of abundance in environment

genome_percentage <- df_all %>%
  group_by(genome) %>%
  summarise(fraction = sum(rel_abundance > 0) / length(unique(accession)) * 100,
            percentage_unknown = mean(percentage_unknown, na.rm = TRUE),
            class = first(class),
            order = first(order)
  )

genome_percentage_order <- merge(genome_percentage, new_order,
                                 by.x = "genome", by.y = "user_genome", all = TRUE)

genome_percentage_all <- genome_percentage_order[!is.na(genome_percentage_order$order),]


# manually change order for EX4484-6 MAGs to newly defined orders

write.table(genome_percentage_order, "genome_percentage.txt", 
            quote = FALSE,
            sep="\t",
            col.names = TRUE,
            row.names = FALSE)

genome_percentage_order <- read.csv("genome_percentage.txt",
                                    sep="\t", 
                                    header=TRUE, 
                                    fill=TRUE, 
                                    quote="")




genome_percentage_all$class <- factor(genome_percentage_all$class, levels = c("E2", "EX4484-6","Poseidoniia","SW-10-69-26","Thermoplasmata","UBA186"), ordered = T)

n_observations_plot2 <- genome_percentage_all %>% count(new.order)

genome_percentage_all$new.order <- factor(genome_percentage_all$new.order, levels = c("DHVEG-1","UBA202","UBA9212",
                                                                                      "EX4484-6_1","EX4484-6_2","EX4484-6_3","EX4484-6_4",
                                                                                      "Poseidoniales","MGIII",
                                                                                      "JACQPN01","SW-10-69-26",
                                                                                      "Aciduliprofundales","Ca. Angelarchaeales","Ca. Gimiplasmatales","Ca. Lutacidiplasmatales","Ca. Sysuiplasmatales","Methanomassiliicoccales",
                                                                                      "Thermoplasmatales","ARK-15","B87-G9","PWKY01","SG8-5","Thermoplasmata_unclassified",
                                                                                      "DTKX01",
                                                                                      "UBA186","UBA287"),
                                          labels = c("DHVEG-1 (N=21)","UBA202 (N=12)","UBA9212 (N=3)",
                                                     "EX4484-6_1 (N=15)","EX4484-6_2 (N=2)","EX4484-6_3 (N=4)","EX4484-6_4 (N=1)",
                                                     "Poseidoniales (N=67)","MGIII (N=10)",
                                                     "JACQPN01 (N=1)","SW-10-69-26 (N=1)",
                                                     "Aciduliprofundales (N=10)","Ca. Angelarchaeales (N=15)","Ca. Gimiplasmatales (N=10)","Ca. Lutacidiplasmatales (N=43)","Ca. Sysuiplasmatales (N=2)","Methanomassiliicoccales (N=96)",
                                                     "Thermoplasmatales (N=64)","ARK-15 (N=3)","B87-G9 (N=1)","PWKY01 (N=12)","SG8-5 (N=7)","Thermoplasmata_unclassified (N=1)",
                                                     "DTKX01 (N=2)",
                                                     "UBA186 (N=1)","UBA287 (N=1)"),ordered = T)



## boxplot fraction abundance ####

cols_fraction <- c("#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf","#2166ac")

plot2 <- ggplot(genome_percentage_all, aes(x=new.order, y=fraction, fill = class)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=1,
               outlier.size=0.9, linewidth = 0.2)+
  theme_bw()+
  theme(plot.title = element_text(face="bold", size=7),
        axis.title.y = element_text(face="bold", size=7),
        axis.title.x = element_text(face="bold", size=7),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7),
        axis.text.y = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill='white', color=NA),
        legend.background = element_rect(fill='white'),
        legend.position = "right")+
  scale_y_sqrt(expand = c(0,0),limits= c(0,50),breaks = c(0,1,5,10,20,30,40))+
  scale_fill_manual(values = cols_fraction)+
  labs(x = "",
       y = "fraction [%] (sqrt)",
       fill = "class") 

ggsave("boxplot_fraction_samples.png", width = 18, height = 16, units = "cm", bg='white')
ggsave("boxplot_fraction_samples.pdf", width = 18, height = 16, units = "cm", bg='white')


plot2 /plot1 + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 7, face = "bold"))

ggsave("../../CoverM/v3/boxplot_fraction_abundance_orders.png", width = 18, height = 20, units = "cm", bg='white')
ggsave("../../CoverM/v3/boxplot_fraction_abundance_orders.pdf", width = 18, height = 20, units = "cm", bg='white')

