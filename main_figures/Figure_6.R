# Combine percentage of unknown genes with relative abundance data

library(tidyverse)
library(ggplot2)
library(reshape2)
library(patchwork)
library(ggpubr)

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
df_all_subset <- df_all[df_all$rel_abundance != 0, ]


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




genome_percentage_all$class <- factor(genome_percentage_all$class, levels = c("E2", "EX4484-6","Poseidoniia","SW-10-69-26","Thermoplasmata_A","Thermoplasmata","UBA186"), ordered = T)

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
                                                     "Thermoplasmatales (N=64)","ARK-15 (N=3)","B87-G9 (N=1)","PWKY01 (N=12)","SG8-5 (N=7)","Thermoplasmata unclassified (N=1)",
                                                     "DTKX01 (N=2)",
                                                     "UBA186 (N=1)","UBA287 (N=1)"),ordered = T)


# ## boxplot fraction abundance ####
# 
cols_fraction <- c("#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf","#2166ac")
# 
# fraction_abundance <- ggplot(genome_percentage_all, aes(x=new.order, y=fraction, fill = class)) + 
#   geom_boxplot(outlier.colour="black", outlier.shape=1,
#                outlier.size=0.9, linewidth = 0.2)+
#   theme_bw()+
#   theme(plot.title = element_text(face="bold", size=7),
#         axis.title.y = element_text(face="bold", size=7),
#         axis.title.x = element_text(face="bold", size=7),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5),
#         axis.text.y = element_text(size = 7),
#         legend.text = element_text(size = 7),
#         legend.title = element_text(size = 7, face = "bold"),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         strip.background = element_blank(),
#         panel.background = element_rect(fill='white'),
#         plot.background = element_rect(fill='white', color=NA),
#         legend.background = element_rect(fill='white'),
#         legend.position = "none")+
#   scale_y_sqrt(breaks = c(0,1,5,10,20,30,40), expand = c(0,0), limits = c(0, 40))+
#   scale_fill_manual(values = cols_fraction)+
#   labs(x = "",
#        y = "fraction [%] (sqrt)",
#        fill = "class") 
# 
# ggsave("boxplot_fraction_samples.png", width = 18, height = 16, units = "cm", bg='white')
# ggsave("boxplot_fraction_samples.pdf", width = 18, height = 16, units = "cm", bg='white')


n_observations_percentage_unknown <- genome_percentage_all %>% count(new.order)

## boxplot percentage unknown ####
percentage_unknown <- ggplot(genome_percentage_all, aes(x=new.order, y=percentage_unknown, fill = class)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=1,
               outlier.size=0.9, linewidth = 0.2)+
  theme_bw()+
  theme(plot.title = element_text(face="bold", size=7),
        axis.title.y = element_text(face="bold", size=7, color = "black"),
        axis.title.x = element_text(face="bold", size=7, color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7, color = "black"),
        legend.title = element_text(size = 7, face = "bold", color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill='white', color=NA),
        legend.background = element_rect(fill='white'),
        legend.position = "none")+
  scale_y_continuous(expand= c(0,0), limits = c(0,79),breaks = scales::pretty_breaks(n = 10))+
  scale_fill_manual(values = cols_fraction)+
  labs(x = "",
       y = "unknown genes [%]",
       fill = "class") 
ggsave("boxplot_unknown_genes_samples.png", width = 18, height = 16, units = "cm", bg='white')
ggsave("boxplot_unknown_genes_samples.pdf", width = 18, height = 16, units = "cm", bg='white')


# define rarity #### 
df_matrix <- select(df_all_subset, c("genome","accession","rel_abundance"))

mat_df <- reshape::cast(df_matrix, genome~accession, value = "rel_abundance", fun.aggregate = "mean", fill = 0)
rownames(mat_df) <- mat_df$genome
mat <- mat_df[, -1]
genome_rarity <- apply(mat, 1, function(x) median(x[x > 0])) * apply(mat, 1, function(x) sum(x > 0)/length(x)) # rarity as the median of relative abundances * fraction in screened samples
names(genome_rarity) <- rownames(mat)
hist(genome_rarity)
all.equal(genome_percentage$genome, names(genome_rarity))
genome_percentage$rarity <- genome_rarity[match(genome_percentage$genome, names(genome_rarity))]
plot(sqrt(genome_percentage$rarity), genome_percentage$percentage_unknown, col = as.numeric(as.factor(genome_percentage$class)), pch = 16)

## add more data about environments of single MAGs ####
habitat <- read.csv("GTDB_classification_habitat.txt", 
                    sep="\t", 
                    header=TRUE, 
                    fill=TRUE, 
                    quote="")

habitat_subset <- select(habitat, c("user_genome","habitat","group","overview"))

genome_percentage_habitat <- full_join(x= genome_percentage, y = habitat_subset, by = c("genome" = "user_genome"))
write.table(genome_percentage_habitat, "genome_percentage_habitat.txt", 
            quote = FALSE,
            sep="\t",
            col.names = TRUE,
            row.names = FALSE)

unique(genome_percentage_habitat$group) # 18 different groups


# boxplot with rarity categories ####
## combine with habitat
## using the median of rarity as threshold for rare vs common species
plot(sqrt(genome_percentage_habitat$rarity), genome_percentage_habitat$percentage_unknown, col = as.numeric(as.factor(genome_percentage_habitat$class)), pch = 16)
abline(v = sqrt(median(genome_percentage_habitat$rarity, na.rm = T)))
genome_percentage_habitat$rarity_group <- ifelse(genome_percentage_habitat$rarity <= median(genome_percentage_habitat$rarity, na.rm = T), "rare", "common")

#genomes with no presence in our screened dataset were colelcted as not detected (nd)
genome_percentage_habitat$rarity_group[is.na(genome_percentage_habitat$rarity_group)] <- "nd"
table(genome_percentage_habitat$rarity_group)

#subset data frame based on genome, group (habitat) and rarity_group
overview_rarity_group <- table(genome_percentage_habitat$group, genome_percentage_habitat$rarity_group)
overview_rarity_group <- as.data.frame.matrix(overview_rarity_group)

write.table(overview_rarity_group, "habitat_overview_rarity_group.txt", 
            quote = FALSE,
            sep="\t",
            col.names = TRUE,
            row.names = TRUE)


## for plotting
## using the median of rarity as threshold for rare vs common species
plot(sqrt(genome_percentage$rarity), genome_percentage$percentage_unknown, col = as.numeric(as.factor(genome_percentage$class)), pch = 16)
abline(v = sqrt(median(genome_percentage$rarity, na.rm = T)))
genome_percentage$rarity_group <- ifelse(genome_percentage$rarity <= median(genome_percentage$rarity, na.rm = T), "rare", "common")

#genomes with no presence in our screened dataset were colelcted as not detected (nd)
genome_percentage$rarity_group[is.na(genome_percentage$rarity_group)] <- "nd"
table(genome_percentage$rarity_group)


## plot unknown genes vs rarity ####

# ggplot does not transform the data as squared, but only displays the axis as squared with scale_x_sqrt
intercept_rarity_threshold <- median(genome_percentage$rarity, na.rm = T)

genome_percentage$class <- factor(genome_percentage$class, levels = c("E2","EX4484-6","Poseidoniia","SW-10-69-26","Thermoplasmata","UBA186"), ordered = T)
genome_percentage$rarity_group <- factor(genome_percentage$rarity_group, levels = c("rare","common","nd"),
                                         labels = c("rare \n(N=137)","common \n(N=136)","nd \n(N=132)"),ordered = T)


n_observations_genome_percentage <- genome_percentage %>% count(rarity_group)

## plot class only ####
rarity_plot <- ggplot()+
  geom_point(data = genome_percentage,
             aes(x = rarity, y = percentage_unknown, color = class), shape = 16, size = 2)+
  labs(x = "rarity (sqrt)",
       y = "unknown genes [%]",
       fill = "environment") +
  theme_bw() +  
  scale_color_manual("class", values = cols_fraction,drop = F) +
  theme(legend.position = "right")+
  theme(plot.title = element_text(face="bold", size=7),
        axis.title.y = element_text(face="bold", size=7, color = "black"),
        axis.title.x = element_text(face="bold", size=7, color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7, color = "black"),
        legend.title = element_text(size = 7, face = "bold", color = "black"),
        legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill='white', color=NA),
        legend.background = element_rect(fill='white'))+
  scale_y_continuous(expand= c(0,0), limits = c(0,79),breaks = scales::pretty_breaks(n = 10))+
  scale_x_sqrt(expand= c(0,0), breaks = c(0.001,0.005, 0.01,0.015, 0.02, 0.025, 0.1, 0.15), limits =c(0,0.035))+
  geom_vline(xintercept = intercept_rarity_threshold, linewidth = 0.3)+
  annotate("text", x=0.00013, y=65, label="median rarity", angle=90, size = 2.5)+
  guides(col = guide_legend(override.aes = list(size=2)))

ggsave("rarity_vs_unknown_genes.png", width = 18, height = 16, units = "cm", bg='white')


my_comparisons <- list( c("rare \n(N=137)", "common \n(N=136)"), c("rare \n(N=137)", "nd \n(N=132)"), c("common \n(N=136)", "nd \n(N=132)") )

##boxplot using rarity_group ####
box_rarity <- ggplot(genome_percentage, aes(x=rarity_group, y=percentage_unknown, fill = rarity_group)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=1,
               outlier.size=1, linewidth = 0.3, fatten = 1)+
  theme_bw()+
  theme(plot.title = element_text(face="bold", size=7),
        axis.title.y = element_text(face="bold", size=7, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7, face = "bold"),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill='white', color=NA),
        legend.background = element_rect(fill='white'))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), expand= c(0,0), limits=c(0,79))+
  scale_fill_brewer("rarity group", palette="Greys")+
  labs(x = "rarity group",
       y = "unknown genes [%]",
       fill = "rarity_group") + 
  stat_compare_means(comparisons = my_comparisons, label.y = c(62, 71, 65),
                     label = "p.signif", method = "wilcox.test")


# boxplot of percentage of unknown genes sorted by rarity group - quick test
## the boxplot and following t-tests show significant differences between the single groups
boxplot(genome_percentage$percentage_unknown ~ genome_percentage$rarity_group)
wilcox.test(percentage_unknown ~ rarity_group, data = genome_percentage[genome_percentage$rarity_group != "nd", ])
wilcox.test(percentage_unknown ~ rarity_group, data = genome_percentage[genome_percentage$rarity_group != "rare", ])
wilcox.test(percentage_unknown ~ rarity_group, data = genome_percentage[genome_percentage$rarity_group != "common", ])
0.05/3 # 0.01666667
# end




plot_2 <- rarity_plot + inset_element(box_rarity, left = 0.3, bottom = 0.3, right = 0.99, top = 0.99) + 
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 7, face = "bold"))

ggsave("Figure6_manuscript.png", width = 18, height = 12, units = "cm", bg='white')


# combined plot ####
combined_plot <- (plot_2 / percentage_unknown) +
  plot_layout(guides = "collect")+
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 7, face = "bold"),
        legend.justification = "left") 

ggsave(plot =  combined_plot, "Figure6_ms_mod.png", width = 18, height = 24, units = "cm", bg='white')
ggsave(plot =  combined_plot, "Figure6_ms_mod.emf", width = 18, height = 24, units = "cm", bg='white')
ggsave(plot =  combined_plot, "Figure6_ms_mod.pdf", width = 18, height = 24, units = "cm", bg='white')



save.image("rarity_unknown_genes.RData")
load("rarity_unknown_genes.RData")
