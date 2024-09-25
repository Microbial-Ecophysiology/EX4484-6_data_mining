
library(reshape)
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)
library(patchwork)

setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Data_mining/Orthofinder/NMDS")

load("gene_cluster_analysis.RData")

### load data obtained from orthofinder ####
Orthogroup_gene_count <- read.table("../Orthogroups/Orthogroups.GeneCount.tsv", header = T, sep="\t", comment.char = "",
                                    quote = "\"",
                                    check.names = FALSE)
Orthogroup_gene_count <-select(Orthogroup_gene_count, -Total)
rownames(Orthogroup_gene_count) <- Orthogroup_gene_count[,1]
Orthogroup_gene_count[,1] <- NULL

# transpose dataframe 
Orthogroup_gene_count_t <- t(Orthogroup_gene_count)

### create a presence absence matrix by converting each field with gene cluster presence >2 to 1 ####
Orthogroup_gene_count_t[Orthogroup_gene_count_t > 0] <- 1


#  All orthogroups ####

orthogroups <- colnames(Orthogroup_gene_count_t)

nmds_all <- metaMDS(Orthogroup_gene_count_t, dist="jaccard", k=2, trymax= 200) # Best solution repeated 2 times
test_all <- as.data.frame(scores(nmds_all, "sites"))
test_all$Genome <- row.names(Orthogroup_gene_count_t)


## plot NMDS ####
ggplot(test_all, aes(NMDS1, NMDS2, label=Genome)) + geom_text_repel()

test_all$Genome <- gsub(".genes","",test_all$Genome)
test_all$Genome

# color based on new taxonomy
new_family <- read.table("new_family_groups.txt", sep="\t", header=T)
names(new_family)[1] <- "Genome"
nmds_plot_family_all <- merge(test_all, new_family, by="Genome")

# color based on habitat
habitat <- read.table("GTDB_classification_habitat.txt", sep="\t", header=F)
names(habitat)[1] <- "Genome"
nmds_plot_habitat_all <- merge(test_all, habitat, by="Genome")

save.image("gene_cluster_analysis.RData")
load("gene_cluster_analysis.RData")


nmds_plot_family_all$Genome <- factor(nmds_plot_family_all$Genome, levels = c("E3_1_d157_spades_bin_16_orig_1_refined-contigs",
                                                                "GCA_003663595.1_ASM366359v1_genomic",
                                                                "GCA_003650025.1_ASM365002v1_genomic",
                                                                "GCA_003663625.1_ASM366362v1_genomic",
                                                                "GCA_011039765.1_ASM1103976v1_genomic",
                                                                "GCA_016928095.1_ASM1692809v1_genomic",
                                                                "GCA_023145185.1_ASM2314518v1_genomic",
                                                                "GCA_021158585.1_ASM2115858v1_genomic",
                                                                "GCA_021160765.1_ASM2116076v1_genomic",
                                                                "PRJNA368391_bin_183_orig-contigs",
                                                                "PRJNA368391_bin_90_orig-contigs",
                                                                "PRJNA531756_bin_72_orig_refined-contigs",
                                                                "PRJNA531756_bin_93_strict_refined-contigs",
                                                                "PRJNA541421_bin_70_strict-contigs",
                                                                "PRJNA541421_bin_72_orig_refined-contigs",
                                                                "PRJNA541421_bin_98_orig_refined-contigs",
                                                                "PRJNA704804_bin_147_orig_refined-contigs",
                                                                # "PRJNA713414_bin_118_orig_1_refined-contigs",
                                                                # "PRJNA713414_bin_178_orig_1_refined-contigs",
                                                                "PRJNA721298_bin_46_orig_refined-contigs",
                                                                "PRJNA889212_bin_1_orig_refined-contigs",
                                                                "PRJNA889212_bin_7_strict_refined-contigs",
                                                                "MSM105_N25025F_bin_277_ori_permissive_1",
                                                                "EMB267_Co_bin_434_strict_1",
                                                                "SCRA20-1_SAMN11854494_MAG_00000079",
                                                                "SUTE22-1_SAMN10231893_MAG_00000067",
                                                                "SUTE22-1_SAMN10231913_MAG_00000219",
                                                                "SUTE22-1_SAMN10231913_MAG_00000177",
                                                                "SUTE22-1_SAMN10231914_MAG_00000186",
                                                                "SUTE22-1_SAMN10231914_MAG_00000292",
                                                                "SUTE22-1_SAMN10231894_MAG_00000173",
                                                                "SUTE22-1_SAMN10231904_MAG_00000012",
                                                                "TARA_SAMEA2620113_MAG_00000097",
                                                                "ZHEN22-1_SAMN22703512_MAG_00000312",
                                                                "ZHEN22-1_SAMN22703513_MAG_00000265",
                                                                "ZHEN22-1_SAMN22703514_MAG_00000188",
                                                                "ZORZ22-1_SAMN30647027_MAG_00000060"),
                               labels = c("E3_1_d157",
                                          "GCA_003663595.1",
                                          "GCA_003650025.1",
                                          "GCA_003663625.1",
                                          "GCA_011039765.1",
                                          "GCA_016928095.1",
                                          "GCA_023145185.1",
                                          "GCA_021158585.1",
                                          "GCA_021160765.1",
                                          "PRJNA368391_1",
                                          "PRJNA368391_2",
                                          "PRJNA531756_1",
                                          "PRJNA531756_2",
                                          "PRJNA541421_1",
                                          "PRJNA541421_2",
                                          "PRJNA541421_3",
                                          "PRJNA704804_1",
                                          # "PRJNA713414_1",
                                          # "PRJNA713414_2",
                                          "PRJNA721298_1",
                                          "PRJNA889212_1",
                                          "PRJNA889212_2",
                                          "MSM105",
                                          "EMB267",
                                          "SAMN11854494",
                                          "SAMN10231893",
                                          "SAMN10231913_1",
                                          "SAMN10231913_2",
                                          "SAMN10231914_1",
                                          "SAMN10231914_2",
                                          "SAMN10231894",
                                          "SAMN10231904",
                                          "SAMEA2620113",
                                          "SAMN22703512",
                                          "SAMN22703513",
                                          "SAMN22703514",
                                          "SAMN30647027"), ordered = TRUE)


nmds_plot_habitat_all$Genome <- factor(nmds_plot_habitat_all$Genome, levels = c("E3_1_d157_spades_bin_16_orig_1_refined-contigs",
                                                                              "GCA_003663595.1_ASM366359v1_genomic",
                                                                              "GCA_003650025.1_ASM365002v1_genomic",
                                                                              "GCA_003663625.1_ASM366362v1_genomic",
                                                                              "GCA_011039765.1_ASM1103976v1_genomic",
                                                                              "GCA_016928095.1_ASM1692809v1_genomic",
                                                                              "GCA_023145185.1_ASM2314518v1_genomic",
                                                                              "GCA_021158585.1_ASM2115858v1_genomic",
                                                                              "GCA_021160765.1_ASM2116076v1_genomic",
                                                                              "PRJNA368391_bin_183_orig-contigs",
                                                                              "PRJNA368391_bin_90_orig-contigs",
                                                                              "PRJNA531756_bin_72_orig_refined-contigs",
                                                                              "PRJNA531756_bin_93_strict_refined-contigs",
                                                                              "PRJNA541421_bin_70_strict-contigs",
                                                                              "PRJNA541421_bin_72_orig_refined-contigs",
                                                                              "PRJNA541421_bin_98_orig_refined-contigs",
                                                                              "PRJNA704804_bin_147_orig_refined-contigs",
                                                                              # "PRJNA713414_bin_118_orig_1_refined-contigs",
                                                                              # "PRJNA713414_bin_178_orig_1_refined-contigs",
                                                                              "PRJNA721298_bin_46_orig_refined-contigs",
                                                                              "PRJNA889212_bin_1_orig_refined-contigs",
                                                                              "PRJNA889212_bin_7_strict_refined-contigs",
                                                                              "MSM105_N25025F_bin_277_ori_permissive_1",
                                                                              "EMB267_Co_bin_434_strict_1",
                                                                              "SCRA20-1_SAMN11854494_MAG_00000079",
                                                                              "SUTE22-1_SAMN10231893_MAG_00000067",
                                                                              "SUTE22-1_SAMN10231913_MAG_00000219",
                                                                              "SUTE22-1_SAMN10231913_MAG_00000177",
                                                                              "SUTE22-1_SAMN10231914_MAG_00000186",
                                                                              "SUTE22-1_SAMN10231914_MAG_00000292",
                                                                              "SUTE22-1_SAMN10231894_MAG_00000173",
                                                                              "SUTE22-1_SAMN10231904_MAG_00000012",
                                                                              "TARA_SAMEA2620113_MAG_00000097",
                                                                              "ZHEN22-1_SAMN22703512_MAG_00000312",
                                                                              "ZHEN22-1_SAMN22703513_MAG_00000265",
                                                                              "ZHEN22-1_SAMN22703514_MAG_00000188",
                                                                              "ZORZ22-1_SAMN30647027_MAG_00000060"),
                                      labels = c("E3_1_d157",
                                                 "GCA_003663595.1",
                                                 "GCA_003650025.1",
                                                 "GCA_003663625.1",
                                                 "GCA_011039765.1",
                                                 "GCA_016928095.1",
                                                 "GCA_023145185.1",
                                                 "GCA_021158585.1",
                                                 "GCA_021160765.1",
                                                 "PRJNA368391_1",
                                                 "PRJNA368391_2",
                                                 "PRJNA531756_1",
                                                 "PRJNA531756_2",
                                                 "PRJNA541421_1",
                                                 "PRJNA541421_2",
                                                 "PRJNA541421_3",
                                                 "PRJNA704804_1",
                                                 # "PRJNA713414_1",
                                                 # "PRJNA713414_2",
                                                 "PRJNA721298_1",
                                                 "PRJNA889212_1",
                                                 "PRJNA889212_2",
                                                 "MSM105",
                                                 "EMB267",
                                                 "SAMN11854494",
                                                 "SAMN10231893",
                                                 "SAMN10231913_1",
                                                 "SAMN10231913_2",
                                                 "SAMN10231914_1",
                                                 "SAMN10231914_2",
                                                 "SAMN10231894",
                                                 "SAMN10231904",
                                                 "SAMEA2620113",
                                                 "SAMN22703512",
                                                 "SAMN22703513",
                                                 "SAMN22703514",
                                                 "SAMN30647027"), ordered = TRUE)



### plot NMDS, using GTDB classification on family level for color coding ####
### NMDS single plot 
ggplot(nmds_plot_family_all, aes(NMDS1, NMDS2, col=Family)) + 
  geom_point(size = 1) + 
  theme_bw() +
  scale_color_manual(values = c("Family 1A" = "#c7e9c0", "Family 1B" = "#74c476","Family 2" = "#006d2c", "Family 3A" = "#08306b", "Family 3B" = "#2171b5", "Family 4" = "#9ecae1"))+
  theme(axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        axis.title.y = element_text(face="bold", size=7),
        axis.title.x = element_text(face="bold", size=7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=7),
        legend.title=element_text(size=7, face="bold"),
        legend.key.size = unit(0.35, "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size = 7, face="bold"))+
  geom_text_repel(label=nmds_plot_family_all$Genome,
                  size = 2,
                  point.size = 2,
                  lineheight = 0.5,
                  show.legend = F,
                  max.overlaps = Inf)+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

ggsave("NMDS_single.png", width = 18, height = 15, units = "cm")
ggsave("NMDS_single.pdf", width = 18, height = 15, units = "cm")


# alternatively by habitat ####
cols <- c("#000000","#e31a1c","#d9d9d9","#74add1","#07316A")
ggplot(nmds_plot_habitat_all, aes(NMDS1, NMDS2, col=V8)) + 
  geom_point(size = 1) + 
  theme_bw() +
  scale_color_manual("environment", values = cols)+
  theme(axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        axis.title.y = element_text(face="bold", size=7),
        axis.title.x = element_text(face="bold", size=7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=7),
        legend.title=element_text(size=7, face="bold"),
        legend.key.size = unit(0.35, "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size = 7, face="bold"))+
  geom_text_repel(label=nmds_plot_habitat_all$Genome,
                  size = 2,
                  point.size = 2,
                  lineheight = 0.5,
                  show.legend = F,
                  max.overlaps = Inf)+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

ggsave("NMDS_single_habitat.png", width = 18, height = 15, units = "cm")
ggsave("NMDS_single_habitat.pdf", width = 18, height = 15, units = "cm")
