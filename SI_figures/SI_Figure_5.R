library(reshape)
library(tidyverse)
library(ggplot2)

setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Data_mining/fastANI/fastANI_EX_MAG_v3/")

# load fastANI data ####
fastANI_data <- read.table("EX_fastANI_out.txt", sep="\t", header=FALSE, col.names=c("target", "query", "ANI"))

fastANI_data$target <- sapply(strsplit(fastANI_data$target, "/\\s*"), tail, 1)
fastANI_data$target <- gsub(".fasta","",fastANI_data$target)
fastANI_data$target <- gsub(".fa","",fastANI_data$target)

fastANI_data$query <- sapply(strsplit(fastANI_data$query, "/\\s*"), tail, 1)
fastANI_data$query <- gsub(".fasta","",fastANI_data$query)
fastANI_data$query <- gsub(".fa","",fastANI_data$query)

# create Matrix ####
# extract unique target genomes
target_genomes <- unique(fastANI_data$target)

# extract unique query genomes
query_genomes <- unique(fastANI_data$query)


# Create an empty matrix with dimensions based on the unique genomes
ANI_matrix <- matrix(nrow = length(target_genomes), ncol = length(query_genomes),
                     dimnames = list(target_genomes, query_genomes))

# Fill in the matrix with average nucleotide identity values
for (i in 1:nrow(fastANI_data)) {
  target_genome <- fastANI_data[i, "target"]
  query_genome <- fastANI_data[i, "query"]
  ani_value <- fastANI_data[i, "ANI"]
  ANI_matrix[target_genome, query_genome] <- ani_value
}

# Convert the matrix to a data frame for tabular representation
#ANI_matrix_df <- as.data.frame(ANI_matrix)

ANI_matrix_df <- melt.matrix(ANI_matrix)

colnames(ANI_matrix_df) <- c("target","query","ANI")


## Fill the whole matrix ####
# to fill the whole matrix the table is also reversed and target and query changed and then merged with initial table
reversed_table <- ANI_matrix_df
reversed_table <- reversed_table[c("target", "query", "ANI")]
colnames(reversed_table) <- c("query", "target", "ANI")
merged_table <- merge(ANI_matrix_df, reversed_table, by = c("query", "target"), all = TRUE)

# Combine the separate AAI columns into a single AAI column
merged_table$ANI <- ifelse(is.na(merged_table$ANI.x), merged_table$ANI.y, merged_table$ANI.x)

# Remove the separate AAI columns if needed
merged_table <- subset(merged_table, select = -c(ANI.x, ANI.y))



# reorder the levels
merged_table_1 <- merged_table

levels(merged_table_1$query)
levels(merged_table_1$target)
merged_table_1$query <- factor(merged_table_1$query, levels = c("SUTE22-1_SAMN10231914_MAG_00000186","SUTE22-1_SAMN10231894_MAG_00000173","SUTE22-1_SAMN10231913_MAG_00000219","SUTE22-1_SAMN10231893_MAG_00000067",
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
                                          "EMB267","PRJNA368391_1","GCA_003663595.1","MSM105","PRJNA541421_1","PRJNA889212_2","E3_1_d157"), ordered = T)

merged_table_1$target <- factor(merged_table_1$target, levels = c("SUTE22-1_SAMN10231914_MAG_00000186","SUTE22-1_SAMN10231894_MAG_00000173","SUTE22-1_SAMN10231913_MAG_00000219","SUTE22-1_SAMN10231893_MAG_00000067",
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
                                           "EMB267","PRJNA368391_1","GCA_003663595.1","MSM105","PRJNA541421_1","PRJNA889212_2","E3_1_d157"), ordered = T)


## plot heatmap ####
ggplot(merged_table_1, aes(x = target, y = query, fill = ANI)) +
  geom_tile(color = "#f0f0f0") +
  scale_fill_gradient2(low = "#ffeda0",
                       mid = "#a6dba0",
                       high = "#053061",
                       midpoint = 88,
                       na.value = "white") +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7, face = "bold"))

ggsave("FastANI_EX_MAGs_Heatmap.png", plot = last_plot(), width = 15, height = 15, units = "cm", device = "png" )


# Only triangle heatmap ####
merged_table_2 <- merged_table

merged_table_2$query <- factor(merged_table_2$query, levels = c("SUTE22-1_SAMN10231914_MAG_00000186","SUTE22-1_SAMN10231894_MAG_00000173","SUTE22-1_SAMN10231913_MAG_00000219","SUTE22-1_SAMN10231893_MAG_00000067",
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
                                          "EMB267","PRJNA368391_1","GCA_003663595.1","MSM105","PRJNA541421_1","PRJNA889212_2","E3_1_d157"), ordered = T)

merged_table_2$target <- factor(merged_table_2$target, levels = c("SUTE22-1_SAMN10231914_MAG_00000186","SUTE22-1_SAMN10231894_MAG_00000173","SUTE22-1_SAMN10231913_MAG_00000219","SUTE22-1_SAMN10231893_MAG_00000067",
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
                                           "EMB267","PRJNA368391_1","GCA_003663595.1","MSM105","PRJNA541421_1","PRJNA889212_2","E3_1_d157"), ordered = T)

# extract unique target genomes
target_genomes_merged <- unique(merged_table_2$target)
# extract unique query genomes
query_genomes_merged <- unique(merged_table_2$query)

merged_table_matrix <- matrix(nrow = length(target_genomes_merged), ncol = length(query_genomes_merged),
                              dimnames = list(target_genomes_merged, query_genomes_merged))



new_order_tree = c("SAMN10231914_1", "SAMN10231894_1","SAMN10231913_1", "SAMN10231893_1",
                   "SAMN10231914_2","SAMN10231913_2", 
                   "SAMN10231904_1","SAMN11854494_1","PRJNA704804_1","SAMEA2620113_1",
                   "PRJNA541421_3","PRJNA531756_1",
                   "GCA_021160765.1","GCA_021158585.1",
                   "GCA_003663625.1","PRJNA368391_2","GCA_003650025.1","GCA_011039765.1",
                   "SAMN22703512_1","SAMN22703514_1","SAMN22703513_1",
                   "SAMN30647027_1","GCA_023145185.1",
                   "PRJNA531756_2","GCA_016928095.1","PRJNA889212_1","PRJNA541421_2","PRJNA721298_1",
                   "EMB267","PRJNA368391_1","GCA_003663595.1","MSM105","PRJNA541421_1","PRJNA889212_2","E3_1_d157")


valid_names <- intersect(new_order_tree, rownames(merged_table_matrix))
merged_table_matrix <- merged_table_matrix[valid_names, valid_names]


# fill in the matrix
for (i in 1:nrow(merged_table_2)) {
  target_genome <- merged_table_2[i, "target"]
  query_genome <- merged_table_2[i, "query"]
  ani_value <- merged_table_2[i, "ANI"]
  merged_table_matrix[target_genome, query_genome] <- ani_value
}

# Get upper triangle of the matrix
get_upper_tri <- function(merged_table_matrix){
  merged_table_matrix[lower.tri(merged_table_matrix)]<- NA
  return(merged_table_matrix)
}

upper_tri <- get_upper_tri(merged_table_matrix)  

     

# Convert the matrix to a data frame for tabular representation
upper_tri_df <- melt(upper_tri)
colnames(upper_tri_df) <- c("target","query","ANI")

# plot heatmap
ggplot(upper_tri_df, aes(x = target, y = query, fill = ANI)) +
  geom_tile(color = "#f0f0f0") +
  geom_text(aes(label = round(ANI)), size = 2.5) +
  scale_fill_gradientn("ANI [%]",colours = c("#5aae61","#d9f0d3","#e7d4e8","#9970ab","#40004b"),
                       na.value = "white",
                       limits = c(75,100)) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7, face = "bold"))

ggsave("FastANI_EX_MAGs_Heatmap_upper_tri.png", plot = last_plot(), width = 18, height = 18, units = "cm", device = "png" )

save.image("ANI_heatmap.RData")
