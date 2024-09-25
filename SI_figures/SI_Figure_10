# Script to analyze shared Orthogroups retrieved through Orthofinder

library(reshape)
library(ggplot2)

setwd("C:/Users/admin/OneDrive/PhD/Data_mining/Orthofinder/")

# Comparative statistics ####
## Orthogroups Species Overlap

Species_overlap <- read.table("Comparative_Genomics_Statistics/Orthogroups_SpeciesOverlaps.tsv", header = T, sep="\t", comment.char = "",
                         quote = "\"",
                         check.names = FALSE)

rownames(Species_overlap) <- Species_overlap[,1]
Species_overlap[,1] <- NULL


# Function to calculate percentage of shared orthologs relative to the reference species
calculate_percentage_relative_to_reference <- function(species_column, reference_species_total_ortho) {
  return(species_column / reference_species_total_ortho * 100)
}

# Get the names of species
species_names <- colnames(Species_overlap)

# Initialize an empty matrix to store the results
result_matrix <- matrix(NA, nrow = length(species_names), ncol = length(species_names), dimnames = list(species_names, species_names))

# Loop through each species as the reference species
for (reference_species in species_names) {
  # Find the index of the reference species in the column names
  reference_species_index <- which(species_names == reference_species)
  
  # Calculate percentage of shared orthologs for each species relative to the reference species
  percentage_relative_to_reference <- calculate_percentage_relative_to_reference(Species_overlap[, reference_species_index], Species_overlap[reference_species, reference_species])
  
  # Store the percentages in the result matrix for both row and column
  result_matrix[, reference_species] <- percentage_relative_to_reference
  result_matrix[reference_species, ] <- percentage_relative_to_reference
}

# View the resulting matrix
print(result_matrix)

orthogroup_shared_df <- melt(result_matrix, varnames = c("Species", "Reference_Species"))

# gives warning
# Warnmeldungen:
#   1: In type.convert.default(X[[i]], ...) :
#   'as.is' should be specified by the caller; using TRUE
# 2: In type.convert.default(X[[i]], ...) :
#   'as.is' should be specified by the caller; using TRUE

orthogroup_shared_df$Species <- factor(orthogroup_shared_df$Species, levels = c("SUTE22-1_SAMN10231914_MAG_00000186.genes","SUTE22-1_SAMN10231894_MAG_00000173.genes","SUTE22-1_SAMN10231913_MAG_00000219.genes","SUTE22-1_SAMN10231893_MAG_00000067.genes",
                                                                            "SUTE22-1_SAMN10231914_MAG_00000292.genes","SUTE22-1_SAMN10231913_MAG_00000177.genes",
                                                                            "SUTE22-1_SAMN10231904_MAG_00000012.genes","SCRA20-1_SAMN11854494_MAG_00000079.genes","PRJNA704804_bin_147_orig_refined-contigs","TARA_SAMEA2620113_MAG_00000097.genes",
                                                                            "PRJNA541421_bin_98_orig_refined-contigs","PRJNA531756_bin_72_orig_refined-contigs",
                                                                            "GCA_021160765.1_ASM2116076v1_genomic","GCA_021158585.1_ASM2115858v1_genomic",
                                                                            "GCA_003663625.1_ASM366362v1_genomic","PRJNA368391_bin_90_orig-contigs","GCA_003650025.1_ASM365002v1_genomic","GCA_011039765.1_ASM1103976v1_genomic",
                                                                            "ZHEN22-1_SAMN22703512_MAG_00000312.genes","ZHEN22-1_SAMN22703514_MAG_00000188.genes","ZHEN22-1_SAMN22703513_MAG_00000265.genes",
                                                                            "ZORZ22-1_SAMN30647027_MAG_00000060.genes","GCA_023145185.1_ASM2314518v1_genomic", 
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

orthogroup_shared_df$Reference_Species <- factor(orthogroup_shared_df$Reference_Species, levels = c("SUTE22-1_SAMN10231914_MAG_00000186.genes","SUTE22-1_SAMN10231894_MAG_00000173.genes","SUTE22-1_SAMN10231913_MAG_00000219.genes","SUTE22-1_SAMN10231893_MAG_00000067.genes",
                                                                            "SUTE22-1_SAMN10231914_MAG_00000292.genes","SUTE22-1_SAMN10231913_MAG_00000177.genes",
                                                                            "SUTE22-1_SAMN10231904_MAG_00000012.genes","SCRA20-1_SAMN11854494_MAG_00000079.genes","PRJNA704804_bin_147_orig_refined-contigs","TARA_SAMEA2620113_MAG_00000097.genes",
                                                                            "PRJNA541421_bin_98_orig_refined-contigs","PRJNA531756_bin_72_orig_refined-contigs",
                                                                            "GCA_021160765.1_ASM2116076v1_genomic","GCA_021158585.1_ASM2115858v1_genomic",
                                                                            "GCA_003663625.1_ASM366362v1_genomic","PRJNA368391_bin_90_orig-contigs","GCA_003650025.1_ASM365002v1_genomic","GCA_011039765.1_ASM1103976v1_genomic",
                                                                            "ZHEN22-1_SAMN22703512_MAG_00000312.genes","ZHEN22-1_SAMN22703514_MAG_00000188.genes","ZHEN22-1_SAMN22703513_MAG_00000265.genes",
                                                                            "ZORZ22-1_SAMN30647027_MAG_00000060.genes","GCA_023145185.1_ASM2314518v1_genomic", 
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

## plot heatmap of the percentages of shared otrthogroups ####

ggplot(orthogroup_shared_df, aes(x = Species, y = Reference_Species, fill = value)) +
  geom_tile(color = "#f0f0f0") +
  geom_text(aes(label = round(value)), size = 2.5) +
  scale_fill_gradientn("percentage shared \northogroups [%]",colours = c("#40004b","#9970ab","#e7d4e8","white", "#d9f0d3","#5aae61","#00441b"),
                       breaks = c(40, 60, 80, 100), labels = c("40", "60",  "80",  "100")) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1, size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7, face = "bold"))


ggsave("shared_orthogroups.png", plot = last_plot(), width = 18, height = 18, units = "cm", device = "png" )
ggsave("shared_orthogroups.pdf", plot = last_plot(), width = 18, height = 18, units = "cm", device = "pdf" )


