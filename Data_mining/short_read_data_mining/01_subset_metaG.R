# Script to combine and search a table for marine metagenome samples

# From ENA we downloaded a list containing all ecological metagenomes (total 44968 MGs):
# tax_tree(410657) AND library_strategy="WGS" AND instrument_platform="ILLUMINA" AND library_source="METAGENOMIC" AND library_selection="RANDOM"

library(dplyr)
require(marmap)
require(maps)
require(R.utils)
require(tidyverse)

setwd("C:/Users/admin/OneDrive/PhD/Data_mining")

Metagenome <- read.table(
  "results_MG_mining.txt", 
  sep="\t",
  header=T, 
  comment.char = "",
  quote = "",
  check.names = FALSE
)
str(Metagenome)

# get a list of all unique values
listMGs <- sort(unique(Metagenome$scientific_name))

write.table(listMGs, 'listMGs.txt', row.names = F, col.names = F, sep = "\t", quote = F)
listMGs_subsampled <- scan(what = "character", file = "listMGs_subsampled.txt", sep = "\n")

table(Metagenome$library_layout)
summary(Metagenome$base_count)
View(Metagenome[is.na(Metagenome$base_count),])
sum(Metagenome$fastq_ftp == "")
sum(is.na(Metagenome$fastq_ftp))


subset_assemblies <- Metagenome %>% 
  filter(
    !is.na(lat),
    !is.na(lon),
    Metagenome$fastq_ftp != "",
    scientific_name %in% listMGs_subsampled,
    Metagenome$base_count >= 3000000000,
    instrument_model != "unspecified"
  ) 

summary(subset_assemblies$read_count)
table(subset_assemblies$library_layout)
summary(subset_assemblies$base_count)
table(subset_assemblies$instrument_model)
View(subset_assemblies[subset_assemblies$instrument_model == "unspecified", c("study_accession", "sample_description", "sample_title", "study_title", "base_count", "collection_date")])
# remove unspecified

sort(table(subset_assemblies$study_accession), decreasing = T)[1:100]
View(subset_assemblies[subset_assemblies$scientific_name == "biofilm metagenome", c("study_accession", "sample_description", "sample_title", "study_title")])

# check file size
hsize(sum(as.numeric(unlist(strsplit(subset_assemblies$fastq_bytes, ";")))))

# check if there are samples with more than one experiment/run
which(table(subset_assemblies$sample_accession) > 1)
View(subset_assemblies[subset_assemblies$sample_accession %in% names(which(table(subset_assemblies$sample_accession) > 1)), ])



# do manual refinement in Excel, then reload file into R

Metagenome2 <- read.table(
  "subset_MG_scientific_name_8391.txt", 
  sep="\t",
  header=T, 
  comment.char = "",
  quote = "\"",
  check.names = FALSE
) %>% 
  filter(study_accession != "PRJEB34634") # farm animal
dim(Metagenome2)

tmp <- table(Metagenome2$sample_accession) 
sum(tmp > 1)
sort(tmp, decreasing = T)[1:50]
View(Metagenome2[Metagenome2$sample_accession == "SAMEA7060539", c("last_updated", "study_title", "study_accession")])

Metagenome2$count <- Metagenome2$base_count/Metagenome2$read_count
sum(Metagenome$count < 150)
tmp <- Metagenome2[Metagenome2$count < 150, c("study_title")]

select_subset <- sample(nrow(Metagenome2), 50)
hsize(sum(as.numeric(unlist(strsplit(Metagenome2$fastq_bytes[select_subset], ";")))))

# define random subsets of ~50 runs each
select_random <- rep(1:150, 56)[shuffle(nrow(Metagenome2))]
write(select_random, "recover_select_random.txt")

for(i in unique(select_random)) {
  write.table(
    Metagenome2[select_random == i, ], 
    file=paste0("ENA_subsets_random/subset_", i, ".txt"),
    sep="\t", 
    row.names=FALSE,
    quote = F
  )
}


# show on map
maps::map(
  "world", # world coastlines
  resolution = 2, # highest resolution
  xlim = c(-180, 180), # area to be plotted
  ylim = c(-90, 90), 
  fill = T, # fill land masses
  col = "grey90", # land color
  border = NA,
  mar = c(0, 0, 0, 0)
)
points(
  Metagenome2$lon,
  Metagenome2$lat,
  pch = 16,
  cex = 0.2
)
