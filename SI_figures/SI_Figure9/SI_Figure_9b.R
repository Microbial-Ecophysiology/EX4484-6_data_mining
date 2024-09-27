
# parse annotation

library(tidyverse)
library(reshape)
library(patchwork)
library(RColorBrewer)
library(vegan)
library(ImpactEffectsize)

setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Data_mining/Annotation_v3/hypothetical_allThermo")


load("hypothetical_unknown_OG_analysis.RData")

# load annotation ####

## kegg ----
kegg <- read.csv("concat_kegg_annotation_best.txt", 
                 sep="\t", 
                 header=TRUE, 
                 fill=TRUE, 
                 quote="",
                 col.names=c("qseqid", "kegg_sseqid","kegg_pident","kegg_length","kegg_mismatch","kegg_gapopen","kegg_qstart",
                             "kegg_qend","kegg_sstart","kegg_send","kegg_evalue","kegg_bitscore","kegg_bsr","kegg_KO","kegg_KO_name","kegg_taxon","kegg_gene"))
kegg <- subset(kegg, kegg_sseqid != "sseqid")

kegg_EX <- read.csv("../full/concat_kegg_annotation_best.txt", 
                    sep="\t", 
                    header=TRUE, 
                    fill=TRUE, 
                    quote="",
                    col.names=c("qseqid", "kegg_sseqid","kegg_pident","kegg_length","kegg_mismatch","kegg_gapopen","kegg_qstart",
                                "kegg_qend","kegg_sstart","kegg_send","kegg_evalue","kegg_bitscore","kegg_bsr","kegg_KO","kegg_KO_name","kegg_taxon","kegg_gene"))
kegg_EX <- subset(kegg_EX, kegg_sseqid != "sseqid")
# kegg_EX$qseqid <- gsub(".genes", "", kegg_EX$qseqid)
# kegg_EX$qseqid <- gsub("refined-contigs", "refined.contigs", kegg_EX$qseqid)


## NR ----
NR <- read.csv("concat_NR_annotation_best.txt", 
               sep="\t", 
               header=TRUE, 
               fill=TRUE, 
               quote="",
               col.names=c("qseqid", "NR_sseqid","NR_pident","NR_length","NR_mismatch","NR_gapopen","NR_qstart",
                           "NR_qend","NR_sstart","NR_send","NR_evalue","NR_bitscore","NR_qlen","NR_slen","NR_stitle","NR_staxids","NR_sscinames","NR_salltitles","NR_bsr"))
NR <- subset(NR, NR_sseqid != "sseqid")

NR_EX <- read.csv("../full/concat_NR_annotation_best.txt", 
                  sep="\t", 
                  header=TRUE, 
                  fill=TRUE, 
                  quote="",
                  col.names=c("qseqid", "NR_sseqid","NR_pident","NR_length","NR_mismatch","NR_gapopen","NR_qstart",
                              "NR_qend","NR_sstart","NR_send","NR_evalue","NR_bitscore","NR_qlen","NR_slen","NR_stitle","NR_staxids","NR_sscinames","NR_salltitles","NR_bsr"))
NR_EX <- subset(NR_EX, NR_sseqid != "sseqid")
# NR_EX$qseqid <- gsub(".genes", "", NR_EX$qseqid)
# NR_EX$qseqid <- gsub("refined-contigs", "refined.contigs", NR_EX$qseqid)


# Combine annotation data of kegg and NR
annotation <- full_join(kegg, NR, by = c("qseqid")) 
annotation_EX <- full_join(kegg_EX, NR_EX, by = c("qseqid")) 
all.equal(colnames(annotation), colnames(annotation_EX))
annotation_all <- rbind(annotation, annotation_EX)

# remove qseqid's of outgroup taxa
# Korarchaeia, Methanomethylicia, Methanosarcinia, Nitrososphaeria, Bathy and 2 EX MAGs, which are duplicated in the data set

# GCA_001918745.1_ASM191874v1_genomic, GCA_001593935.1_ASM159393v1_genomic, GCA_018396415.1_ASM1839641v1_genomic, GCA_003601775.1_ASM360177v1_genomic
# GCF_000019605.1_ASM1960v1_genomic, GCA_001717035.1_ASM171703v1_genomic, GCA_000200715.1_ASM20071v1_genomic,
# GCF_000007345.1_ASM734v1_genomic, GCF_000970085.1_ASM97008v1_genomic, GCF_000970205.1_ASM97020v1_genomic,GCF_000970265.1_ASM97026v1_genomic,
# GCF_000970285.1_ASM97028v1_genomic, GCF_001304615.2_ASM130461v2_genomic, GCF_002287235.1_ASM228723v1_genomic, 
# PRJNA713414_bin_118_orig_1_refined-contigs,PRJNA713414_bin_178_orig_1_refined-contigs

accession_numbers_to_exclude <- c(
  "GCA_001918745.1_ASM191874v1_genomic", "GCA_001593935.1_ASM159393v1_genomic", 
  "GCA_018396415.1_ASM1839641v1_genomic", "GCA_003601775.1_ASM360177v1_genomic",
  "GCF_000019605.1_ASM1960v1_genomic", "GCA_001717035.1_ASM171703v1_genomic", 
  "GCA_000200715.1_ASM20071v1_genomic", "GCF_000007345.1_ASM734v1_genomic", 
  "GCF_000970085.1_ASM97008v1_genomic", "GCF_000970205.1_ASM97020v1_genomic", 
  "GCF_000970265.1_ASM97026v1_genomic", "GCF_000970285.1_ASM97028v1_genomic", 
  "GCF_001304615.2_ASM130461v2_genomic", "GCF_002287235.1_ASM228723v1_genomic",
  "PRJNA713414_bin_118_orig_1_refined-contigs","PRJNA713414_bin_178_orig_1_refined-contigs"
)

# Create a pattern that matches any qseqid starting with any of the accession numbers
pattern_exclude <- paste0("^(", paste(accession_numbers_to_exclude, collapse = "|"), ")")

# Filter out rows with qseqid matching the pattern
annotation_filtered <- annotation_all %>%
  filter(!grepl(pattern_exclude, qseqid))

annotation_thermo <- annotation
annotation <- annotation_all

annotation$kegg_KO[is.na(annotation$kegg_KO)] <- ""
annotation$kegg_gene[is.na(annotation$kegg_gene)] <- ""
annotation$NR_stitle[is.na(annotation$NR_stitle)] <- ""

uniq_genes <- unique(annotation$kegg_gene[!is.na(annotation$kegg_KO)])
uniq_stitle <- unique(annotation$NR_stitle)
grep("hypothetical|[pP]rotein of unknown function", uniq_genes, value = T)
head(grep("hypothetical|[pP]rotein of unknown function", uniq_stitle, value = T))


tmp <- data.frame(
  tmp_KO = ifelse(annotation$kegg_KO != "", "annotated", ""),
  tmp_gene = ifelse(!grepl("hypothetical|[pP]rotein of unknown function", annotation$kegg_gene) & annotation$kegg_gene != "", "annotated", ""),
  tmp_stitle = ifelse(!grepl("hypothetical|[pP]rotein of unknown function", annotation$NR_stitle) & annotation$NR_stitle != "", "annotated", "")
)
annotation$annotation_status <- ifelse(apply(tmp, 1, function(x) sum(x == "") == 3), "hypothetical", "annotated")
tmp <- annotation[, c("qseqid", "annotation_status")]

# load Orthogroups.txt
orthogroups <- scan("../../Orthofinder_allThermo/Orthogroups.txt", what = "character", sep = "\n")
length(orthogroups)
tmp_df <- data.frame(
  do.call(
    "rbind",
    strsplit(orthogroups, ": ")
  )
)
tmp_list <- strsplit(tmp_df$X2, " ", fixed = T)
tmp_df_filt <- tmp_df[sapply(tmp_list, length) > 1, ]
tmp_list_filt <- tmp_list[sapply(tmp_list, length) > 1]
orthogroups_map <- data.frame(
  OG = rep(tmp_df_filt$X1, sapply(tmp_list_filt, length)),
  qseqid = unlist(tmp_list_filt)
)
table(table(orthogroups_map$OG))
orthogroups_annotation_status <- left_join(orthogroups_map, tmp, by = c("qseqid")) 
orthogroups_annotation_status$annotation_status[is.na(orthogroups_annotation_status$annotation_status)] <- "unknown"
orthogroups_annotation_status$annotation_status <- as.factor(orthogroups_annotation_status$annotation_status)
table(orthogroups_annotation_status$annotation_status)
OG_annotation_status <- data.frame(do.call("rbind", by(orthogroups_annotation_status$annotation_status, orthogroups_annotation_status$OG, table)))
dim(OG_annotation_status)
OG_annotation_status$OG_status <- ifelse(OG_annotation_status$annotated > 0, "annotated", ifelse(OG_annotation_status$hypothetical > 0, "hypothetical", "unknown"))
table(OG_annotation_status$OG_status)

OG_counts <- read.table("../../Orthofinder_allThermo/Orthogroups.GeneCount.tsv", h = T, sep = "\t", check.names = F)
OG_counts$OG_status <- OG_annotation_status[OG_counts$Orthogroup, "OG_status"]

# remove outgroup MAGs as above 
OG_counts_filtered <- OG_counts %>%
  select(-all_of(accession_numbers_to_exclude)) # excluding outgroup MAGs with accession_numbers_to_exclude

rm(OG_counts)

# filter rows only present in EX genomes
# vegan decostand pa

OG_counts_pa <- decostand(OG_counts_filtered[, 2:419], "pa")

# rowSums will give you number of genomes present
OG_counts_filtered$n_genomes_shared <- rowSums(OG_counts_pa)

# filter rows only present in EX genomes
EX_MAGs_unique_tmp <- annotation_EX %>% 
  separate(col = qseqid, into = c("genome","cluster"), sep = "___")
EX_MAGs_unique <- unique(EX_MAGs_unique_tmp$genome)
# remove PRJNA713414 from list
remove <- c("PRJNA713414_bin_118_orig_1_refined-contigs","PRJNA713414_bin_178_orig_1_refined-contigs")
EX_MAGs_unique <- setdiff(EX_MAGs_unique, remove)
EX_MAGs_unique <- gsub(".genes$", "", EX_MAGs_unique)

# Select the OG column and unique EX MAGs columns
subset_OG_counts_EX <- OG_counts_filtered %>%
  select(Orthogroup, OG_status, all_of(EX_MAGs_unique))

# Create a logical vector indicating rows with counts in EX_MAGs_unique
has_counts_in_EX_MAGs <- rowSums(subset_OG_counts_EX[, c(3:37)]) > 0 # 6231 OG in total

## Orthogroups of EX ####
subset_og_counts_EX_only <- subset_OG_counts_EX[has_counts_in_EX_MAGs, ]

subset_og_counts_EX_only_df <- melt(subset_og_counts_EX_only)

## load family information and join with orthogroup gene count table ####
families <- read.table("new_family_groups.txt", header = T, sep="\t", comment.char = "",
                       quote = "\"",
                       check.names = FALSE)


# combine both tables
og_counts_EX_only_merged_df <- full_join(subset_og_counts_EX_only_df, families, by=c("variable"="MAG"))


# find orthogroups, which have counts in sediment and water ####
summary_df <- og_counts_EX_only_merged_df %>%
  group_by(Orthogroup) %>%
  summarize(
    water_present = any(environment == "water" & value > 0),
    sediment_present = any(environment == "sediment" & value > 0)
  )

## Orthogroups only in water ####
only_water <- summary_df %>%
  filter(water_present == TRUE & sediment_present == FALSE) %>%
  pull(Orthogroup) # 2116

## Orthogroups only in sediment ####
only_sediment <- summary_df %>%
  filter(water_present == FALSE & sediment_present == TRUE) %>%
  pull(Orthogroup) # 2308

## Orthogroups in both environments ####
both_environments <- summary_df %>%
  filter(water_present == TRUE & sediment_present == TRUE) %>%
  pull(Orthogroup) # 1807


## combine info in table ####
og_counts_EX_only_merged_df_2 <- og_counts_EX_only_merged_df %>% 
  mutate(env_group = case_when(Orthogroup %in% only_water ~ "water",
                               Orthogroup %in% only_sediment ~ "sediment",
                               Orthogroup %in% both_environments ~ "both"))



# OGs per different environment_detailed ####
summary_env_detailes_df <- og_counts_EX_only_merged_df %>%
  group_by(Orthogroup) %>%
  summarize(
    water_present = any(environment_detailed == "marine water" & value > 0),
    sediment_present = any(environment_detailed == "marine sediment" & value > 0),
    lake_present = any(environment_detailed == "lake sediment" & value > 0),
    hydrothermal_present = any(environment_detailed == "marine hydrothermal sediment" & value > 0),
    cold_seep_present = any(environment_detailed == "cold seep fluid" & value > 0)
  )

## Orthogroups in water ####
only_water_env_detailed <- summary_env_detailes_df %>%
  filter(water_present == T & sediment_present == F, lake_present == F, hydrothermal_present == F, cold_seep_present ==F) %>%
  pull(Orthogroup) # 2116

## Orthogroups in sediment ####
only_sediment_env_detailed <- summary_env_detailes_df %>%
  filter(water_present == F & sediment_present == T, lake_present == F, hydrothermal_present == F, cold_seep_present ==F) %>%
  pull(Orthogroup) # 911

## Orthogroups in lake ####
only_lake_env_detailed <- summary_env_detailes_df %>%
  filter(water_present == F & sediment_present == F, lake_present == T, hydrothermal_present == F, cold_seep_present ==F) %>%
  pull(Orthogroup) # 168

## Orthogroups in hydrothermal sediment ####
only_hydrothermal_env_detailed <- summary_env_detailes_df %>%
  filter(water_present == F & sediment_present == F, lake_present == F, hydrothermal_present == T, cold_seep_present ==F) %>%
  pull(Orthogroup) # 431

## Orthogroups in cold seep ####
only_coldseep_env_detailed <- summary_env_detailes_df %>%
  filter(water_present == F & sediment_present == F, lake_present == F, hydrothermal_present == F, cold_seep_present ==T) %>%
  pull(Orthogroup) # 23

all_env_detailed <- summary_env_detailes_df %>%
  filter(water_present == T & sediment_present == T, lake_present == T, hydrothermal_present == T, cold_seep_present ==T) %>%
  pull(Orthogroup) # 944



# per Family ####

## Summarize counts by Orthogroup and environment ####
summary_family_df <- og_counts_EX_only_merged_df %>%
  group_by(Orthogroup) %>%
  summarize(
    Family_1A = any(Family == "1A" & value > 0),
    Family_1B = any(Family == "1B" & value > 0),
    Family_2 = any(Family == "2" & value > 0),
    Family_3A = any(Family == "3A" & value > 0),
    Family_3B = any(Family == "3B" & value > 0),
    Family_4 = any(Family == "4" & value > 0)
  )

## Orthogroups only in single families ####
only_Family_1A <- summary_family_df %>%
  filter(Family_1A == TRUE & Family_1B == FALSE & Family_2 == FALSE & Family_3A == FALSE & Family_3B == FALSE & Family_4 == FALSE) %>%
  pull(Orthogroup) # 1179
Family_1A_total <- summary_family_df %>%
  filter(Family_1A == TRUE)%>%
  pull(Orthogroup) #2961
# 39.8 % of orthogroups within family 1A are family-specific

only_Family_1B <- summary_family_df %>%
  filter(Family_1A == FALSE & Family_1B == TRUE & Family_2 == FALSE & Family_3A == FALSE & Family_3B == FALSE & Family_4 == FALSE) %>%
  pull(Orthogroup) # 377
Family_1B_total <- summary_family_df %>%
  filter(Family_1B == TRUE)%>%
  pull(Orthogroup) #1596
# 23.6 % of orthogroups within family 1B are family-specific

only_Family_2 <- summary_family_df %>%
  filter(Family_1A == FALSE & Family_1B == FALSE & Family_2 == TRUE & Family_3A == FALSE & Family_3B == FALSE & Family_4 == FALSE) %>%
  pull(Orthogroup) # 427
Family_2_total <- summary_family_df %>%
  filter(Family_2 == TRUE)%>%
  pull(Orthogroup) #2011
# 21.2 % of orthogroups within family 2 are family-specific

only_Family_3A <- summary_family_df %>%
  filter(Family_1A == FALSE & Family_1B == FALSE & Family_2 == FALSE & Family_3A == TRUE & Family_3B == FALSE & Family_4 == FALSE) %>%
  pull(Orthogroup) # 848
Family_3A_total <- summary_family_df %>%
  filter(Family_3A == TRUE)%>%
  pull(Orthogroup) #2680
# 31.6 % of orthogroups within family 3A are family-specific

only_Family_3B <- summary_family_df %>%
  filter(Family_1A == FALSE & Family_1B == FALSE & Family_2 == FALSE & Family_3A == FALSE & Family_3B == TRUE & Family_4 == FALSE) %>%
  pull(Orthogroup) # 144
Family_3B_total <- summary_family_df %>%
  filter(Family_3B == TRUE)%>%
  pull(Orthogroup) #1543
# 9.3 % of orthogroups within family 3B are family-specific

only_Family_4 <- summary_family_df %>%
  filter(Family_1A == FALSE & Family_1B == FALSE & Family_2 == FALSE & Family_3A == FALSE & Family_3B == FALSE & Family_4 == TRUE) %>%
  pull(Orthogroup) # 800
Family_4_total <- summary_family_df %>%
  filter(Family_4 == TRUE)%>%
  pull(Orthogroup) #2027
# 39.5 % of orthogroups within family 4 are family-specific

# Orthogroups in all families
all_families <- summary_family_df %>%
  filter(Family_1A == TRUE & Family_1B == TRUE & Family_2 == TRUE & Family_3A == TRUE & Family_3B == TRUE & Family_4 == TRUE) %>%
  pull(Orthogroup) # 586


## add categories to table ####
## combine info in table 
og_counts_EX_only_merged_df_3 <- og_counts_EX_only_merged_df_2 %>% 
  mutate(family_group = case_when(Orthogroup %in% only_Family_1A ~ "1A",
                                  Orthogroup %in% only_Family_1B ~ "1B",
                                  Orthogroup %in% only_Family_2 ~ "2",
                                  Orthogroup %in% only_Family_3A ~ "3A",
                                  Orthogroup %in% only_Family_3B ~ "3B",
                                  Orthogroup %in% only_Family_4 ~ "4",
                                  Orthogroup %in% all_families ~ "all_families"))
og_counts_EX_only_merged_df_3$family_group <- og_counts_EX_only_merged_df_3$family_group %>% replace_na("multiple_families")

# 4361 OG are in categories, remaining 1870 OG are in multiple families but not all


# Barplot for environment ####
# create a table with OG counts for each environmental group and OG_status
unique_EX_OG_merged_df <- og_counts_EX_only_merged_df_3 %>%
  distinct(Orthogroup, OG_status, env_group)

count_EX_OG_merged_df <- unique_EX_OG_merged_df %>%
  group_by(env_group, OG_status) %>%
  summarise(count = n(), .groups = 'drop')


# Create the stacked barplot
ggplot(count_EX_OG_merged_df, aes(x = env_group, y = count, fill = OG_status)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Orthogroups Count by Environment Group and OG Status",
       x = "Environment Group",
       y = "Count of Orthogroups",
       fill = "OG Status") +
  theme_minimal()


# Barplot for family ####
# create a table with OG counts for each family group and OG_status

unique_EX_OG_merged_family_df <- og_counts_EX_only_merged_df_3 %>%
  distinct(Orthogroup, OG_status, family_group)

count_EX_OG_merged_family_df <- unique_EX_OG_merged_family_df %>%
  group_by(family_group, OG_status) %>%
  summarise(count = n(), .groups = 'drop')

# Create the stacked barplot
ggplot(count_EX_OG_merged_family_df, aes(x = family_group, y = count, fill = OG_status)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Orthogroups Count by Family Group and OG Status",
       x = "Family Group",
       y = "Count of Orthogroups",
       fill = "OG Status") +
  theme_minimal()





# Upset plot
head(og_counts_EX_only_merged_df_3)
require(reshape)
OG_fam_pa <- cast(og_counts_EX_only_merged_df_3, "Orthogroup ~ Family", value = "value", fun.aggregate = sum)
rownames(OG_fam_pa) <- OG_fam_pa$Orthogroup
OG_fam_pa <- OG_fam_pa[, -1]
OG_fam_pa[OG_fam_pa > 0] <- 1
colSums(OG_fam_pa)
OG_fam_pa_list <- lapply(as.list(OG_fam_pa), function(x) rownames(OG_fam_pa)[ x == 1])


require(UpSetR)

# nintersects = 30, as these OGs represent 95% of all OGs in data (5919,45)

png(file="upset_plot_OG.png",
    width=18, height=10, units = "cm", res = 300)
upset(
  data = upset_data,
  nsets = 6,
  nintersects = 30,
  #sets = colnames(OG_fam_pa),
  sets = c("4","3B","3A","2","1B","1A"),
  keep.order = T,
  order.by = "freq",
  mainbar.y.label = "Orthogroup Intersections",
  sets.x.label = "Orthogroups per genome",
  text.scale = c(0.8, 0.7,0.8, 0.7, 0.8, 0.7),
  mb.ratio = c(0.65, 0.35),
  query.legend = "bottom",
  queries = list(
    list(query = intersects, params = list("3A"), 
         color = "#35978f", active = T, query.name = "water specific Orthogroups"),
    list(query = intersects, params = list("3B"), 
         color = "#35978f", active = T, query.name = ""),
    list(query = intersects, params = list("4"), 
         color = "#35978f", active = T, query.name = ""),
    list(query = intersects, params = list("3A","3B"), 
         color = "#35978f", active = T, query.name = ""),
    list(query = intersects, params = list("3A","4"), 
         color = "#35978f", active = T, query.name = ""),
    # list(query = intersects, params = list("3B","4"), 
    #      color = "#2b8cbe", active = T, query.name = ""),
    list(query = intersects, params = list("3A","3B","4"), 
         color = "#35978f", active = T, query.name = ""),
    list(query = intersects, params = list("1A"),
         color = "#bf812d", active = T, query.name = "sediment specific Orthogroups"),
    list(query = intersects, params = list("1B"),
         color = "#bf812d", active = T, query.name = ""),
    list(query = intersects, params = list("2"),
         color = "#bf812d", active = T, query.name = ""),
    list(query = intersects, params = list("1A","1B"),
         color = "#bf812d", active = T, query.name = ""),
    list(query = intersects, params = list("1A","2"),
         color = "#bf812d", active = T, query.name = ""),
    # list(query = intersects, params = list("1B","2"),
    #      color = "#74c476", active = T, query.name = ""),
    list(query = intersects, params = list("1A","1B","2"),
         color = "#bf812d", active = T, query.name = "")))
dev.off()


pdf(file="upset_plot_OG.pdf",
    width=7, height=4)
upset(
  data = upset_data,
  nsets = 6,
  nintersects = 30,
  #sets = colnames(OG_fam_pa),
  sets = c("4","3B","3A","2","1B","1A"),
  keep.order = T,
  order.by = "freq",
  mainbar.y.label = "Orthogroup Intersections",
  sets.x.label = "Orthogroups per genome",
  text.scale = c(0.8, 0.7,0.8, 0.7, 0.8, 0.7),
  mb.ratio = c(0.65, 0.35),
  query.legend = "bottom",
  queries = list(
    list(query = intersects, params = list("3A"), 
         color = "#35978f", active = T, query.name = "water specific Orthogroups"),
    list(query = intersects, params = list("3B"), 
         color = "#35978f", active = T, query.name = ""),
    list(query = intersects, params = list("4"), 
         color = "#35978f", active = T, query.name = ""),
    list(query = intersects, params = list("3A","3B"), 
         color = "#35978f", active = T, query.name = ""),
    list(query = intersects, params = list("3A","4"), 
         color = "#35978f", active = T, query.name = ""),
    # list(query = intersects, params = list("3B","4"), 
    #      color = "#2b8cbe", active = T, query.name = ""),
    list(query = intersects, params = list("3A","3B","4"), 
         color = "#35978f", active = T, query.name = ""),
    list(query = intersects, params = list("1A"),
         color = "#bf812d", active = T, query.name = "sediment specific Orthogroups"),
    list(query = intersects, params = list("1B"),
         color = "#bf812d", active = T, query.name = ""),
    list(query = intersects, params = list("2"),
         color = "#bf812d", active = T, query.name = ""),
    list(query = intersects, params = list("1A","1B"),
         color = "#bf812d", active = T, query.name = ""),
    list(query = intersects, params = list("1A","2"),
         color = "#bf812d", active = T, query.name = ""),
    # list(query = intersects, params = list("1B","2"),
    #      color = "#74c476", active = T, query.name = ""),
    list(query = intersects, params = list("1A","1B","2"),
         color = "#bf812d", active = T, query.name = "")))
dev.off()
