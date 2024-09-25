library(tidyverse)

# combine tables of coverm coverage search with min aligned percentage 50 and similarity 50
setwd("/storage/hdd1/chh/Mara_metaG_data_mining/v1/coverm/EX")
EX_mean_cov_sim50_aln50 <- list.files("mean_cov_sim50_aln50",pattern = "txt",full.names = T, recursive = TRUE)

EX_mean_cov_sim50_aln50_list <- map(
  EX_mean_cov_sim50_aln50,
  function(x) {
    tmp <- read.table(
      x,
      sep="\t",
      header=T,
      check.names = FALSE,
      row.names = 1
    ) %>% t() %>% as.data.frame()
    if(ncol(tmp) > 0) {
      return(tmp)
    }
  }
)
EX_mean_cov_sim50_aln50_df <- bind_rows(EX_mean_cov_sim50_aln50_list[!sapply(EX_mean_cov_sim50_aln50_list, is.null)])
rownames(EX_mean_cov_sim50_aln50_df) <- gsub("_filt Mean", "", rownames(EX_mean_cov_sim50_aln50_df))

write.table(EX_mean_cov_sim50_aln50_df, 'EX_mean_cov_sim50_aln50.txt', row.names = T, col.names = T, sep = "\t", quote = F)


# combine coverage for all studies
setwd("C:/Users/admin/OneDrive/PhD/Data_mining")
bins <- read.table("CoverM/EX/results_coverm/EX_mean_cov_sim50_aln50.txt", header = T, row.names = 1)
meta <- read.table("subset_MG_scientific_name_8391.txt", header = T, sep="\t", comment.char = "",
                   quote = "\"",
                   check.names = FALSE,
                   row.names = 79)

bins$mean_cov_sum <- rowSums(bins[ , c(1:5)], na.rm=TRUE)

bins_mean_cov_ordered <- bins[order(-bins$mean_cov_sum),]
merged <- merge(bins_mean_cov_ordered, meta, by=0)

row.names(merged) <- merged$Row.names

# only select coverage values and study_accession column
merged_study <- merged[c(2:7,112)]

tmp3 <- merged_study %>% 
  group_by(study_accession) %>% summarise_each(list(sum))

tmp3_ordered <- tmp3[order(-tmp3$mean_cov_sum),]

write.table(tmp3_ordered, file="CoverM/EX/metaG_studies_coverage_50_50",row.names = F, col.names = T, sep = "\t", quote = F)

