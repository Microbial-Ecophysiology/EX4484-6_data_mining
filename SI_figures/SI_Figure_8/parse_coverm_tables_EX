# Parse coverm results for coverage, relative abundance and breadth of coverage


require(tidyverse)

setwd("/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs/coverm_v2/EX_quantified/")

# parse for EX4484-6 MAGs only ####

EX_cov_95 <- list.files(pattern = "cov_95.txt",full.names = T, recursive = TRUE)
EX_cov_90 <- list.files(pattern = "cov_90.txt",full.names = T, recursive = TRUE)

EX_breadth_95 <- list.files(pattern = "breadth_95.txt",full.names = T, recursive = TRUE)
EX_breadth_90 <- list.files(pattern = "breadth_90.txt",full.names = T, recursive = TRUE)

EX_rel_95 <- list.files(pattern = "rel_95.txt",full.names = T, recursive = TRUE)
EX_rel_90 <- list.files(pattern = "rel_90.txt",full.names = T, recursive = TRUE)




# EX coverage ####
## 95% ####
EX_cov_95_list <- map(
  EX_cov_95,
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
EX_cov_95_df <- bind_rows(EX_cov_95_list[!sapply(EX_cov_95_list, is.null)])
rownames(EX_cov_95_df) <- gsub(" Mean", "", rownames(EX_cov_95_df))

write.table(EX_cov_95_df, "../EX_cov_95.txt", row.names = T, col.names = T, sep = "\t", quote = F)

## 90% ####
EX_cov_90_list <- map(
  EX_cov_90,
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
EX_cov_90_df <- bind_rows(EX_cov_90_list[!sapply(EX_cov_90_list, is.null)])
rownames(EX_cov_90_df) <- gsub(" Mean", "", rownames(EX_cov_90_df))

write.table(EX_cov_90_df, "../EX_cov_90.txt", row.names = T, col.names = T, sep = "\t", quote = F)



# EX breadth of coverage ####
## 95% ####
EX_breadth_95_list <- map(
  EX_breadth_95,
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
EX_breadth_95_df <- bind_rows(EX_breadth_95_list[!sapply(EX_breadth_95_list, is.null)])
rownames(EX_breadth_95_df) <- gsub(" Covered Fraction", "", rownames(EX_breadth_95_df))

write.table(EX_breadth_95_df, "../EX_breadth_95.txt", row.names = T, col.names = T, sep = "\t", quote = F)

## 90% ####
EX_breadth_90_list <- map(
  EX_breadth_90,
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
EX_breadth_90_df <- bind_rows(EX_breadth_90_list[!sapply(EX_breadth_90_list, is.null)])
rownames(EX_breadth_90_df) <- gsub(" Covered Fraction", "", rownames(EX_breadth_90_df))

write.table(EX_breadth_90_df, "../EX_breadth_90.txt", row.names = T, col.names = T, sep = "\t", quote = F)


# EX relative abundance ####
## 95% ####
EX_rel_95_list <- map(
  EX_rel_95,
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
EX_rel_95_df <- bind_rows(EX_rel_95_list[!sapply(EX_rel_95_list, is.null)])
rownames(EX_rel_95_df) <- gsub(" Relative Abundance (%)", "", rownames(EX_rel_95_df))

write.table(EX_rel_95_df, "../EX_rel_95.txt", row.names = T, col.names = T, sep = "\t", quote = F)

## 90% ####
EX_rel_90_list <- map(
  EX_rel_90,
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
EX_rel_90_df <- bind_rows(EX_rel_90_list[!sapply(EX_rel_90_list, is.null)])
rownames(EX_rel_90_df) <- gsub(" Relative Abundance (%)", "", rownames(EX_rel_90_df))

write.table(EX_rel_90_df, "../EX_rel_90.txt", row.names = T, col.names = T, sep = "\t", quote = F)




# Parse coverm results of all Thermoplasmatota for coverage, relative abundance and breadth of coverage ####


require(tidyverse)

setwd("/storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs/coverm_v2/Thermo_quantified/")


EX_cov_95 <- list.files(pattern = "cov_95.txt",full.names = T, recursive = TRUE)
EX_cov_90 <- list.files(pattern = "cov_90.txt",full.names = T, recursive = TRUE)

EX_breadth_95 <- list.files(pattern = "breadth_95.txt",full.names = T, recursive = TRUE)
EX_breadth_90 <- list.files(pattern = "breadth_90.txt",full.names = T, recursive = TRUE)

EX_rel_95 <- list.files(pattern = "rel_95.txt",full.names = T, recursive = TRUE)
EX_rel_90 <- list.files(pattern = "rel_90.txt",full.names = T, recursive = TRUE)




# EX coverage ####
## 95% ####
EX_cov_95_list <- map(
  EX_cov_95,
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
EX_cov_95_df <- bind_rows(EX_cov_95_list[!sapply(EX_cov_95_list, is.null)])
rownames(EX_cov_95_df) <- gsub(" Mean", "", rownames(EX_cov_95_df))

write.table(EX_cov_95_df, "../Thermo_cov_95.txt", row.names = T, col.names = T, sep = "\t", quote = F)

## 90% ####
EX_cov_90_list <- map(
  EX_cov_90,
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
EX_cov_90_df <- bind_rows(EX_cov_90_list[!sapply(EX_cov_90_list, is.null)])
rownames(EX_cov_90_df) <- gsub(" Mean", "", rownames(EX_cov_90_df))

write.table(EX_cov_90_df, "../Thermo_cov_90.txt", row.names = T, col.names = T, sep = "\t", quote = F)



# EX breadth of coverage ####
## 95% ####
EX_breadth_95_list <- map(
  EX_breadth_95,
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
EX_breadth_95_df <- bind_rows(EX_breadth_95_list[!sapply(EX_breadth_95_list, is.null)])
rownames(EX_breadth_95_df) <- gsub(" Covered Fraction", "", rownames(EX_breadth_95_df))

write.table(EX_breadth_95_df, "../Thermo_breadth_95.txt", row.names = T, col.names = T, sep = "\t", quote = F)

## 90% ####
EX_breadth_90_list <- map(
  EX_breadth_90,
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
EX_breadth_90_df <- bind_rows(EX_breadth_90_list[!sapply(EX_breadth_90_list, is.null)])
rownames(EX_breadth_90_df) <- gsub(" Covered Fraction", "", rownames(EX_breadth_90_df))

write.table(EX_breadth_90_df, "../Thermo_breadth_90.txt", row.names = T, col.names = T, sep = "\t", quote = F)


# EX relative abundance ####
## 95% ####
EX_rel_95_list <- map(
  EX_rel_95,
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
EX_rel_95_df <- bind_rows(EX_rel_95_list[!sapply(EX_rel_95_list, is.null)])
rownames(EX_rel_95_df) <- gsub(" Relative Abundance (%)", "", rownames(EX_rel_95_df))

write.table(EX_rel_95_df, "../Thermo_rel_95.txt", row.names = T, col.names = T, sep = "\t", quote = F)

## 90% ####
EX_rel_90_list <- map(
  EX_rel_90,
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
EX_rel_90_df <- bind_rows(EX_rel_90_list[!sapply(EX_rel_90_list, is.null)])
rownames(EX_rel_90_df) <- gsub(" Relative Abundance (%)", "", rownames(EX_rel_90_df))

write.table(EX_rel_90_df, "../Thermo_rel_90.txt", row.names = T, col.names = T, sep = "\t", quote = F)
