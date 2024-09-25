## this script is meant to be sourced during the dada2 pipeline to provide functions to create nSeq files

# function to count unique sequences in dadaFR files
 # used in functions below
getN <- function(x) sum(getUniques(x))


# function to create nSeq for each lib using merging by rownames to make sure nothing is mixed up
 # output includes columns for read numbers after demultiplexing, primer clipping, filtering, denoising and merging
 # and all separate numbers for fr and rf reads used to sum them up

create_nSeq <- function(lib_no_input, nSeq_file_input) { # lib_no_input needs to have '_' in front, makes it possible to give also NULL if only 1 library is analysed
  nSeq <-
    list(
  # 1. Demux read no.
  (
    read.table(nSeq_file_input, h = T, stringsAsFactors = F) %>% 
      tidyr::separate(SID, into = c("SID_n"), sep = "_")
  ),
  # 2. clipped and filtered read no.
  (
    get(paste0("filt_FR.out", lib_no_input)) %>%                            # file with numbers for clipped and filtered reads of fr reads
      as.data.frame() %>%                                             # reformat from some kind of list to data.frame
      tibble::rownames_to_column() %>%                                        # new column with former rownames
      tidyr::separate(rowname, into = c("SID_n"), sep = "_") %>%             # create from columns new column with sample IDs without ending
      dplyr::rename("clipped.fr" = "reads.in", "filtered.fr" = "reads.out")  # make new unique colnames to be able to distinguish with rf file
  ),
  (
    get(paste0("filt_RF.out", lib_no_input)) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column() %>% 
      tidyr::separate(rowname, into = c("SID_n"), sep = "_") %>% 
      dplyr::rename("clipped.rf" = "reads.in", "filtered.rf" = "reads.out")
  ),
  # 3. denoised read no.
  (
    sapply(get(paste0("dadaFR_R1", lib_no_input)), getN) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column() %>%
      dplyr::rename("SID_n" = "rowname", "denoised_fwd_fr" = ".")
  ),
  (
    sapply(get(paste0("dadaRF_R1", lib_no_input)), getN) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column() %>%
      dplyr::rename("SID_n" = "rowname", "denoised_fwd_rf" = ".")
  ),
  (
    sapply(get(paste0("dadaFR_R2", lib_no_input)), getN) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column() %>%
      dplyr::rename("SID_n" = "rowname", "denoised_rev_fr" = ".")
  ),
  (
    sapply(get(paste0("dadaRF_R2", lib_no_input)), getN) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column() %>%
      dplyr::rename("SID_n" = "rowname", "denoised_rev_rf" = ".")
  ),
  # 4. merged
  (
    sapply(get(paste0("mergers_FR", lib_no_input)), getN) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column() %>%
      dplyr::rename("SID_n" = "rowname", "merged_fr" = ".")
  ),
  (
    sapply(get(paste0("mergers_RF", lib_no_input)), getN) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column() %>%
      dplyr::rename("SID_n" = "rowname", "merged_rf" = ".")
  )) %>% 
  # 5. join all tables
  purrr::reduce(full_join, by = "SID_n") %>% 
  # 6. sum fr and rf reads 
  dplyr::mutate(Clipped = clipped.fr + clipped.rf,                   # sum fr and rf read numbers for clipped
         Filtered = filtered.fr + filtered.rf,                # ...for filtered            
         Denoised_fwd = denoised_fwd_fr + denoised_fwd_rf,    # ...for denoised fwd (=R1)
         Denoised_rev = denoised_rev_fr + denoised_rev_rf,    # ...for denoised rev (=R2)
         Merged = merged_fr + merged_rf)                      # ...for merged
  return(nSeq)
}


# this function merges the nSeq files created above
 # and adds read numbers after Chimera removal, removal of too long or short sequences 
 # and number of reads after classification and removal of unwanted taxa and with 
 # different bootstrap values
 # it uses a named list with all the different nSeq objects produced before as input

amend_nSeq <- function(lib_list_input) {
  # combine nSeq files from different libraries
  nSeq2 <- bind_rows(                       
    lib_list_input,
    .id = "lib"
    ) %>% 
    list(
      (
        # read no. after chimera removal
        rowSums(seqtab.nochim) %>% 
          as.data.frame() %>%
          tibble::rownames_to_column() %>%
          dplyr::rename("SID_n" = "rowname", "Nochim" = ".")
        ),
      (
        # read no. after discarding too long or too short sequences
        rowSums(seqtab.nochim2) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column() %>%
          dplyr::rename("SID_n" = "rowname", "Opt.length" = ".")
        ),
      (
        # read no. after classification and removing unwanted lineages
        rowSums(otu.filt) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column() %>%
          dplyr::rename("SID_n" = "rowname", "Classified_wanted" = ".")
      ),
      (
        # read no. after classification with final selected bootstrap value
        rowSums(otu.filt_f) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column() %>%
          dplyr::rename("SID_n" = "rowname", "Classified_final" = ".")
        ),
      (
        # read no. after classification with minimum bootstrap of 60
        rowSums(otu.filt_m) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column() %>%
          dplyr::rename("SID_n" = "rowname", "Classified_mod" = ".")
        )) %>% 
    # join tables
    purrr::reduce(full_join, by = "SID_n")
  return(nSeq2)
}


# this function calculates percent of reads lost in each step
 # it uses the nSeq table from above as input
perc_nSeq <- function(nSeq_input) {
  nSeq3 <- nSeq_input %>% 
    dplyr::select(lib, SID_n, Demux, Clipped, Filtered, Denoised_fwd, Denoised_rev,   # select only columns with combined reads
                  Merged, Nochim, Opt.length, Classified_wanted, Classified_final, Classified_mod) %>% 
    dplyr::mutate(across(where(is.numeric), ~ (. / Demux) * 100)) %>%                 # calculate percent of reads retained in each step compared to reads after demultiplexing
    dplyr::mutate(across(where(is.numeric), round, 2))                                # round to 2 digits after decimal
 return(nSeq3) 
}
