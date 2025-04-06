#Script for overall statistics based on meta data from screened metagenomic runs

library(tidyverse)
library(reshape2)
require(R.utils)
library(patchwork)
library(ggplot2)

setwd("/Users/mara/Library/CloudStorage/OneDrive-Personal/PhD/Data_mining/")
load("Overview_stats/Figure1_overview_stats.RData")

# studies from OMD
OMD <- read.table("Overview_stats/genome_export_OMD.txt", header = T, sep="\t", comment.char = "",
                  quote = "\"",
                  check.names = FALSE)

#subset table and only choose biosamples with ENA sample accession
OMD_filtered <- OMD %>%
  filter(str_starts(biosample, "SAM")) %>%
  mutate(biosample = str_replace(biosample,"-.*", ""))

unique(OMD_filtered$biosample)

# Get unique biosample values and create a comma-separated list
biosample_list <- OMD_filtered %>%
  distinct(biosample) %>%
  pull(biosample) 

split_lists <- split(biosample_list, cut(seq_along(biosample_list), 20, labels = FALSE))

for (i in seq_along(split_lists)) {
  writeLines(paste(split_lists[[i]], collapse = ","), paste0("biosample_list_part_", i, ".txt"))
}

# Set the directory path
file_list <- list.files(path = "Overview_stats/Biosample_OMD_meta", pattern = "\\.txt$", full.names = TRUE)

# Read and combine tables
OMD_meta <- do.call(rbind, 
                    lapply(
                      file_list, 
                      read.table, 
                      header = TRUE, 
                      fill = TRUE, 
                      quote = "",
                      sep = "\t"))

# 40 samples are missing. Likely the ones with additional strings
missing_biosamples <- OMD_filtered %>%
  filter(!biosample %in% OMD_meta$sample_accession)  # Check which biosample is missing

unique(missing_biosamples$biosample)
# [1] "SAMEA115094352" "SAMEA115094353" "SAMEA115094354" "SAMEA115094355" "SAMEA115094356" "SAMEA115094357" "SAMEA115094358" "SAMEA115094359" "SAMEA115094360" "SAMEA115094361" "SAMEA115094362"
# [12] "SAMEA115094363" "SAMEA115094364" "SAMEA115094365" "SAMEA115094366" "SAMEA115094367" "SAMEA115094368" "SAMEA115094369" "SAMEA115094370" "SAMEA115094371" "SAMEA115094372" "SAMEA115094373"
# [23] "SAMEA115094374" "SAMEA115094375" "SAMEA115094376" "SAMEA115094377" "SAMEA115094378" "SAMEA115094379" "SAMEA115094380" "SAMEA115094381" "SAMEA115094382" "SAMEA115094383" "SAMEA115094384"
# [34] "SAMEA115094385" "SAMEA115094386" "SAMEA115094387" "SAMEA115094388" "SAMEA115094389" "SAMEA115094390" "SAMEA115094391"

# All unique missing biosamples did not have longitude and latitude data, therefore can be ignored.

# subset OMD to biosample, longitude and latitude
OMD_short <- unique(OMD_filtered[c("biosample", "latitude","longitude")])

OMD_combined <- left_join(OMD_meta,
                          OMD_short,
                          by = c("sample_accession"="biosample"))

OMD_combined <- OMD_combined %>%
  mutate(DB = ifelse(sample_accession %in% c("SAMN11854494", "SAMN10231893", "SAMN10231913", "SAMN10231914", 
                                      "SAMN10231894", "SAMN10231904", "SAMEA2620113", "SAMN22703512",
                                      "SAMN22703513","SAMN22703514","SAMN30647027"), 
                     "EX", "OMD"))

# subset data frame
OMD_combined <- OMD_combined[c("sample_accession","scientific_name","latitude","longitude","DB")]
colnames(OMD_combined) <- c("sample_accession","scientific_name","lat","lon","DB")


# own screened data ####
metaG_mining <- read.table(
  "subset_MG_scientific_name_8287_final.txt", 
  sep="\t",
  header=T, 
  comment.char = "",
  quote = "\"",
  check.names = FALSE
)

metaG_mining_short <- unique(metaG_mining[c("sample_accession", "study_accession","lat","lon", "scientific_name")])

metaG_mining_short <- metaG_mining_short %>%
  mutate(DB = ifelse(study_accession %in% c("PRJNA362212", "PRJNA480137", 
                                             "PRJNA690107", "PRJNA707313", "PRJNA713414", 
                                            "PRJNA368391",
                                             "PRJNA531756","PRJNA541421","PRJNA704804",
                                            "PRJNA721298","PRJNA889212"), 
                     "EX", "metaG_mining"))

# remove column study accession

metaG_mining_short <- metaG_mining_short[c("sample_accession", "lat","lon", "scientific_name","DB")]

new_rows <- tibble(
  sample_accession = c("SAMEA116057598", "MSM", "EMB"),
  lat = as.numeric(c("54.087639", "-25", "54.2573")),  # Convert to numeric
  lon = as.numeric(c("7.968194", "14.38333333", "14.328883")),  # Convert to numeric
  scientific_name = rep("marine sediment metagenome", 3),
  DB = rep("EX",3)
)

metaG_mining_short <- bind_rows(metaG_mining_short, new_rows)

# Combine OMD_combined and metaG_mining_short. Select all from metaG mining and then only add those sample_accessions that are not yet included in metaG

# Filter out rows from OMD_combined where sample_accession is already in metaG_mining_short
OMD_tmp <- OMD_combined %>%
  filter(!sample_accession %in% metaG_mining_short$sample_accession)

# Combine the original metaG_mining_short with the filtered OMD_combined
OMD_metaG <- bind_rows(metaG_mining_short, OMD_tmp)

scientific_name_plot <- OMD_metaG
scientific_name_plot$scientific_name <- gsub(" metagenome", "", scientific_name_plot$scientific_name)

unique(scientific_name_plot$scientific_name)

keywords <- c("marine", "freshwater", "seawater", 
              "sediment", "lake water", "aquatic", 
              "marine plankton", "freshwater sediment", "marine sediment",
              "sponge", "salt lake", "hot springs",
              "microbial mat", "hydrothermal vent", "estuary",
              "riverine", "subsurface", "mine drainage","cold seep")

scientific_name_plot <- scientific_name_plot %>%
  mutate(environment = case_when(
    scientific_name %in% keywords ~ scientific_name, 
    TRUE ~ "other"
  ))

scientific_name_plot$DB <- factor(
  scientific_name_plot$DB, 
  levels = c("EX","OMD","metaG_mining") ,ordered = T)

scientific_name_plot_df <- scientific_name_plot %>%
  group_by(environment, DB) %>%
  summarise(sample_count = n(), .groups = "drop")

scientific_name_plot_df$environment <- factor(
  scientific_name_plot_df$environment, 
  levels = c("marine", "freshwater", "seawater", 
             "sediment", "lake water", "aquatic", 
             "marine plankton", "freshwater sediment", "marine sediment",
             "sponge", "salt lake", "hot springs",
             "microbial mat", "hydrothermal vent", "estuary",
             "riverine", "subsurface", "mine drainage","cold seep","other") ,ordered = T)

scientific_name_plot_df$DB <- factor(
  scientific_name_plot_df$DB, 
  levels = c("EX","OMD","metaG_mining") ,ordered = T)

cols <- c("#78c679","#b3b3b3","#797979")
# plot statistics of runs and reads ####
ggplot(scientific_name_plot_df, aes(x = environment, y=sample_count,fill = DB)) +
  geom_bar(stat = "identity",width = 0.6)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 7, angle=45, hjust = 1, color = "black"),
        legend.title = element_text(size=7, face="bold", color="black"),
        legend.text = element_text(size=7, color="black"),
        axis.title = element_text(size = 7, face = "bold", color="black"),
        axis.text.y = element_text(size = 7, color="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'))+
  scale_fill_manual(values = cols)+
  scale_y_continuous(expand = c(0,0), limits = c(0,7500), breaks = c(0,1000,2000,3000,4000,5000,6000,7000))+
  xlab("")+ ylab("metagenomic samples screened")

ggsave("Overview_stats/Bar_plot_searched_samples_revised.pdf", width = 150, height = 80, device="pdf", units = "mm")

# Plot world map with all samples in which EX MAGs were present
library(maps)
# plot 
World <- map_data("world")
ggplot() +
  geom_polygon(data = World, aes(x=long, y =lat, group = group), fill="#e8e8e8")+
  geom_point(data=scientific_name_plot[scientific_name_plot$DB =="OMD",],
             aes(x=lon, y=lat),color ="#b3b3b3", size = 0.4)+
  geom_point(data=scientific_name_plot[scientific_name_plot$DB =="metaG_mining",],
             aes(x=lon, y=lat),color ="#797979", size = 0.4)+
  geom_point(data=scientific_name_plot[scientific_name_plot$DB =="EX",],
             aes(x=lon, y=lat),color ="#238443",fill = "#78c679", size = 2, shape=24, stroke = 0.3)+
  labs(x="longitude", y="latitude")+
  theme(plot.title = element_text(size=7, face="bold"),
        legend.title = element_text(size=7, face="bold", color="black"),
        legend.text = element_text(size=7, color="black"),
        legend.position = "bottom",
        axis.title = element_text(size = 7, face = "bold", color="black"),
        axis.text = element_text(size = 7, color="black"),
        axis.ticks = element_line(linewidth = 0.3),
        legend.key=element_rect(fill="white"),
        panel.background = element_rect(fill = "white", color ="black", linewidth = 0.3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

ggsave("Overview_stats/Overview_all_studies_EX4484-6_MAGs_revised.pdf", width = 180, height = 100, device="pdf", units = "mm")


save.image("Overview_stats/Figure2_overview_stats_revised.RData")
