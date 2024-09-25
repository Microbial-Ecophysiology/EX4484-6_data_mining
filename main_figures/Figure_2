#Script for overall statistics based on meta data from screened metagenomic runs

library(tidyverse)
library(reshape2)
require(R.utils)
library(patchwork)
library(ggplot2)

setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Data_mining/")
load("Overview_stats/Figure1_overview_stats.RData")

# plots for statistics screened ####
df <- read.table(
  "subset_MG_scientific_name_8287_final.txt", 
  sep="\t",
  header=T, 
  comment.char = "",
  quote = "\"",
  check.names = FALSE
)

# check file size
hsize(sum(as.numeric(unlist(strsplit(df$fastq_bytes, ";")))))
# "57.8 TiB"

unique(df$scientific_name)

scientific_name <- df %>% 
  group_by(scientific_name) %>%
  summarise(count = length(scientific_name),
            base_count = sum(base_count)) %>%
  mutate(base_count_Gbp = base_count/1000000)

scientific_name <- select(scientific_name, c("scientific_name","count","base_count_Gbp"))
scientific_name_plot <- melt(scientific_name)

scientific_name_plot$scientific_name <- gsub(" metagenome", "", scientific_name_plot$scientific_name)
scientific_name_plot$scientific_name <- factor(scientific_name_plot$scientific_name, levels = c("marine","freshwater",
                                                                                                "sediment","lake water",
                                                                                                "aquatic","freshwater sediment",
                                                                                                "seawater","marine sediment",
                                                                                                "salt lake","hot springs",
                                                                                                "estuary","microbial mat",
                                                                                                "riverine","mine drainage",
                                                                                                "subsurface","hydrothermal vent",
                                                                                                "hydrocarbon","pond",
                                                                                                "cold seep","biofilm",
                                                                                                "soda lake","hypersaline lake",
                                                                                                "mine tailings","lagoon",
                                                                                                "cold spring","glacier",
                                                                                                "mine","surface",
                                                                                                "alkali sediment","hypolithon") ,ordered = T)

# plot statistics of runs and reads ####
plot_1 <- ggplot(scientific_name_plot[scientific_name_plot$variable == "count",], aes(x = scientific_name, y=value)) +
  geom_bar(stat = "identity", position = position_dodge(1.0),  width = 0.6)+
  theme_bw()+
  theme(axis.text.x= element_blank(),
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
  scale_y_continuous(expand = c(0,0), limits = c(0,3200), breaks = c(0,500,1000,1500,2000,2500,3000))+
  xlab("")+ ylab("metagenomic runs \nscreened")


plot_2 <-ggplot(scientific_name_plot[scientific_name_plot$variable == "base_count_Gbp",], aes(x = scientific_name, y=value)) +
  geom_bar(stat = "identity", position = position_dodge(1.0),  width = 0.6)+
  theme_bw()+
  theme(axis.text.x= element_text(size = 7, angle=45, hjust = 1),
        legend.title = element_text(size=7, face="bold", color="black"),
        legend.text = element_text(size=7, color="black"),
        axis.title = element_text(size = 7,face = "bold", color="black"),
        axis.text.y = element_text(size = 7, color="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'))+
  scale_y_continuous(expand = c(0,0), limits = c(0,34000000))+
  xlab("")+ ylab("base count [Gbp]")



# plots for overview screened ####
# Starting table of all MGs we obtained through ENA
initial_MG <- read.table("results_MG_mining.txt", header = T, sep="\t", comment.char = "",
                         quote = "\"",
                         check.names = FALSE) 
#distinct samples within initial_MG
n_distinct(initial_MG$run_accession) # 44968
# distinct studies within initial_MG
n_distinct(initial_MG$study_accession) # 6816


# after subsampling
subsampled_MG <- read.table("subset_MG_scientific_name_8287_final.txt", header = T, sep="\t", comment.char = "",
                            quote = "\"",
                            check.names = FALSE)
# distinct samples within subsampled_MG
n_distinct(subsampled_MG$run_accession) # 8287
#distinct studies within subsampled_MG
n_distinct(subsampled_MG$study_accession) #2184


# include read counts with similarity 50 / aligned 50
bins50 <- read.table("CoverM/EX/results_coverm/EX_bin_count_aln_50.txt", header = T, row.names = 1)

# select all samples with >1000 reads in at least one MAG
bins50_subset <- bins50 %>%
  filter(E3_1_d157_spades_bin_16_orig_1_refined.contigs >= 1000 | GCA_003650025.1_ASM365002v1_genomic >= 1000 | GCA_003663595.1_ASM366359v1_genomic >= 1000 | GCA_016928095.1_ASM1692809v1_genomic >= 1000 | GCA_023145185.1_ASM2314518v1_genomic >= 1000)

bins50_merged <- merge(bins50_subset, subsampled_MG, by.x="row.names",by.y="run_accession")

#distinct samples within coverm_reads_MG
n_distinct(bins50_merged$Row.names) # 1758
#distinct studies within coverm_reads_MG
n_distinct(bins50_merged$study_accession) #392 of which 8 have MAGs of EX4484-6


# using coverage / study >2
cov50 <- read.table("CoverM/EX/results_coverm/EX_mean_cov_sim50_aln50.txt", header = T, row.names = 1)

# calculate coverage sum for each run_accession
cov50$mean_cov_sum <- rowSums(cov50[ , c(1:5)], na.rm=TRUE)
#merge tables by run_accessions
cov50_merged <- merge(cov50, subsampled_MG, by.x="row.names",by.y="run_accession")
#merge table by study_accession and select only studies with coverage >2
tmp <- cov50_merged[c(2:7,112)]
tmp2 <- tmp %>% group_by(study_accession) %>% summarise_each(list(sum))

cov50_subset <- tmp2 %>%
  filter(mean_cov_sum >= 2)
# 30 studies remain

#subsample all samples, which are included in subsampled studies
study_ids <- cov50_subset$study_accession
cov50_samples <- subset(cov50_merged, study_accession %in% study_ids) # 1972 
# not all samples from these studies had read counts >1000 and were therefore excluded ins bins50_subset
# for screening the single studies throughout our metagenomic analysis we included all samples from these studies for coassembly
# total number of samples: 2588 (see also Sl table 1)

# studies from OMD
OMD <- read.table("Overview_stats/genome_export_OMD.txt", header = T, sep="\t", comment.char = "",
              quote = "\"",
              check.names = FALSE)

# studies with EX_MAGs from metaG mining
EX_MAGs_metaG <- c("PRJNA368391","PRJNA531756","PRJNA541421","PRJNA704804","PRJNA721298","PRJNA889212")
EX_samples <- subset(cov50_merged, study_accession %in% EX_MAGs_metaG)


# studies with EX MAGs from MAG mining and from co authors
EX_MAGs_MAG <- data.frame(study = c("PRJNA362212","PRJNA480137","PRJNA690107","PRJNA707313","PRJNA713414","PRJNA544741","PRJNA326482","PRJEB4352","PRJNA776043","PRJNA875933","MSM105","EMB267","HMA"),
                          lat = c("29.52416667", "27.0078","34.99562", "22.16", "23.957043","10.617", "10.617","18.7241", "36.68", "42.15971392","-25", "54.2573","54.087627"),
                          lon =c("-113.57", "-111.4071","-98.68895", "119.29", "-108.862273","-65.167", "-65.167","66.3896", "100.47", "-62.360103","14.38333333", "14.328883","	7.968191"))


##############################
#Plot

df <- subsampled_MG[, c("lat", "lon", "run_accession")]
df <- df[complete.cases(df),]                                                     
df$group <- factor("sreened")


df2 <- subsampled_MG %>% 
  select(lat, lon, run_accession) %>% 
  filter(complete.cases(.)) %>% 
  mutate(group = factor(if_else(run_accession %in% EX_samples$Row.names, "metaG studies including EX MAGs", 
                                if_else(run_accession %in% cov50_samples$Row.names, "studies with coverage >2", "screened")),
                        # levels = c("EX", "coverage", "subsampled", "initial")))
                        levels = c("screened", "studies with coverage >2", "metaG studies including EX MAGs"),
                        labels = c("screened metagenomic runs","metagenomic runs with coverage of EX4484-6 MAGs >2", "metagenomic runs including EX4484-6 MAGs"))) %>% 
  arrange(group, .by_group = T)

# plot 

# Plot world map with samples
library(maps)

World <- map_data("world")

cols = c("#969696","#6baed6","#084594")

ggplot() +
  geom_polygon(data = World, aes(x=long, y =lat, group = group), fill="grey")+
  geom_point(data=df2,  
             aes(x=lon, y=lat, color = group, size = group)) +
  scale_size_manual(values=c(0.5,0.5,1.5))+
  scale_color_manual(values = cols)+
  labs(x="longitude", y="latitude", size= "metagenomic runs", color="metagenomic runs")+
  ggtitle("Overview of searched metagenomic runs")+
  theme(plot.title = element_text(size=7, face="bold"),
        legend.title = element_text(size=7, face="bold", color="black"),
        legend.text = element_text(size=7, color="black"),
        axis.title = element_text(size = 7, color="black"),
        axis.text = element_text(size = 7, color="black"),
        legend.key=element_rect(fill="white"))

ggsave("Overview_stats/verview_searched_MGs_subsampled.png", width = 180, height = 80, device="png", units = "mm")

# only subsample in screened and samples with found EX MAGs
df3 <- subsampled_MG %>% 
  select(lat, lon, run_accession) %>% 
  filter(complete.cases(.)) %>% 
  mutate(group = factor(if_else(run_accession %in% EX_samples$Row.names, "metaG EX MAGs","screened"),
                        levels = c("screened", "metaG EX MAGs"))) %>% 
  arrange(group, .by_group = T)


# Plot world map with all samples in which EX MAGs were present

cols_screening = c("#737373", "#ce1256")


# plot 

EX_MAGs_map <- ggplot() +
  geom_polygon(data = World, aes(x=long, y =lat, group = group), fill="grey")+
  geom_point(data=OMD,
             aes(x=longitude, y=latitude), color = "#737373", size = 0.4) +
  geom_point(data=df3,  # these contain sites of studies with coverage >2, in which EX MAGs could be found
             aes(x=lon, y=lat, color = group, size = group)) +
  geom_point(data = EX_MAGs_MAG, # this is the EX_MAGs defined previously containing MAGs from ENA assembly mining and co-authors
             aes(x = as.numeric(lon), y = as.numeric(lat), color = "metaG EX MAGs", size = "metaG EX MAGs")) +
  scale_size_manual(values=c(0.4, 1.25),
                    breaks = c("screened","metaG EX MAGs"),
                    labels = c("screened", "containing medium and \nhigh quality EX4484-6 MAGs"))+
  scale_color_manual(values = cols_screening, 
                     breaks = c("screened","metaG EX MAGs"),
                     labels = c("screened", "containing medium and \nhigh quality EX4484-6 MAGs"))+
  labs(x="longitude", y="latitude", size= "metagenomic runs", color="metagenomic runs")+
  #ggtitle("Location of studies containing EX4484-6 MAGs")+
  theme(plot.title = element_text(size=7, face="bold"),
        legend.title = element_text(size=7, face="bold", color="black"),
        legend.text = element_text(size=7, color="black"),
        legend.position = "bottom",
        axis.title = element_text(size = 7, face = "bold", color="black"),
        axis.text = element_text(size = 7, color="black"),
        legend.key=element_rect(fill="white"),
        panel.background = element_rect(fill = "#f0f0f0"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )+ 
  guides(colour = guide_legend(override.aes = list(size=2)))

ggsave("Overview_stats/Overview_all_studies_EX4484-6_MAGs.png", width = 180, height = 80, device="png", units = "mm")



# Include data about single steps for stacked bar plots
# from files we can see: total number of screened metagenomic runs: 8287
# screened MAGs for studies with cumulative coverage >2: 2588 (Sl Tab1)
# redundant EX4484-6 MAGs found from screening: 17 MAGs, of which 4 were already present through MAG mining,
# in 8 studies containing 141 samples

# Create input table for plot ####

data_metaG <- data.frame(
  group=c("screened","studies with cumulative coverage >2", "studies containing EX4484-6 MAGs"),
  total=c(8287, 2588, 141)
)

# Calculate the proportion of screened MAGs at each step (all screened: 8287-2588, etc.)
data_metaG$proportion_screened <- c(5699,2447,141)


data_metaG$group <- factor(data_metaG$group, levels=c("screened","studies with cumulative coverage >2", "studies containing EX4484-6 MAGs"),
                           labels = c("screened","runs of studies with \ncumulative coverage >2", "runs of studies containing \nEX4484-6 MAGs"),
                           ordered = T)
# barplot metaG

cols2= c( "#969696","#6baed6","#084594")

metaG <- ggplot(data_metaG, aes(x = "", y = proportion_screened, fill = group)) +
  geom_bar(stat = "identity") +
  labs(y = "metagenomic runs \n(n = 8287)",
       fill = "")+
  theme_bw() +
  scale_fill_manual("", values = cols2) +
  theme(legend.position = "right")+
  theme(axis.title.y = element_text(face="bold", size=7, color="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1, size = 7, color="black"),
        axis.text.y = element_text(size = 7, color="black"),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 7, color="black"),
        legend.title = element_text(size = 7, face = "bold", color="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.key.size = unit(0.5,'cm'))+
  scale_y_continuous(expand= c(0,0),  breaks = scales::pretty_breaks(n = 6))+ 
  guides(colour = guide_legend(override.aes = list(size=2)))


# barplot MAG mining

# Create input table for plot ####
# for this plot the GDTB retrieved MAGs are added to each of the groups: screened, quality filtered, taxonomy filtered, as these have been screened as well, 
# see count_data_mining_steps.R script

data_MAG <- data.frame(
  group=c("screened","quality filtered", "taxonomy filtered", "dereplicated", "EX4484-6 MAGs"),
  total=c(11479, 1734, 844, 388, 5)
)

# Calculate the proportion of screened MAGs at each step (all screened: 11479-1734, etc.)
data_MAG$proportion_screened <- c(9745, 890, 456, 383, 5)


data_MAG$group <- factor(data_MAG$group, levels=c("screened","quality filtered", "taxonomy filtered", "dereplicated", "EX4484-6 MAGs"),
                         ordered = T)

# barplot

cols_MAG= c( "#969696","#238443","#78c679","#c2e699","#ffffcc")

MAG <- ggplot(data_MAG, aes(x = "", y = proportion_screened, fill = group)) +
  geom_bar(stat = "identity") +
  labs(y = "MAGs (n = 11479)",
       fill = "")+
  theme_bw() +
  scale_fill_manual("", values = cols_MAG) +
  theme(legend.position = "right")+
  theme(axis.title.y = element_text(face="bold", size=7, color="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.25, hjust=1, size = 7, color="black"),
        axis.text.y = element_text(size = 7, color="black"),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 7, color="black"),
        legend.title = element_text(size = 7, face = "bold", color="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.key.size = unit(0.5,'cm'))+
  scale_y_continuous(expand= c(0,0),  breaks = scales::pretty_breaks(n = 8))

#combine plots to one ####

plot_high <- ((MAG / metaG ) | ( plot_1  / plot_2))+
  plot_layout(widths = c("1","5"))&  
  theme(legend.justification = "left")

plot_high / EX_MAGs_map &  
  plot_annotation(tag_levels = "a")& 
  theme(plot.tag = element_text(size = 7, face = "bold"))

ggsave("Overview_stats/Plot_all_screened_map_plusOMD.png", width = 180, height = 210, device="png", units = "mm")
ggsave("Overview_stats/Plot_all_screened_map_plusOMD.emf", width = 180, height = 210, device="emf", units = "mm")
ggsave("Overview_stats/Plot_all_screened_map_plusOMD.pdf", width = 180, height = 210, device="pdf", units = "mm")


save.image("Overview_stats/Figure1_overview_stats.RData")
