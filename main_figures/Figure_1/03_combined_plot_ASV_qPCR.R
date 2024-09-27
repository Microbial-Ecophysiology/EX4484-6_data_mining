# Create combined plot of 16S Amplicon data and qPCR data

library(tidyverse)
library(patchwork)
library(ggplot2)

setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Manuscript_Thermo_ENA/Plots/")

# load qPCR data ####
load("Fig_1_qPCR_data/qPCR_plot.RData")

# Plot Thermoplasmatota and Lokiarchaeia

cols2=c("#4393c3","#a6dba0")

qPCR.plot <- ggplot(qPCR[qPCR$Sample %in% c("Control","Protein") & qPCR$Taxa %in% c("EX4484-6",
                                                                                    "Lokiarchaeia (Loki-2b)") & qPCR$Day != "372",], 
                    aes(x = Day, y = value)) + 
  geom_bar(aes(color = Taxa, fill = Taxa),       
           stat = "identity", position = position_dodge(0.9),  width = 0.80) +
  geom_errorbar(aes(x=Day, y=value, ymax=value+SD, ymin=value-SD, group=interaction(Sample, Taxa)), color="black", 
                lty="solid",position = position_dodge(0.9), width=0.4, linewidth=0.25)+
  ggtitle("")+
  labs( x = "Day", y = "gene copies / mL slurry")+
  scale_color_manual("Taxa",values=cols2)+
  scale_y_continuous(expand=c(0,0),limits = c(0, 15000000), labels = scales::scientific)+
  scale_fill_manual("Taxa", values=cols2)+
  guides(color = guide_legend(ncol=1))+
  theme(title = element_text(size=7, face="bold"),
        axis.text.x = element_text(size=7, color = "black"),
        axis.text.y =element_text(size=7, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size=7, face="bold", color = "black"),
        axis.title.y = element_text(size=7, face="bold", color = "black"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size=7, face="bold"),
        legend.text = element_text(size=6),
        strip.background = element_blank(), 
        strip.text = element_text(size = 7, face="bold"))+
  facet_grid(~ Sample)



# load 16S Amplicon data ####
load("/Dokumente und Einstellungen/admin/OneDrive/PhD/Nadja_BSc/16S_Amplicon/dada2_analysis_new_14102023/Plot_rel_abundance_Archaea/arc_mod_GTDB_long_rel_abd_table_bel1perc.RData")


### select which classes will be plotted
unique(ta1$ASV_table_taxa_abv_1_perc$Class_mod)

### subset to only column with class_mod as single taxonomy and Phylum_mod to group them by
ta1_class <- select(ta1$ASV_table_taxa_abv_1_perc, OTU, Sample, Class_mod, Phylum_mod, Abundance, Incubation_time, Sediment_depth, AB,  BES, C_source) %>% 
  mutate(Class_mod=factor(Class_mod, levels=c("Heimdallarchaeia","Lokiarchaeia","Thorarchaeia","other_p_Asgardarchaeota_<1%", "Methanosarcinia", "Syntropharchaeia", "other_p_Halobacteriota_<1%",
                                              "Nanoarchaeia","E2","EX4484-6","Thermoplasmata","other_p_Thermoplasmatota_<1%","Bathyarchaeia","other_p_Thermoproteota_<1%" , "other_Archaea_<1%"), 
                          labels = c("Heimdallarchaeia","Lokiarchaeia","Thorarchaeia","other Asgardarchaeota", "Methanosarcinia", "Syntropharchaeia", "other Halobacteriota",
                                     "Nanoarchaeia","E2","EX4484-6","Thermoplasmata","other Thermoplasmatota","Bathyarchaeia","other Thermoproteota" , "other Archaea"),
                          ordered = T), 
         Phylum_mod=factor(Phylum_mod)) %>% 
  group_by(Class_mod)



### sum ASVs of same taxonomy
ta1_class_s <- ta1_class %>% 
  group_by(Sample, Class_mod, Phylum_mod, Incubation_time, Sediment_depth, AB,  BES, C_source) %>% 
  summarise(Abundance = sum(Abundance), .groups = "keep") %>% 
  ungroup()


# Relabel the samples based on their day and replicate *(replicate 2)
ta1_class_s$Sample <- factor(ta1_class_s$Sample, levels = c("A3.1.day98_1","A3.2.day98_1","A3.1.day157_1","A3.2.day157_1",
                                                            "E3.1.day98","E3.2.day98","E3.1.day157","E3.2.day157"),
                             labels = c("98", "98*", "157", "157*",
                                        "98", "98*", "157", "157*"), ordered = T)



# color vector, manually selected
cols_all = c("#00441b","#a6dba0","#1b7837","#e0e0e0","#cab2d6","#40004b","#e0e0e0","#fde0dd","#abd9e9","#4393c3","#081d58","#e0e0e0","#ffffb3","#e0e0e0","#969696")


# plot only replicate one of Control and Protein amended samples for manuscript
Amplicon.plot <- ggplot(ta1_class_s[ta1_class_s$C_source %in% c("Control","Protein") & ta1_class_s$BES == "noBES" & ta1_class_s$AB == "AB" & ta1_class_s$Sediment_depth =="upper" & ta1_class_s$Sample %in% c("98","157"),], aes(x = Sample, y = Abundance, fill = Class_mod)) + 
  geom_bar(position = "fill", stat = "identity", col = "black", linewidth = 0.1) + 
  scale_fill_manual("Taxa", values = cols_all) + theme_bw() +
  theme(plot.title = element_text(size = 7, face = "bold",hjust=0),
        axis.text.x=element_text(size=7, hjust = 1, vjust = 0.5, angle = 90, color = "black"),
        axis.text.y=element_text(size=7, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face="bold", size=7, color = "black"),
        axis.title.x = element_text(face="bold", size=7, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=6),
        legend.title=element_text(size=7, face="bold"),
        legend.key.size = unit(0.35, "cm"),
        strip.background = element_blank(), 
        strip.text = element_text(size = 7, face="bold")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  labs(x="Day", y = "relative abundance [%]")+
  guides(fill=guide_legend(ncol=1)) +
  facet_grid(~ C_source,scales = "free")


combined_EX_abundance <- Amplicon.plot + qPCR.plot + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 7, face = "bold"))

ggsave("Figure1_AbundanceEX.png", width = 18, height = 9, units = "cm")
ggsave("Figure1_AbundanceEX.pdf", width = 18, height = 9, units = "cm")


save.image("Figure1_AbundanceEX.RData")
