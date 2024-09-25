library(ggplot2)
library(reshape2)
library(tidyr)
library(patchwork)
library(gridExtra)
library(ggpubr)


setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/qPCR/Protein_degradation_Thermo_Loki/")

# load data ####
qPCR <- read.table("Generation_2_qpcr_results.txt", sep="\t", header=T, check.names = FALSE)
names(qPCR) <- c("Sample","Taxa","Day","value", "SD")

str(qPCR)

## order factors ####
qPCR$Sample <- factor(qPCR$Sample, levels=c("A3_1","A3_2",
                                            "E3_1","E3_2"),
                      labels=c("Control","Control*", "Protein","Protein*"))

qPCR$Day <- factor(qPCR$Day, levels=c("98",
                                      "157",
                                      "372"), ordered= T)

qPCR$Taxa <- factor(qPCR$Taxa, levels=c("Bacteria",
                                        "Archaea",
                                        "Thermoplasmatota",
                                        "Loki-2b"),
                    labels = c("Bacteria",
                               "Archaea",
                               "EX4484-6",
                               "Lokiarchaeia (Loki-2b)"),ordered= T)


# Plot all taxa for both replicates

cols=c("#e0e0e0","#969696","#4393c3","#a6dba0")

ggplot(qPCR, 
       aes(x = Day, y = value)) + 
  geom_bar(aes(color = Taxa, fill = Taxa),       
           stat = "identity", position = position_dodge(0.9),  width = 0.80) +
  geom_errorbar(aes(x=Day, y=value, ymax=value+SD, ymin=value-SD, group=interaction(Sample, Taxa)), color="black", 
                lty="solid",position = position_dodge(0.9), width=0.4, linewidth=0.25)+
  ggtitle("")+
  labs( x = "Day", y = "gene copies / mL slurry (log10 scale)")+
  scale_color_manual("Sample",values=cols)+
  scale_y_continuous(expand=c(0,0),limits = c(1, 1000000000), trans="log10",breaks=c(1,10,100,1e03,1e04,1e05,1e06,1e07,1e08,1e09), labels = scales::scientific)+
  scale_fill_manual("Sample", values=cols)+
  guides(color = guide_legend(ncol=1))+
  theme(title = element_text(size=7, face="bold"),
        axis.text.x = element_text(size=7,colour = "black"),
        axis.text.y =element_text(size=7,colour = "black"),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size=7, face="bold",colour = "black"),
        axis.title.y = element_text(size=7, face="bold",colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(size=7, face="bold"),
        legend.text = element_text(size=7),
        strip.background = element_blank(), 
        strip.text = element_text(size = 7, face="bold"))+
  facet_grid(~ Sample)

ggsave("../../Manuscript_Thermo_ENA/Plots/Fig_1_qPCR_data/qPCR_all_samples_Sl.png", width = 18, height = 9, device="png", units = "cm")
ggsave("../../Manuscript_Thermo_ENA/Plots/Fig_1_qPCR_data/qPCR_all_samples_Sl.pdf", width = 18, height = 9, device="pdf", units = "cm")
