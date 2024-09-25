# Script to investigate CAZY counts

library(tidyverse)
library(reshape2)

setwd("/Dokumente und Einstellungen/admin/OneDrive/PhD/Data_mining/Annotation_v3")

annotation <- read.table("annotation_parsed.txt", 
                         sep="\t", 
                         header=TRUE, 
                         fill=TRUE, 
                         quote="")

#count different CAZys within MAGs
CAZy_count <- aggregate(annotation$dbCAN_HMMER ~ annotation$genome, data = annotation, FUN = function(x) sum(!is.na(x)))
print(CAZy_count)


## kegg ----
kegg <- read.csv("full/concat_kegg_annotation_best.txt", 
                 sep="\t", 
                 header=TRUE, 
                 fill=TRUE, 
                 quote="",
                 col.names=c("qseqid", "kegg_sseqid","kegg_pident","kegg_length","kegg_mismatch","kegg_gapopen","kegg_qstart",
                             "kegg_qend","kegg_sstart","kegg_send","kegg_evalue","kegg_bitscore","kegg_bsr","kegg_KO","kegg_KO_name","kegg_taxon","kegg_gene"))
kegg <- subset(kegg, kegg_sseqid != "sseqid")

## NR ----
NR <- read.csv("full/concat_NR_annotation_best.txt", 
               sep="\t", 
               header=TRUE, 
               fill=TRUE, 
               quote="",
               col.names=c("qseqid", "NR_sseqid","NR_pident","NR_length","NR_mismatch","NR_gapopen","NR_qstart",
                           "NR_qend","NR_sstart","NR_send","NR_evalue","NR_bitscore","NR_qlen","NR_slen","NR_stitle","NR_staxids","NR_sscinames","NR_salltitles","NR_bsr"))
NR <- subset(NR, NR_sseqid != "sseqid")

## dbCAN ----
dbCAN <- read.csv("CAZy/concat_dbCAN_annotation.txt", 
                  sep="\t", 
                  header=TRUE, 
                  fill=TRUE, 
                  quote="",
                  col.names=c("Gene.ID","dbCAN_EC","dbCAN_HMMER","dbCAN_eCAMI","dbCAN_DIAMOND","dbCAN_no.Tools"))
dbCAN$Gene.ID <- gsub(".faa","",as.character(dbCAN$Gene.ID))





# combine table with Cazymes in first few columns
CAZy <- full_join(dbCAN, kegg, by = c("Gene.ID"="qseqid")) %>%
  full_join(., NR, by = c("Gene.ID"="qseqid")) %>%
  separate(col = Gene.ID, into = c("genome","cluster"), sep = "___")

# remove all rows from dataframe which don't have any annotation through NR or kegg
CAZy_short <- CAZy[complete.cases(CAZy[, 3:7]) & complete.cases(CAZy[, 8:41]), ]
# 852 CAZymes remain

write.table(CAZy_short, "CAZy/CAZy_short_parsed.txt", 
            quote = FALSE,
            sep="\t",
            col.names = TRUE,
            row.names = FALSE)
# manually refine CAzy by combining kegg, NR and CAZy counts. Only those CAZy's are used, that have at least one annotation through another annotation program

# use manually refined table to plot CAZy heatmap
CAZy_refined <- read.table("CAZy/CAZy_manual_refined.txt", 
                         sep="\t", 
                         header=TRUE, 
                         fill=TRUE, 
                         quote="")

# combine single CAZy families
CAZy_summed <- CAZy_refined %>%
  group_by(CAZy) %>%
  summarize_all(sum)

CAZy_summed_melt <- melt(CAZy_summed)
colnames(CAZy_summed_melt) <- c("CAZy","genome","value")

unique(CAZy_summed_melt$genome)
CAZy_summed_melt$genome <- factor(CAZy_summed_melt$genome, levels = c("SUTE22.1_SAMN10231914_MAG_00000186","SUTE22.1_SAMN10231894_MAG_00000173","SUTE22.1_SAMN10231913_MAG_00000219","SUTE22.1_SAMN10231893_MAG_00000067",
                                                                    "SUTE22.1_SAMN10231914_MAG_00000292","SUTE22.1_SAMN10231913_MAG_00000177",
                                                                    "SUTE22.1_SAMN10231904_MAG_00000012","SCRA20.1_SAMN11854494_MAG_00000079","PRJNA704804_bin_147_orig_refined.contigs","TARA_SAMEA2620113_MAG_00000097",
                                                                    "PRJNA541421_bin_98_orig_refined.contigs","PRJNA531756_bin_72_orig_refined.contigs",
                                                                    "GCA_021160765.1_ASM2116076v1_genomic","GCA_021158585.1_ASM2115858v1_genomic",
                                                                    "GCA_003663625.1_ASM366362v1_genomic","PRJNA368391_bin_90_orig.contigs","GCA_003650025.1_ASM365002v1_genomic","GCA_011039765.1_ASM1103976v1_genomic",
                                                                    "ZHEN22.1_SAMN22703512_MAG_00000312","ZHEN22.1_SAMN22703514_MAG_00000188","ZHEN22.1_SAMN22703513_MAG_00000265",
                                                                    "ZORZ22.1_SAMN30647027_MAG_00000060","GCA_023145185.1_ASM2314518v1_genomic", 
                                                                    "PRJNA531756_bin_93_strict_refined.contigs","GCA_016928095.1_ASM1692809v1_genomic","PRJNA889212_bin1_orig_refined.contigs", "PRJNA541421_bin_72_orig_refined.contigs","PRJNA721298_bin_46_orig_refined.contigs",
                                                                    "EMB267_Co_bin_434_strict_1","PRJNA368391_bin_183_orig.contigs","GCA_003663595.1_ASM366359v1_genomic",
                                                                    "MSM105_N25025F_bin_277_ori_permissive_1","PRJNA541421_bin_70_strict.contigs","PRJNA889212_bin7_strict_refined.contigs","E3_MAG"),
                                  labels = c("SAMN10231914_1", "SAMN10231894_1","SAMN10231913_1", "SAMN10231893_1",
                                            "SAMN10231914_2","SAMN10231913_2", 
                                            "SAMN10231904_1","SAMN11854494_1","PRJNA704804_1","SAMEA2620113_1",
                                            "PRJNA541421_3","PRJNA531756_1",
                                            "GCA_021160765.1","GCA_021158585.1",
                                            "GCA_003663625.1","PRJNA368391_2","GCA_003650025.1","GCA_011039765.1",
                                            "SAMN22703512_1","SAMN22703514_1","SAMN22703513_1",
                                            "SAMN30647027_1","GCA_023145185.1",
                                            "PRJNA531756_2","GCA_016928095.1","PRJNA889212_1","PRJNA541421_2","PRJNA721298_1",
                                            "EMB267","PRJNA368391_1","GCA_003663595.1","MSM105","PRJNA541421_1","PRJNA889212_2","E3_1_d157"),ordered =T)

# plot heatmap
ggplot(CAZy_summed_melt, aes(x = genome, y = CAZy)) +
  geom_point(aes(size=value, fill=value, col=value), shape = 21)+
  theme_bw() +
  coord_flip() +
  scale_y_discrete(position="right")+
  xlab("genome") + ylab("CAZy family")+
  theme(axis.text.x = element_text(angle=-90, hjust=1, vjust=-0.8, color='grey40', size=7),
        axis.text.y =element_text(color='grey40',size=7),
        legend.position = 'right',
        legend.title = element_text(size=7),
        legend.text = element_text(size=7),
        axis.title=element_text(face="bold", size=7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_color_viridis_c("CAZyme count",guide = "legend", option="mako", direction=-1)+
  scale_fill_viridis_c("CAZyme count",guide = "legend", option="mako", direction=-1)+
  scale_size_area("CAZyme count",max_size = 6)

ggsave("CAZy/CAZy_families_overview.png", plot = last_plot(), width = 18, height = 21, units = "cm", device = "png" )

