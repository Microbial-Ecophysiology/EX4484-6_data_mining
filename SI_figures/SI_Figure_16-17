# Overview of predicted and hypothetical genes within enrichments

# ## on olorin
# printf "predicted_genes\tgenome\n" > /storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs/predicted_genes_count.txt
# # #
# ls -1 /storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs/mmseq2/prodigal_out/*.faa | xargs -n1 basename | sed "s/\.faa//" | while read line
# do
#  grep '^>' -c /storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs/mmseq2/prodigal_out/${line}".faa"  | sed -e "s/$/\t${line}/"
# done >> /storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs/predicted_genes_count.txt
# # 
# ls -1 /storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs/ETH_MAGs/ETH_selection_final/prodigal/*.faa | xargs -n1 basename | sed "s/\.faa//" | while read line
# do
#    grep '^>' -c /storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs/ETH_MAGs/ETH_selection_final/prodigal/${line}".faa"  | sed -e "s/$/\t${line}/"
# done >> /storage/hdd2/mmaeke/Metagenomics/PhD/ENA_data_mining/EX_MAGs/predicted_genes_count.txt

# copy predicted_genes_count to local 

library(dplyr)
library(ggplot2)
library(patchwork)

setwd("C:/Dokumente und Einstellungen/admin/OneDrive/PhD/Data_mining/Annotation_v3/hypothetical/")

predicted_genes <- read.table("predicted_genes_count.txt",
                              header = T,
                              sep = "\t")
predicted_genes$genome <- gsub(".genes","",as.character(predicted_genes$genome))

annotation <- read.table("../annotation_parsed.txt", 
                         sep="\t", 
                         header=TRUE, 
                         fill=TRUE, 
                         quote="")
annotation <- annotation %>% filter(cluster != "qseqid")
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

# subset table with unknown and hypothetical genes
annotation_hypothetical <- annotation[annotation$annotation_status == "hypothetical",]
write.table(annotation_hypothetical, "hypothetical_genes_EX4484-6.txt", row.names = F, col.names = T, sep = "\t", quote = F)


annotation_summary <- data.frame(do.call("rbind", by(annotation$annotation_status, annotation$genome, table)))
all.equal(predicted_genes$genome, rownames(annotation_summary))
annotation_summary$unknown <- predicted_genes$predicted_genes - rowSums(annotation_summary[, 1:2])

unknown_percent <- annotation_summary

mean(unknown_percent$percentage_unknown, na.rm = TRUE) # 26.47821

unknown_percent <- unknown_percent %>%
  mutate(percentage_annotated = (annotated / (annotated + hypothetical + unknown)) * 100)%>%
  mutate(percentage_hypothetical = (hypothetical / (annotated + hypothetical + unknown)) * 100)%>%
  mutate(percentage_unknown = (unknown / (annotated + hypothetical + unknown)) * 100)

write.table(unknown_percent, "unknown_percent_EX4484-6.txt", row.names = T, col.names = T, sep = "\t", quote = F)




# Stacked bar plot
library(tidyr)

# Add a new column for genome names
annotation_summary$genome <- rownames(annotation_summary)

# Reset row names to NULL to avoid issues with gather()
rownames(annotation_summary) <- NULL

#write table with data
write.table(annotation_summary, "summary_genes_EX4484-6.txt", 
            quote = FALSE,
            sep="\t",
            col.names = TRUE,
            row.names = FALSE)

# Reshape the data to a long format for counts
annotation_summary_long <- gather(annotation_summary, key = "group", value = "count", -genome)


annotation_summary_long$group <- factor(annotation_summary_long$group, 
                                        levels = c("unknown","hypothetical","annotated"),
                                        ordered = T)


annotation_summary_long$genome <- factor(annotation_summary_long$genome, levels = c("E3_1_d157_spades_bin_16_orig_1_refined-contigs","PRJNA889212_bin_7_strict_refined-contigs","PRJNA541421_bin_70_strict-contigs","MSM105_N25025F_bin_277_ori_permissive_1",
                                                                                    "GCA_003663595.1_ASM366359v1_genomic","PRJNA368391_bin_183_orig-contigs","EMB267_Co_bin_434_strict_1",
                                                                                    "PRJNA721298_bin_46_orig_refined-contigs","PRJNA541421_bin_72_orig_refined-contigs","PRJNA889212_bin_1_orig_refined-contigs","GCA_016928095.1_ASM1692809v1_genomic",
                                                                                    "PRJNA531756_bin_93_strict_refined-contigs","GCA_023145185.1_ASM2314518v1_genomic","ZORZ22-1_SAMN30647027_MAG_00000060",
                                                                                    "ZHEN22-1_SAMN22703513_MAG_00000265","ZHEN22-1_SAMN22703514_MAG_00000188","ZHEN22-1_SAMN22703512_MAG_00000312",
                                                                                    "GCA_011039765.1_ASM1103976v1_genomic","GCA_003650025.1_ASM365002v1_genomic","PRJNA368391_bin_90_orig-contigs","GCA_003663625.1_ASM366362v1_genomic",
                                                                                    "GCA_021158585.1_ASM2115858v1_genomic","GCA_021160765.1_ASM2116076v1_genomic","PRJNA531756_bin_72_orig_refined-contigs","PRJNA541421_bin_98_orig_refined-contigs",
                                                                                    "TARA_SAMEA2620113_MAG_00000097","PRJNA704804_bin_147_orig_refined-contigs","SCRA20-1_SAMN11854494_MAG_00000079","SUTE22-1_SAMN10231904_MAG_00000012",
                                                                                    "SUTE22-1_SAMN10231913_MAG_00000177","SUTE22-1_SAMN10231914_MAG_00000292","SUTE22-1_SAMN10231893_MAG_00000067","SUTE22-1_SAMN10231913_MAG_00000219",
                                                                                    "SUTE22-1_SAMN10231894_MAG_00000173","SUTE22-1_SAMN10231914_MAG_00000186"),
                                         labels = c("E3_1_d157","PRJNA889212_2","PRJNA541421_1","MSM105","GCA_003663595.1","PRJNA368391_1","EMB267",
                                                    "PRJNA721298_1","PRJNA541421_2","PRJNA889212_1","GCA_016928095.1","PRJNA531756_2",
                                                    "GCA_023145185.1","SAMN30647027_1","SAMN22703513_1","SAMN22703514_1","SAMN22703512_1",
                                                    "GCA_011039765.1","GCA_003650025.1","PRJNA368391_2","GCA_003663625.1",
                                                    "GCA_021158585.1","GCA_021160765.1","PRJNA531756_1","PRJNA541421_3",
                                                    "SAMEA2620113_1","PRJNA704804_1","SAMN11854494_1","SAMN10231904_1",
                                                    "SAMN10231913_2","SAMN10231914_2", "SAMN10231893_1","SAMN10231913_1", "SAMN10231894_1","SAMN10231914_1"), ordered = TRUE)  


# plot data
cols = c("#e0e0e0","#92c5de","#053061")


GENES <- ggplot(annotation_summary_long, aes(x = genome, y = count, fill = group)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "",
       y = "number of genes",
       fill = "predicted genes") +
  theme_bw() +
  scale_fill_manual("predicted genes", values = cols) +
  theme(legend.position = "right")+
  theme(axis.title.y = element_text(face="bold", size=7, color="black"),
        axis.title.x = element_text(face="bold", size=7, color="black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7, color="black"),
        axis.text.y = element_text(size = 7, color="black"),
        legend.text = element_text(size = 7, color="black"),
        legend.title = element_text(size = 7, face = "bold", color="black"),
        legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'))+
  scale_y_continuous(expand= c(0,0), breaks = scales::pretty_breaks(n = 10))

GENES_rel_abd <- ggplot(annotation_summary_long, aes(x = genome, y = count, fill = group)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "",
       y = "percentage of genes [%]",
       fill = "predicted genes") +
  theme_bw() +
  scale_fill_manual("predicted genes", values = cols) +
  theme(legend.position = "right")+
  theme(axis.title.y = element_text(face="bold", size=7, color="black"),
        axis.title.x = element_text(face="bold", size=7, color="black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7, color="black"),
        axis.text.y = element_text(size = 7, color="black"),
        legend.text = element_text(size = 7, color="black"),
        legend.title = element_text(size = 7, face = "bold", color="black"),
        legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'))+
  scale_y_continuous(expand= c(0,0), breaks = scales::pretty_breaks(n = 10), labels = scales::percent)


GENES / GENES_rel_abd + 
  plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 7, face = "bold"))& 
  theme(legend.justification = "left")

ggsave("combined_gene_prediction_kegg_NR.png", width = 18, height = 22, units = "cm", bg='transparent')
ggsave("combined_gene_prediction_kegg_NR.emf", width = 18, height = 22, units = "cm", bg='transparent')
ggsave("combined_gene_prediction_kegg_NR.pdf", width = 18, height = 22, units = "cm", bg='transparent')




# search category for all hypothetical and unknown genes by using agnostos integraded DB summary
agnostos <- read.table("agnostos_integrated_DB_summary_information.tsv", 
                       sep="\t", 
                       header=TRUE, 
                       fill=TRUE, 
                       quote="")

# Find non unique genes in dataframe
non_unique_genes <- agnostos$gene_name[duplicated(agnostos$gene_name) | duplicated(agnostos$gene_name, fromLast = TRUE)]

# Filter the dataframe for the non-unique gene names
non_unique_df <- agnostos %>%
  filter(gene_name %in% non_unique_genes)

# check if same gene_name has the same category
result_df <- non_unique_df %>%
  group_by(gene_name) %>%
  summarise(same_category = n_distinct(category) == 1)

# Check if same_category always is TRUE
all_true <- all(result_df$same_category) #TRUE


# Sum up the counts for each unique genome and category
# Subset the dataframe for columns genome, gene_name, and category
agnostos_sub_df <- agnostos %>%
  select(genome, gene_name, category)

# Combine rows with duplicate gene_names
agnostos_combined_df <- agnostos_sub_df %>%
  group_by(genome, gene_name, category) %>%
  summarise()

# summarise single categories for each genome
agnostos_summed <- agnostos_combined_df %>%
  group_by(genome,  category) %>%
  summarise(count = n()) %>%
  ungroup()


agnostos_summed$genome <- factor(agnostos_summed$genome, levels = c("E3_1_d157_spades_bin_16_orig_1_refined-contigs","PRJNA889212_bin_7_strict_refined-contigs","PRJNA541421_bin_70_strict-contigs","MSM105_N25025F_bin_277_ori_permissive_1",
                                                                    "GCA_003663595.1_ASM366359v1_genomic","PRJNA368391_bin_183_orig-contigs","EMB267_Co_bin_434_strict_1",
                                                                    "PRJNA721298_bin_46_orig_refined-contigs","PRJNA541421_bin_72_orig_refined-contigs","PRJNA889212_bin_1_orig_refined-contigs","GCA_016928095.1_ASM1692809v1_genomic",
                                                                    "PRJNA531756_bin_93_strict_refined-contigs","GCA_023145185.1_ASM2314518v1_genomic","ZORZ22-1_SAMN30647027_MAG_00000060",
                                                                    "ZHEN22-1_SAMN22703513_MAG_00000265","ZHEN22-1_SAMN22703514_MAG_00000188","ZHEN22-1_SAMN22703512_MAG_00000312",
                                                                    "GCA_011039765.1_ASM1103976v1_genomic","GCA_003650025.1_ASM365002v1_genomic","PRJNA368391_bin_90_orig-contigs","GCA_003663625.1_ASM366362v1_genomic",
                                                                    "GCA_021158585.1_ASM2115858v1_genomic","GCA_021160765.1_ASM2116076v1_genomic","PRJNA531756_bin_72_orig_refined-contigs","PRJNA541421_bin_98_orig_refined-contigs",
                                                                    "TARA_SAMEA2620113_MAG_00000097","PRJNA704804_bin_147_orig_refined-contigs","SCRA20-1_SAMN11854494_MAG_00000079","SUTE22-1_SAMN10231904_MAG_00000012",
                                                                    "SUTE22-1_SAMN10231913_MAG_00000177","SUTE22-1_SAMN10231914_MAG_00000292","SUTE22-1_SAMN10231893_MAG_00000067","SUTE22-1_SAMN10231913_MAG_00000219",
                                                                    "SUTE22-1_SAMN10231894_MAG_00000173","SUTE22-1_SAMN10231914_MAG_00000186"),
                                 labels = c("E3_1_d157","PRJNA889212_2","PRJNA541421_1","MSM105","GCA_003663595.1","PRJNA368391_1","EMB267",
                                            "PRJNA721298_1","PRJNA541421_2","PRJNA889212_1","GCA_016928095.1","PRJNA531756_2",
                                            "GCA_023145185.1","SAMN30647027_1","SAMN22703513_1","SAMN22703514_1","SAMN22703512_1",
                                            "GCA_011039765.1","GCA_003650025.1","PRJNA368391_2","GCA_003663625.1",
                                            "GCA_021158585.1","GCA_021160765.1","PRJNA531756_1","PRJNA541421_3",
                                            "SAMEA2620113_1","PRJNA704804_1","SAMN11854494_1","SAMN10231904_1",
                                            "SAMN10231913_2","SAMN10231914_2", "SAMN10231893_1","SAMN10231913_1", "SAMN10231894_1","SAMN10231914_1"), ordered = TRUE)   

agnostos_summed$category <- factor(agnostos_summed$category, levels = c("DISC","EU","GU","KWP","K"),
                                   labels = c("discarded","environmental unknown","genomic unknown","known without Pfam \nannotation","known with Pfam \nannotation"), ordered = T)

# Create a stacked barplot
cols2 = c("#67001f","#d6604d","#fddbc7","#92c5de","#053061")

AGNOSTOS <- ggplot(agnostos_summed, aes(x = genome, y = count, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "",
       y = "number of genes",
       fill = "AGNOSTOS category") +
  theme_bw() +
  scale_fill_manual("AGNOSTOS category", values = cols2) +
  theme(legend.position = "right")+
  theme(axis.title.y = element_text(face="bold", size=7, color="black"),
        axis.title.x = element_text(face="bold", size=7, color="black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7, color="black"),
        axis.text.y = element_text(size = 7, color="black"),
        legend.text = element_text(size = 7, color="black"),
        legend.title = element_text(size = 7, face = "bold", color="black"),
        legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'))+
  scale_y_continuous(expand= c(0,0), breaks = scales::pretty_breaks(n = 10))

AGNOSTOS_rel_abd <- ggplot(agnostos_summed, aes(x = genome, y = count, fill = category)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "",
       y = "percentage of genes [%]",
       fill = "AGNOSTOS category") +
  theme_bw() +
  scale_fill_manual("AGNOSTOS category", values = cols2) +
  theme(legend.position = "right")+
  theme(axis.title.y = element_text(face="bold", size=7, color="black"),
        axis.title.x = element_text(face="bold", size=7, color="black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 7, color="black"),
        axis.text.y = element_text(size = 7, color="black"),
        legend.text = element_text(size = 7, color="black"),
        legend.title = element_text(size = 7, face = "bold", color="black"),
        legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'))+
  scale_y_continuous(expand= c(0,0), breaks = scales::pretty_breaks(n = 10), labels = scales::percent)

AGNOSTOS / AGNOSTOS_rel_abd + 
  plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 7, face = "bold"))& 
  theme(legend.justification = "left")

ggsave("combined_gene_prediction_AGNOSTOS.png", width = 18, height = 22, units = "cm", bg='transparent')
ggsave("combined_gene_prediction_AGNOSTOS.emf", width = 18, height = 22, units = "cm", bg='transparent')
ggsave("combined_gene_prediction_AGNOSTOS.pdf", width = 18, height = 22, units = "cm", bg='transparent')



load("predicted_genes.RData")

save.image("predicted_genes.RData")
