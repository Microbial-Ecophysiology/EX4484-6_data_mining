library(castor)
library(ape)

setwd("C:/Dokumente und Einstellungen/admin/OneDrive/PhD/Data_mining/phylogenetic_tree/marker_gene/marker_gene_v3/")

tree <- read.tree("marker_gene_v3_newick_rooted.txt")
is.rooted(tree) #TRUE

#calculate relative evolutionary divergence (RED)
REDs = get_reds(tree)

tree_plot <- tree
tree_plot$node.label <- ifelse(
  tree$node.label == "",
  round(REDs, 3),
  paste0(tree$node.label, ": ", round(REDs, 3))
)

pdf("marker_gene_RED.pdf", width = 10, height = 30)
plot(tree_plot, show.node.label = T, show.tip.label = F, cex = 0.5, edge.color = "grey", edge.width = 0.5)
dev.off()
