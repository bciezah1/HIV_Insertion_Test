
# reading in data from deseq2
df = read.table("C:/Users/bciez/Documents/Basilio/UCSF/HIV_Insertion_Test/data/a1_table_ranked_by_most_freq_using60ENCONDES.csv", header=TRUE)
head(df)

# we want the log2 fold change 
original_gene_list <- df$span_count

# name the vector
names(original_gene_list) <- df$gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list[1:20]


# Load libraries
library(fgsea)
library(msigdbr)


####################################################
################ for HIV ###################
###########################################

# Load Reactome pathways from MSigDB (for the human genome)
msigdbr_genesets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

# Prepare the gene sets for fgsea
pathways <- split(msigdbr_genesets$gene_symbol, msigdbr_genesets$gs_name)

# Run fgsea
gene_list
fgsea_results <- fgsea(pathways = pathways, stats = gene_list, minSize = 15, maxSize = 500, nperm = 10000)

# Check results
fgsea_results

str(fgsea_results)

# Convert list columns (e.g., leadingEdge) into strings
fgsea_results$leadingEdge <- sapply(fgsea_results$leadingEdge, function(x) paste(x, collapse = ";"))

# Now write the result to a CSV
write.csv(fgsea_results, "C:/Users/bciez/Documents/Basilio/UCSF/HIV_Insertion_Test/data/gsea_results_span_count_C2_60encondes.csv", row.names = FALSE)
fgsea_results

####################################
############## plots ###############
####################################

# Extract the top pathways (e.g., top 10)
topPathways <- fgsea_results[order(pval), ][1:10, ]

# Plot a barplot of the top pathways
library(ggplot2)
ggplot(topPathways, aes(reorder(pathway, NES), NES)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 Enriched Pathways", x = "Pathway", y = "Normalized Enrichment Score (NES)") +
  theme_minimal()

# Plot enrichment for a specific pathway (e.g., "REACTOME_HIV_INFECTION")
plotEnrichment(pathway = pathways[["REACTOME_HIV_INFECTION"]], stats = gene_list) +
  labs(title = "Enrichment of HIV Infection Pathway")

# Save the results as a CSV file
write.csv(fgsea_results, "C:/Users/bciez/Documents/Basilio/UCSF/output_test.csv", row.names = FALSE)

# Print top pathways
head(fgsea_results[order(padj), ][, c("pathway", "NES", "pval", "padj", "size")])

#NES vs p-value Scatter Plot
ggplot(fgsea_results, aes(x = NES, y = -log10(pval))) +
  geom_point(aes(color = pval < 0.05)) +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "NES vs -log10(p-value)", x = "Normalized Enrichment Score", y = "-log10(p-value)") +
  theme_minimal()

# summary table

# Select top pathways with padj < 0.05
significant_pathways <- fgsea_results[fgsea_results$pval < 0.1, ]

# Print a table with the top significant pathways
head(significant_pathways[, c("pathway", "NES", "pval", "padj", "size")])

##############
ggplot(fgsea_results, aes(x = size, y = NES)) +
  geom_point(aes(color = pval < 0.1)) +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "Pathway Size vs NES", x = "Pathway Size", y = "Normalized Enrichment Score") +
  theme_minimal()

###### gene rank vs pathway enrichement

# Select a pathway (e.g., "REACTOME_HIV_INFECTION")
gene_ranks <- gene_list[names(pathways[["REACTOME_INFECTIOUS_DISEASE"]])]

# Plot the ranks
plot(gene_ranks, type = "h", main = "Gene Ranks for HIV Infection Pathway",
     xlab = "Gene Rank", ylab = "Ranked Score")


#enrichment heatmap
# Extract top 20 pathways based on padj
library(pheatmap)
top20_pathways <- fgsea_results[order(pval), ][1:10, ]

top20_pathways

# Create a heatmap
pheatmap(matrix(top20_pathways$NES, nrow = 1), 
         cluster_cols = FALSE, 
         labels_col = top20_pathways$pathway, 
         main = "NES for Top 20 Enriched Pathways",
         color = colorRampPalette(c("blue", "white", "red"))(100))

