# Load necessary libraries
library(fgsea)
library(msigdbr)
library(ggplot2)
library(pheatmap)

##############################################
#
##############################################
# Set the working directory to your desired path
setwd("C:/Users/bciez/Documents/Basilio/UCSF/HIV_Insertion_Test/output")

# Verify the current working directory
getwd()


##############################################
# Step 1: Data Preparation
##############################################

# Read the input data from the DESeq2 output (gene list ranked by frequency)
df <- read.table("table_ranked_by_the_most_frequently_observed_gene_in_the_highest_number_of_patient_samples.txt", header=TRUE)
head(df)  # Preview the first few rows of the data

# Extract the span_count column (log2 fold change) to use as ranking values for GSEA
original_gene_list <- df$Npart

# Assign gene names to the ranking vector
names(original_gene_list) <- df$gene

# Remove any NA values from the gene list
gene_list <- na.omit(original_gene_list)

# Sort the gene list in decreasing order (required for GSEA)
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list[1:20]  # Preview the top 20 genes

##############################################
# Step 2: Load Pathways and Run GSEA
##############################################

# Load Reactome pathways for Homo sapiens from MSigDB (subcategory: CP:REACTOME)
msigdbr_genesets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

# Prepare the pathways for fgsea by splitting the data by gene set name (gs_name)
pathways <- split(msigdbr_genesets$gene_symbol, msigdbr_genesets$gs_name)

# Run fgsea on the gene list, with specified pathway size and permutations
fgsea_results <- fgsea(pathways = pathways, stats = gene_list, minSize = 15, maxSize = 500, nperm = 10000)

# Check the structure of the fgsea results
str(fgsea_results)

# Convert the leadingEdge column (which is a list) to a single string for easier readability
fgsea_results$leadingEdge <- sapply(fgsea_results$leadingEdge, function(x) paste(x, collapse = ";"))

# Save the fgsea results to a CSV file for future use
write.table(fgsea_results, "GSEA_MostFreqGenes_C2_60encondes_09_23_2024.txt", row.names = FALSE)

##############################################
# Step 3: Visualizing the Results
##############################################

# Extract the top 10 pathways based on p-values
topPathways <- fgsea_results[order(fgsea_results$pval), ][1:10, ]

# Plot a barplot showing the top 10 enriched pathways by their Normalized Enrichment Score (NES)
ggplot(topPathways, aes(reorder(pathway, NES), NES)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 Enriched Pathways", x = "Pathway", y = "Normalized Enrichment Score (NES)") +
  theme_minimal()

# Plot enrichment for a specific pathway (e.g., "REACTOME_HIV_INFECTION")
plotEnrichment(pathway = pathways[["REACTOME_HIV_INFECTION"]], stats = gene_list) +
  labs(title = "Enrichment of HIV Infection Pathway")

# Save the top results to a CSV file
#write.csv(fgsea_results, "C:/Users/bciez/Documents/Basilio/UCSF/output_test.csv", row.names = FALSE)

# Display the top pathways sorted by adjusted p-value (padj)
head(fgsea_results[order(fgsea_results$padj), ][, c("pathway", "NES", "pval", "padj", "size")],10)

##############################################
# Step 4: Additional Plots
##############################################

# Scatter plot: NES vs -log10(p-value) with color-coded significance
ggplot(fgsea_results, aes(x = NES, y = -log10(pval))) +
  geom_point(aes(color = pval < 0.05)) +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "NES vs -log10(p-value)", x = "Normalized Enrichment Score", y = "-log10(p-value)") +
  theme_minimal()

# Plot Pathway Size vs NES, color-coding based on p-value significance
ggplot(fgsea_results, aes(x = size, y = NES)) +
  geom_point(aes(color = pval < 0.05)) +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "Pathway Size vs NES", x = "Pathway Size", y = "Normalized Enrichment Score") +
  theme_minimal()

##############################################
# Step 5: Gene Ranks and Enrichment Analysis
##############################################

# Select the "REACTOME_INFECTIOUS_DISEASE" pathway and extract gene ranks for visualization
infectious_disease_genes <- pathways[["REACTOME_INFECTIOUS_DISEASE"]]

# Extract gene ranks from your gene list
gene_ranks <- gene_list[names(gene_list) %in% infectious_disease_genes]

# Display the ranked genes
gene_ranks

# Plot the ranked gene scores for the pathway
plot(gene_ranks, type = "h", main = "Gene Ranks for HIV Infection Pathway",
     xlab = "Gene Rank", ylab = "Ranked Score")

##############################################
# Step 6: Heatmap of Enriched Pathways
##############################################

# Extract the top 20 pathways for visualization (based on p-values)
top20_pathways <- fgsea_results[order(fgsea_results$pval), ][1:20, ]

# Create a heatmap of the NES (Normalized Enrichment Score) for the top 20 pathways
pheatmap(matrix(top20_pathways$NES, nrow = 1), 
         cluster_rows = FALSE, # Disable row clustering since it's a single row
         cluster_cols = FALSE, # Disable column clustering
         labels_col = top20_pathways$pathway, 
         main = "NES for Top 20 Enriched Pathways",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         angle_col = 45) # Rotate column names for better readability



