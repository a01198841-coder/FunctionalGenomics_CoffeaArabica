# Load ggplot2 library for data visualization and plotting
library(ggplot2)
# Load pheatmap library to generate customized heatmaps
library(pheatmap)
# Load dendsort library to reorder dendrograms in heatmaps for better visualization
library(dendsort)
# Load DESeq2 library to perform differential gene expression analysis
library(DESeq2)
# Load edgeR library (specifically to use its filterByExpr function for filtering low counts)
library(edgeR)
# Load GO.db library to access Gene Ontology database information
library(GO.db)
# Load AnnotationDbi library to query database objects like GO.db
library(AnnotationDbi)
# Load ashr library for empirical Bayes shrinkage of fold changes
library(ashr)
# Load EnhancedVolcano library to generate publication-ready volcano plots
library(EnhancedVolcano)
# Load readxl library to read data from Excel files
library(readxl)

# =========================================================================
# --- DATA IMPORT & GLOBAL GO DICTIONARY PREPARATION ---
# =========================================================================

# Read the count matrix text file, set the first column as row names, and prevent R from modifying column names
counts <- read.delim("Matrix_hisat2.txt", row.names = 1, check.names = FALSE)

# Import the gene annotation data from the specified Excel file
annot <- read_excel("GeneAnnotation.xlsx")

# Convert the imported tibble/Excel object into a standard R data frame
annot <- as.data.frame(annot)

# Extract and format the TERM2GENE mapping early for annotations
# Subset the annotation data frame to keep only the gene_id and GOs columns
gene2go <- annot[, c("gene_id", "GOs")]

# Remove rows where the GOs column is empty or contains NA values
gene2go <- gene2go[gene2go$GOs != "" & !is.na(gene2go$GOs), ]

# Create a new data frame to hold one-to-one gene-to-GO term mappings
term2gene <- data.frame(
# Unlist and split the comma-separated GO terms into individual strings for the TERM column
  TERM = unlist(strsplit(as.character(gene2go$GOs), ",\\s*")),
# Repeat the corresponding gene_id for each split GO term to populate the GENE column
  GENE = rep(gene2go$gene_id, sapply(strsplit(as.character(gene2go$GOs), ",\\s*"), length))
)

# Trim leading and trailing whitespace from the GO terms in the TERM column
term2gene$TERM <- trimws(term2gene$TERM)

# Create a TERM2NAME mapping using GO.db and AnnotationDbi
term2name <- AnnotationDbi::select(GO.db,
# Provide the unique GO IDs from our term2gene dataset as keys to search
                                   keys = unique(term2gene$TERM),
# Specify that we want to retrieve the actual term description (TERM column in GO.db)
                                   columns = "TERM",
# Tell the database that the keys we provided are GO IDs
                                   keytype = "GOID")

# Rename the columns of the output to TERM (for the ID) and NAME (for the description)
colnames(term2name) <- c("TERM", "NAME") 

# Merge the term2gene mappings with the term2name descriptions using the shared TERM column
go_merged <- merge(term2gene, term2name, by="TERM", all.x=TRUE)

# Group by GENE and collapse multiple GO IDs into a single semicolon-separated string per gene
collapsed_go <- aggregate(TERM ~ GENE, data = go_merged, FUN = function(x) paste(unique(x), collapse = "; "))

# Group by GENE and collapse multiple GO Descriptions into a single semicolon-separated string per gene
collapsed_name <- aggregate(NAME ~ GENE, data = go_merged, FUN = function(x) paste(unique(x), collapse = "; "))

# Merge the collapsed GO IDs and collapsed descriptions back together by the shared GENE column
go_dict <- merge(collapsed_go, collapsed_name, by="GENE")

# Assign clear column names to the final dictionary data frame
colnames(go_dict) <- c("LOC_Gene_ID", "GO_IDs", "GO_Descriptions")

# =========================================================================
# --- DESEQ2 OBJECT SETUP ---
# =========================================================================

# Convert the counts data frame into a numeric matrix for DESeq2
counts <- as.matrix(counts)

# Round the counts to nearest integers since DESeq2 requires unnormalized integer data
counts <- round(counts)

# Read the metadata text file, setting the first column as row names without altering column names
coldata <- read.delim("Coffee_metadata.txt", row.names = 1, check.names = FALSE)

# Correct a typo in a specific sample column name within the counts matrix (changes "Cax7" to "Ca7x")
colnames(counts)[colnames(counts) == "Cax7"] <- "Ca7x"

# Verify that all column names in the count matrix exist as row names in the metadata
all(colnames(counts) %in% rownames(coldata))

# Reorder the rows of the metadata to perfectly match the column order of the count matrix
coldata <- coldata[colnames(counts), ]

# Double-check that count columns and metadata rows are now identically aligned
all(colnames(counts) == rownames(coldata))

# Convert the cultivar column in the metadata to a factor class
coldata$cultivar <- factor(coldata$cultivar)

# Convert the treatment column in the metadata to a factor class
coldata$treatment <- factor(coldata$treatment)

# Create a new 'group' factor variable by concatenating cultivar and treatment strings
coldata$group <- factor(paste0(coldata$cultivar, "_", coldata$treatment))

# Initialize the DESeqDataSet object using the prepared data
dds <- DESeqDataSetFromMatrix(countData = counts, 
# Supply the sample metadata data frame
                              colData = coldata, 
# Define the experimental design formula to model the expression based on the 'group' factor
                              design = ~ group)

# Calculate a logical vector indicating which genes have sufficient read counts across groups to be kept
  keep <- filterByExpr(counts(dds), group = coldata$group)

# Subset the DESeqDataSet to retain only the genes that passed the filter threshold
  dds.1 <- dds[keep, ]

# Estimate the size factors for the filtered dataset to normalize for sequencing depth
  dds.1 <- estimateSizeFactors(dds.1)

# Apply a Variance Stabilizing Transformation (VST) to the data for clustering and visualization
vst.1 <- vst(dds.1)

# =========================================================================
# --- PCA ---
# =========================================================================

# Generate a PCA plot from the VST data, grouping points by cultivar and treatment
  pca_plot <- plotPCA(vst.1, intgroup = c("cultivar", "treatment")) +
# Add a clean black-and-white theme to the plot
    theme_bw() +
# Increase the size of the scatter plot points to 3 for better visibility
    geom_point(size = 3) +
# Add a main title to the plot describing the data used
    ggtitle("PCA of VST-transformed, filtered counts")
  
# Save the generated PCA plot as a PNG image file with dimensions 6x5 inches
  ggsave("PCA_Plot.png", plot = pca_plot, width = 6, height = 5)

# =========================================================================
# --- DIFFERENTIAL EXPRESSION ANALYSIS & ANNOTATED EXPORT ---
# =========================================================================

# Define the adjusted p-value threshold for statistical significance
padj_threshold <- 0.05

# Define the absolute log2 fold change threshold for biological significance
lfc_threshold <- 1

# Execute the main differential expression analysis steps (estimation of dispersions and model fitting)
dds.1 <- DESeq(dds.1)

# Print the levels of the group factor to the console to verify available contrast groups
print(levels(dds.1$group)) 

# Extract standard results for the contrast between Catuai xylella and Catuai saline groups
res_Cultivar1 <- results(dds.1, contrast = c("group", "Catuai_xylella", "Catuai_saline"))

# Apply shrinkage to the log2 fold changes for the Catuai contrast to reduce noise in low-count genes
res_Cultivar1_shrunk <- lfcShrink(dds.1, contrast = c("group", "Catuai_xylella", "Catuai_saline"), 
# Specify the use of the adaptive shrinkage ('ashr') estimator
                                  type = "ashr")

# Subset the Catuai shrunk results to isolate genes passing both padj and log2FoldChange thresholds
sig_Cultivar1 <- subset(res_Cultivar1_shrunk, padj < padj_threshold & abs(log2FoldChange) >= lfc_threshold)

# Print the total count of significant DEGs identified in the Catuai cultivar to the console
print(paste("Significant DEGs in Catuai:", nrow(sig_Cultivar1)))

# Convert the significant Catuai DESeqResults object into a standard data frame
catuai_df <- as.data.frame(sig_Cultivar1)

# Extract the row names (gene IDs) and place them into a new column called LOC_Gene_ID
catuai_df$LOC_Gene_ID <- rownames(catuai_df)

# Merge the Catuai significant genes with the GO dictionary to add functional annotations
catuai_annotated <- merge(catuai_df, go_dict, by = "LOC_Gene_ID", all.x = TRUE)

# Export the fully annotated significant Catuai genes to a CSV file, embedding the LFC threshold in the name
write.csv(catuai_annotated, file=paste0("Catuai_Xylella_vs_Control_DEGs_LFC", lfc_threshold, ".csv"), row.names=FALSE)

# Extract standard results for the contrast between CR95 xylella and CR95 saline groups
res_Cultivar2 <- results(dds.1, contrast = c("group", "CR95_xylella", "CR95_saline"))

# Apply shrinkage to the log2 fold changes for the CR95 contrast to reduce noise in low-count genes
res_Cultivar2_shrunk <- lfcShrink(dds.1, contrast = c("group", "CR95_xylella", "CR95_saline"), 
# Specify the use of the adaptive shrinkage ('ashr') estimator
                                  type = "ashr")

# Subset the CR95 shrunk results to isolate genes passing both padj and log2FoldChange thresholds
sig_Cultivar2 <- subset(res_Cultivar2_shrunk, padj < padj_threshold & abs(log2FoldChange) >= lfc_threshold)

# Print the total count of significant DEGs identified in the CR95 cultivar to the console
print(paste("Significant DEGs in CR95:", nrow(sig_Cultivar2)))

# Convert the significant CR95 DESeqResults object into a standard data frame
cr95_df <- as.data.frame(sig_Cultivar2)

# Extract the row names (gene IDs) and place them into a new column called LOC_Gene_ID
cr95_df$LOC_Gene_ID <- rownames(cr95_df)

# Merge the CR95 significant genes with the GO dictionary to add functional annotations
cr95_annotated <- merge(cr95_df, go_dict, by = "LOC_Gene_ID", all.x = TRUE)

# Export the fully annotated significant CR95 genes to a CSV file, embedding the LFC threshold in the name
write.csv(cr95_annotated, file=paste0("CR95_Xylella_vs_Control_DEGs_LFC", lfc_threshold, ".csv"), row.names=FALSE)

# =========================================================================
# --- VOLCANO PLOTS ---
# =========================================================================

# ---------------------------------------------------------
# Catuai Volcano Plot Data Prep
# ---------------------------------------------------------
# Extract the list of all analyzed gene IDs from the Catuai shrunk results
catuai_genes <- rownames(res_Cultivar1_shrunk)

# Look up the preferred names for these genes by matching them against the main annotation file
catuai_match <- annot$Preferred_name[match(catuai_genes, annot$gene_id)]

# Fallback to gene_id if Preferred_name is NA, empty, or "0"
catuai_labels <- ifelse(is.na(catuai_match) | catuai_match == "" | catuai_match == "0", 
# Use the original gene ID as the fallback label
                        catuai_genes, 
# Otherwise, use the valid preferred name
                        catuai_match)

# Define custom colors for Up/Down/NS
keyvals_catuai <- ifelse(
# If the adjusted p-value is NA, assign grey50
  is.na(res_Cultivar1_shrunk$padj), 'grey50',
# If padj is significant and LFC is highly positive, assign red (#E41A1C)
  ifelse(res_Cultivar1_shrunk$padj < padj_threshold & res_Cultivar1_shrunk$log2FoldChange >= lfc_threshold, '#E41A1C', # Upregulated (Red)
# If padj is significant and LFC is highly negative, assign blue (#377EB8)
         ifelse(res_Cultivar1_shrunk$padj < padj_threshold & res_Cultivar1_shrunk$log2FoldChange <= -lfc_threshold, '#377EB8', # Downregulated (Blue)
# If none of the conditions apply (not significant), assign grey50
                'grey50'))) # Not Significant (Grey)

# Assign the string label 'Upregulated' to the red hex color in the named vector
names(keyvals_catuai)[keyvals_catuai == '#E41A1C'] <- 'Upregulated'

# Assign the string label 'Downregulated' to the blue hex color in the named vector
names(keyvals_catuai)[keyvals_catuai == '#377EB8'] <- 'Downregulated'

# Assign the string label 'Not Significant' to the grey color in the named vector
names(keyvals_catuai)[keyvals_catuai == 'grey50'] <- 'Not Significant'

# Plot Catuai
vol_catuai <- EnhancedVolcano(res_Cultivar1_shrunk,
# Pass the custom labels array prepared earlier
                              lab = catuai_labels,
# Specify the column name representing the x-axis (fold change)
                              x = 'log2FoldChange',
# Specify the column name representing the y-axis (significance)
                              y = 'padj',
# Define the main title for the volcano plot
                              title = 'Catuai: Xylella vs Control',
# Define the subtitle dynamically displaying the active thresholds
                              subtitle = paste0('Thresholds: p-adj < ', padj_threshold, ', |Log2FC| >= ', lfc_threshold),
# Instruct the plot where to draw the horizontal p-value cutoff line
                              pCutoff = padj_threshold,
# Instruct the plot where to draw the vertical fold change cutoff lines
                              FCcutoff = lfc_threshold,
# Set the size of the individual scatter points
                              pointSize = 2.0,
# Set the size of the text labels applied to top genes
                              labSize = 3.0,
# Pass the custom color keyvals vector mapped to the genes
                              colCustom = keyvals_catuai,
# Set the transparency level of the points
                              colAlpha = 0.7,
# Position the legend on the right side of the plot
                              legendPosition = 'right'
)

# Save the configured Catuai volcano plot to a high-resolution PNG file
ggsave(paste0("Volcano_Catuai_LFC", lfc_threshold, ".png"), plot = vol_catuai, width = 10, height = 8, dpi = 300)

# ---------------------------------------------------------
# CR95 Volcano Plot Data Prep
# ---------------------------------------------------------
# Extract the list of all analyzed gene IDs from the CR95 shrunk results
cr95_genes <- rownames(res_Cultivar2_shrunk)

# Look up the preferred names for these genes by matching them against the main annotation file
cr95_match <- annot$Preferred_name[match(cr95_genes, annot$gene_id)]

# Fallback to gene_id if Preferred_name is NA, empty, or "0"
cr95_labels <- ifelse(is.na(cr95_match) | cr95_match == "" | cr95_match == "0", 
# Use the original gene ID as the fallback label
                      cr95_genes, 
# Otherwise, use the valid preferred name
                      cr95_match)

# Define custom colors for Up/Down/NS
keyvals_cr95 <- ifelse(
# If the adjusted p-value is NA, assign grey50
  is.na(res_Cultivar2_shrunk$padj), 'grey50',
# If padj is significant and LFC is highly positive, assign red (#E41A1C)
  ifelse(res_Cultivar2_shrunk$padj < padj_threshold & res_Cultivar2_shrunk$log2FoldChange >= lfc_threshold, '#E41A1C', # Upregulated (Red)
# If padj is significant and LFC is highly negative, assign blue (#377EB8)
         ifelse(res_Cultivar2_shrunk$padj < padj_threshold & res_Cultivar2_shrunk$log2FoldChange <= -lfc_threshold, '#377EB8', # Downregulated (Blue)
# If none of the conditions apply (not significant), assign grey50
                'grey50')))

# Assign the string label 'Upregulated' to the red hex color in the named vector
names(keyvals_cr95)[keyvals_cr95 == '#E41A1C'] <- 'Upregulated'

# Assign the string label 'Downregulated' to the blue hex color in the named vector
names(keyvals_cr95)[keyvals_cr95 == '#377EB8'] <- 'Downregulated'

# Assign the string label 'Not Significant' to the grey color in the named vector
names(keyvals_cr95)[keyvals_cr95 == 'grey50'] <- 'Not Significant'

# Plot CR95
vol_cr95 <- EnhancedVolcano(res_Cultivar2_shrunk,
# Pass the custom labels array prepared earlier
                            lab = cr95_labels,
# Specify the column name representing the x-axis (fold change)
                            x = 'log2FoldChange',
# Specify the column name representing the y-axis (significance)
                            y = 'padj',
# Define the main title for the volcano plot
                            title = 'CR95: Xylella vs Control',
# Define the subtitle dynamically displaying the active thresholds
                            subtitle = paste0('Thresholds: p-adj < ', padj_threshold, ', |Log2FC| >= ', lfc_threshold),
# Instruct the plot where to draw the horizontal p-value cutoff line
                            pCutoff = padj_threshold,
# Instruct the plot where to draw the vertical fold change cutoff lines
                            FCcutoff = lfc_threshold,
# Set the size of the individual scatter points
                            pointSize = 2.0,
# Set the size of the text labels applied to top genes
                            labSize = 3.0,
# Pass the custom color keyvals vector mapped to the genes
                            colCustom = keyvals_cr95,
# Set the transparency level of the points
                            colAlpha = 0.7,
# Position the legend on the right side of the plot
                            legendPosition = 'right'
)

# Save the configured CR95 volcano plot to a high-resolution PNG file
ggsave(paste0("Volcano_CR95_LFC", lfc_threshold, ".png"), plot = vol_cr95, width = 10, height = 8, dpi = 300)

# =========================================================================
# --- UPREGULATED VS DOWNREGULATED BARPLOT ---
# =========================================================================

# Count the exact number of significantly upregulated genes in Catuai using logical conditions
catuai_up <- sum(res_Cultivar1_shrunk$padj < padj_threshold & res_Cultivar1_shrunk$log2FoldChange >= lfc_threshold, na.rm = TRUE)

# Count the exact number of significantly downregulated genes in Catuai using logical conditions
catuai_down <- sum(res_Cultivar1_shrunk$padj < padj_threshold & res_Cultivar1_shrunk$log2FoldChange <= -lfc_threshold, na.rm = TRUE)

# Count the exact number of significantly upregulated genes in CR95 using logical conditions
cr95_up <- sum(res_Cultivar2_shrunk$padj < padj_threshold & res_Cultivar2_shrunk$log2FoldChange >= lfc_threshold, na.rm = TRUE)

# Count the exact number of significantly downregulated genes in CR95 using logical conditions
cr95_down <- sum(res_Cultivar2_shrunk$padj < padj_threshold & res_Cultivar2_shrunk$log2FoldChange <= -lfc_threshold, na.rm = TRUE)

# Construct a summary data frame combining the counts for plotting
deg_counts <- data.frame(
# Create a column designating the cultivar for each count metric
  Cultivar = c("Catuai", "Catuai", "CR95", "CR95"),
# Create a column designating the direction of regulation for each count metric
  Direction = c("Upregulated", "Downregulated", "Upregulated", "Downregulated"),
# Populate the counts derived from the summation logic above
  Count = c(catuai_up, catuai_down, cr95_up, cr95_down)
)

# Convert the Direction column to a factor with explicitly ordered levels for consistent plotting order
deg_counts$Direction <- factor(deg_counts$Direction, levels = c("Upregulated", "Downregulated"))

# Initialize a ggplot bar chart mapped to the summary data frame variables
p_bar <- ggplot(deg_counts, aes(x = Cultivar, y = Count, fill = Direction)) +
# Add the bar geometry, positioning bars side-by-side (dodge) and giving them a black border
  geom_bar(stat = "identity", position = "dodge", color = "black") +
# Add text labels displaying the specific count number directly over each respective bar
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5, size = 4) +
# Manually apply distinct colors to the bar fills based on regulation direction
  scale_fill_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8")) +
# Apply a clean black-and-white visual theme to the plot
  theme_bw() +
# Set the text details for the plot
  labs(
# Assign the main title for the bar plot
    title = "Differentially Expressed Genes by Cultivar",
# Assign a subtitle documenting the significance thresholds used
    subtitle = paste0("Filters: padj < ", padj_threshold, ", |log2FC| >= ", lfc_threshold),
# Assign the x-axis label
    x = "Cultivar",
# Assign the y-axis label
    y = "Number of Genes"
  ) +
# Expand the y-axis upper limit slightly to ensure the text labels don't get cut off
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) 

# Save the bar chart plot to a high-resolution PNG file
ggsave(paste0("DEG_Counts_Barplot_LFC", lfc_threshold, ".png"), plot = p_bar, width = 6, height = 5, dpi = 300)

# =========================================================================
# --- HEATMAP OF SIGNIFICANT DEGs ---
# =========================================================================

# Extract the normalized data matrix from the VST transformed DESeq object
vst_assay <- assay(vst.1)

# Combine the significant gene IDs from both Catuai and CR95 into a single unique vector
all_sig_genes <- unique(c(genes_Cultivar1, genes_Cultivar2))

# Subset the VST matrix to retain only the rows corresponding to those significant genes
heatmap_matrix <- vst_assay[all_sig_genes, ]

# Create a column annotation data frame from the metadata for mapping to the heatmap header
anno_col <- as.data.frame(coldata[, c("cultivar", "treatment")])

# Perform hierarchical clustering on the rows (genes) and apply dendsort to organize the dendrogram structure
sorted_row_cluster <- dendsort(hclust(dist(heatmap_matrix)))

# Perform hierarchical clustering on the columns (samples) and apply dendsort to organize the dendrogram structure
sorted_col_cluster <- dendsort(hclust(dist(t(heatmap_matrix))))

# Generate a complex heatmap using the subsetted matrix
pheatmap(
# Supply the matrix containing normalized values for significant genes
  mat = heatmap_matrix,
# Scale values by row (gene) to highlight relative expression changes across samples rather than absolute values
  scale = "row",                     
# Attach the metadata annotation block to visually label the sample columns
  annotation_col = anno_col,         
# Provide the custom-sorted row dendrogram to organize genes
  cluster_rows = sorted_row_cluster, 
# Provide the custom-sorted column dendrogram to organize samples
  cluster_cols = sorted_col_cluster, 
# Disable the display of row names (gene IDs) to prevent visual clutter
  show_rownames = FALSE,             
# Enable the display of column names (sample names)
  show_colnames = TRUE,              
# Define a 50-step color gradient palette going from navy blue (low) to firebrick red (high)
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
# Save the resulting heatmap directly to a PNG file
  filename = "Heatmap_AllSigDEGs.png",
# Set the width of the generated heatmap image
  width = 8, 
# Set the height of the generated heatmap image
  height = 10
)

# =========================================================================
# --- KEGG MAPPER COLOR EXPORT ---
# =========================================================================

# Subset the primary annotation data frame to keep only the gene_id and KEGG_ko columns
kegg_mapping <- annot[, c("gene_id", "KEGG_ko")]

# Trim leading and trailing whitespace from the KO term entries
kegg_mapping$KEGG_ko <- trimws(kegg_mapping$KEGG_ko)

# Filter out any rows containing invalid, empty, missing, or placeholder KO mappings
kegg_mapping <- kegg_mapping[kegg_mapping$KEGG_ko != "" & 
# Exclude explicit NA values
                               !is.na(kegg_mapping$KEGG_ko) & 
# Exclude hyphen placeholders
                               kegg_mapping$KEGG_ko != "-" & 
# Exclude literal string representations of NA
                               tolower(kegg_mapping$KEGG_ko) != "na", ]

# -------------------------------------------------------------------------
# --- CATUAI EXPORT
# -------------------------------------------------------------------------
# Convert the significant Catuai DEGs object to a standard data frame
df_catuai <- as.data.frame(sig_Cultivar1)

# Explicitly create a gene_id column from the row names to permit merging
df_catuai$gene_id <- rownames(df_catuai)

# Perform an inner merge to join the significant Catuai genes with their corresponding KEGG KO terms
df_catuai_kegg <- merge(df_catuai, kegg_mapping, by = "gene_id")

# Create a new column assigning 'red' to upregulated genes and 'blue' to downregulated ones
df_catuai_kegg$color <- ifelse(df_catuai_kegg$log2FoldChange > 0, "red", "blue")

# Aggregate colors by KO term to check for mixed regulation
kegg_catuai_agg <- aggregate(color ~ KEGG_ko, data = df_catuai_kegg, FUN = function(x) {
# Determine the unique colors present for this specific KO term
  unique_colors <- unique(x)
# Evaluate if more than one color (both red and blue) maps to the exact same KO term
  if (length(unique_colors) > 1) {
# Return 'yellow' to represent a KO pathway containing both up and down-regulated genes
    return("yellow") # Represents a KO containing both up and down-regulated genes
# If all mapped genes move in the same direction...
  } else {
# Return that single specific color (red or blue)
    return(unique_colors)
  }
})

# Export the formatted KEGG data to a text file compatible with the KEGG Mapper Color web tool
write.table(kegg_catuai_agg, file = paste0("Catuai_KEGG_Upload_LFC", lfc_threshold, ".txt"), 
# Enforce output formatting by omitting row names, headers, and quote marks while using tabs as separators
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
# Print a confirmation message indicating successful export of the Catuai KEGG data
print("Exported clean Catuai dataset for KEGG Mapper Color.")

# -------------------------------------------------------------------------
# --- CR95 EXPORT
# -------------------------------------------------------------------------
# Convert the significant CR95 DEGs object to a standard data frame
df_cr95 <- as.data.frame(sig_Cultivar2)

# Explicitly create a gene_id column from the row names to permit merging
df_cr95$gene_id <- rownames(df_cr95)

# Perform an inner merge to join the significant CR95 genes with their corresponding KEGG KO terms
df_cr95_kegg <- merge(df_cr95, kegg_mapping, by = "gene_id")

# Create a new column assigning 'red' to upregulated genes and 'blue' to downregulated ones
df_cr95_kegg$color <- ifelse(df_cr95_kegg$log2FoldChange > 0, "red", "blue")

# Aggregate colors by KO term to check for mixed regulation
kegg_cr95_agg <- aggregate(color ~ KEGG_ko, data = df_cr95_kegg, FUN = function(x) {
# Determine the unique colors present for this specific KO term
  unique_colors <- unique(x)
# Evaluate if more than one color (both red and blue) maps to the exact same KO term
  if (length(unique_colors) > 1) {
# Return 'yellow' to represent a KO pathway containing both up and down-regulated genes
    return("yellow") # Represents a KO containing both up and down-regulated genes
# If all mapped genes move in the same direction...
  } else {
# Return that single specific color (red or blue)
    return(unique_colors)
  }
})

# Export the formatted KEGG data to a text file compatible with the KEGG Mapper Color web tool
write.table(kegg_cr95_agg, file = paste0("CR95_KEGG_Upload_LFC", lfc_threshold, ".txt"), 
# Enforce output formatting by omitting row names, headers, and quote marks while using tabs as separators
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
# Print a confirmation message indicating successful export of the CR95 KEGG data
print("Exported clean CR95 dataset for KEGG Mapper Color.")
