# Load edgeR for RNA-seq count filtering, normalization, and DGEList construction
library(edgeR)

# Load limma for voom transformation, linear modeling, gene-set analysis, and CAMERA
library(limma)

# Load ggplot2 for visualization of enriched GO terms
library(ggplot2)

# Load GO.db to retrieve Gene Ontology term names and ontology classifications
library(GO.db)

# Load AnnotationDbi for querying and retrieving GO annotations
library(AnnotationDbi)

# Load igraph for constructing and manipulating the GO-term similarity network
library(igraph)

# Load ggraph for visualizing the GO-term enrichment network
library(ggraph)


# =========================================================================
# --- DATA IMPORT & PREPARATION ---
# =========================================================================

# Import the RNA-seq count matrix; rows correspond to genes and columns to samples
counts <- read.table("Matrix_hisat2.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Import the experimental metadata containing sample, cultivar, and treatment information
metadata <- read.table("Coffee_metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Correct a sample-name discrepancy so that the count matrix and metadata use identical sample identifiers
colnames(counts)[colnames(counts) == "Cax7"] <- "Ca7x"

# Reorder the count matrix columns to exactly match the sample order specified in the metadata
counts <- counts[, match(metadata$sample, colnames(counts))]

# Round counts to integer values because RNA-seq count-based methods require integer-like read counts
counts <- round(counts)


# NOTE: Using the TSV here as in your original script.
# Import the gene annotation table containing Gene Ontology (GO) assignments for each gene
annotation <- read.table("Coffee_fullAnnotation.tsv.txt", header = TRUE, sep = "\t",
                         quote = "\"", comment.char = "", fill = TRUE,
                         stringsAsFactors = FALSE, na.strings = c("-", "NA", ""))

# Identify genes for which at least one GO annotation is available
has_go <- !is.na(annotation$GOs)

# Separate comma-delimited GO identifiers so each gene can be associated with individual GO terms
gene2go <- strsplit(annotation$GOs[has_go], ",")

# Assign each gene identifier as the name of its corresponding GO-term vector
names(gene2go) <- annotation$gene_id[has_go]

# Convert the gene-to-GO relationships into a GO-term-to-gene mapping.
# This structure is required by limma::ids2indices for gene-set analysis
term2gene <- split(rep(names(gene2go), lengths(gene2go)), unlist(gene2go))

# Define the minimum number of genes required for a GO term to be analyzed
# This reduces the influence of very small gene sets on the enrichment statistics
min_set_size <- 5


# =========================================================================
# --- CORE FUNCTION: CAMERA WORKFLOW ---
# =========================================================================

# Define a reusable function that performs the complete CAMERA analysis
# independently for each cultivar
run_camera_analysis <- function(cultivar_name, counts, metadata, term2gene, min_size) {
  
  # -----------------------------------------------------------------------
  # 1. SUBSET DATA
  # -----------------------------------------------------------------------
  
  # Select metadata corresponding only to the cultivar being analyzed
  meta_sub <- metadata[metadata$cultivar == cultivar_name, ]
  
  # Select the corresponding RNA-seq count data for that cultivar
  counts_sub <- counts[, meta_sub$sample]
  
  # Define treatment as a categorical variable with saline as the reference level
  # and Xylella as the comparison treatment
  treat <- factor(meta_sub$treatment, levels = c("saline", "xylella"))
  
  
  # -----------------------------------------------------------------------
  # 2. DGEList & NORMALIZATION
  # -----------------------------------------------------------------------
  
  # Create an edgeR DGEList object containing counts and treatment information
  # This object provides the framework for RNA-seq library-size normalization and filtering
  dge <- DGEList(counts = counts_sub, group = treat)
  
  # Identify genes with sufficient expression to provide reliable statistical information
  # filterByExpr uses the experimental design/group structure to remove very lowly expressed genes
  keep <- filterByExpr(dge)
  
  # Retain only genes passing the expression filter
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Calculate normalization factors to account for differences in sequencing depth
  # and RNA composition between libraries
  dge <- calcNormFactors(dge)
  
  
  # -----------------------------------------------------------------------
  # 3. VOOM TRANSFORMATION
  # -----------------------------------------------------------------------
  
  # Construct the linear-model design matrix representing the treatment comparison
  # The intercept corresponds to saline and the treatment coefficient represents Xylella
  design <- model.matrix(~treat)
  
  # Apply limma::voom to transform count data to log2-expression values
  # while estimating the mean-variance relationship of the RNA-seq data.
  # Voom assigns precision weights to observations, allowing heteroscedastic
  # RNA-seq expression data to be analyzed using linear models
  v <- voom(dge, design, plot = FALSE)
  
  
  # -----------------------------------------------------------------------
  # 4. MAP GENES TO GENE-SET INDICES
  # -----------------------------------------------------------------------
  
  # Convert GO term-to-gene relationships into gene indices corresponding
  # to the genes retained after expression filtering
  idx <- ids2indices(term2gene, rownames(v))
  
  # Retain only GO terms containing at least the specified minimum number of genes
  # to avoid unstable enrichment estimates from very small gene sets
  idx <- idx[lengths(idx) >= min_size]
  
  
  # -----------------------------------------------------------------------
  # 5. RUN CAMERA
  # -----------------------------------------------------------------------
  
  # Perform CAMERA competitive gene-set testing.
  # CAMERA evaluates whether genes belonging to each GO set show coordinated
  # expression changes relative to genes outside the set.
  # Importantly, CAMERA accounts for correlation between genes within each set,
  # which prevents treating correlated genes as statistically independent.
  # contrast = 2 tests the treatment coefficient, corresponding to Xylella vs saline
  res <- camera(v, index = idx, design = design, contrast = 2)
  
  # Store the GO identifier as an explicit column for downstream annotation and comparison
  res$GOID <- rownames(res)
  
  
  # -----------------------------------------------------------------------
  # 6. FETCH GO ANNOTATIONS
  # -----------------------------------------------------------------------
  
  # Retrieve the biological description (TERM) and ontology category (BP, CC, or MF)
  # corresponding to each GO identifier using the GO.db annotation database
  term_info <- suppressMessages(
    select(GO.db, keys = res$GOID, columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
  )
  
  # Combine CAMERA statistics with their corresponding GO annotations
  res <- merge(res, term_info, by = "GOID", all.x = TRUE)
  
  # Use the GO identifier when a descriptive GO term is unavailable
  res$TERM[is.na(res$TERM)] <- res$GOID[is.na(res$TERM)]
  
  # Assign "Unknown" when an ontology classification is unavailable
  res$ONTOLOGY[is.na(res$ONTOLOGY)] <- "Unknown"
  
  # Sort GO terms according to their unadjusted CAMERA P-value
  res <- res[order(res$PValue), ]
  
  
  # -----------------------------------------------------------------------
  # 7. EXPORT FULL RESULTS
  # -----------------------------------------------------------------------
  
  # Export the complete CAMERA result table, including enrichment statistics,
  # P-values, FDR values, GO terms, and ontology information
  write.table(res, file = paste0("CAMERA_", cultivar_name, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  
  # -----------------------------------------------------------------------
  # 8. FILTER AND PLOT SIGNIFICANT TERMS
  # -----------------------------------------------------------------------
  
  # Select GO terms with a false discovery rate below 5%.
  # FDR controls for multiple testing across the large number of GO terms evaluated
  sig <- res[res$FDR < 0.05, ]
  
  # Initialize an empty object for the top significant GO terms
  top <- data.frame()
  
  
  # Only generate a plot when at least one GO term passes the FDR threshold
  if (nrow(sig) > 0) {
    
    # Select the 20 GO terms with the smallest FDR values
    top <- head(sig[order(sig$FDR), ], 20)
    
    # Reverse factor order so the most significant terms appear at the top
    # of the horizontal bar plot
    top$TERM <- factor(top$TERM, levels = rev(unique(top$TERM)))
    
    
    # Generate a bar plot of the 20 most significant GO terms.
    # -log10(FDR) transforms small FDR values into larger positive values,
    # making statistical significance easier to visualize
    p_bar <- ggplot(top, aes(x = TERM, y = -log10(FDR), fill = Direction)) +
      
      # Draw one bar for each enriched GO term
      geom_col() +
      
      # Rotate the coordinate system so GO terms are displayed horizontally
      coord_flip() +
      
      # Assign colors according to the direction of enrichment
      scale_fill_manual(values = c(Up = "firebrick", Down = "steelblue")) +
      
      # Define plot title and axis labels
      labs(title = cultivar_name, x = "GO term", y = "-log10(FDR)") +
      
      # Apply a clean black-and-white plotting theme
      theme_bw() +
      
      # Reduce GO-term label size to improve readability
      theme(axis.text.y = element_text(size = 6))
    
    
    # Save the GO enrichment bar plot as a PNG image
    ggsave(paste0("barplot_", tolower(cultivar_name), "_hisat2_camera.png"), p_bar, width = 8, height = 6)
    
  } else {
    
    # Report when no GO term meets the predefined FDR significance threshold
    message("No significant GO terms (FDR < 0.05) found for ", cultivar_name)
  }
  
  
  # Attach cultivar information to facilitate comparisons between analyses
  res$cultivar <- cultivar_name
  sig$cultivar <- cultivar_name
  
  # Attach cultivar information to the top-term table when significant terms are present
  if(nrow(top) > 0) top$cultivar <- cultivar_name
  
  # Return complete results, significant results, and the top 20 terms
  return(list(res = res, sig = sig, top = top))
}


# =========================================================================
# --- EXECUTE ANALYSIS ---
# =========================================================================

# Run the complete CAMERA workflow independently for Catuai
results_Catuai <- run_camera_analysis("Catuai", counts, metadata, term2gene, min_set_size)

# Run the complete CAMERA workflow independently for CR95
results_CR95 <- run_camera_analysis("CR95", counts, metadata, term2gene, min_set_size)


# =========================================================================
# --- CROSS-CULTIVAR COMPARISONS (ENRICHMENT NETWORK MAP) ---
# =========================================================================

# Extract the top 20 significant GO terms identified for each cultivar
top_Catuai <- results_Catuai$top
top_CR95 <- results_CR95$top

# Create a non-redundant list containing GO terms significant in either cultivar
# These terms constitute the nodes considered for the cross-cultivar network
comparison_terms <- unique(c(top_Catuai$GOID, top_CR95$GOID))


# Proceed only when at least one significant GO term is available for comparison
if (length(comparison_terms) > 0) {
  
  # -----------------------------------------------------------------------
  # 1. EXTRACT GENES FOR THE TOP TERMS
  # -----------------------------------------------------------------------
  
  # Retrieve the genes belonging to each GO term selected for comparison
  term_genes <- term2gene[comparison_terms]
  
  # Remove GO terms for which no genes are available
  term_genes <- term_genes[!sapply(term_genes, is.null)]
  
  # Store the GO identifiers of terms with valid gene mappings
  valid_terms <- names(term_genes)
  
  
  # Proceed only when at least two GO terms are available for comparison
  if (length(valid_terms) > 1) {
    
    # ---------------------------------------------------------------------
    # 2. CALCULATE JACCARD SIMILARITY MATRIX
    # ---------------------------------------------------------------------
    
    # Initialize a matrix to store pairwise gene-set similarity values
    sim_matrix <- matrix(0, nrow = length(valid_terms), ncol = length(valid_terms),
                         dimnames = list(valid_terms, valid_terms))
    
    
    # Calculate pairwise Jaccard similarity between all GO-term gene sets
    for(i in 1:(length(valid_terms)-1)) {
      for(j in (i+1):length(valid_terms)) {
        
        # Retrieve the genes associated with each pair of GO terms
        g1 <- term_genes[[i]]
        g2 <- term_genes[[j]]
        
        # Calculate the Jaccard Index:
        # number of genes shared by both terms / total number of unique genes
        # represented across the two terms
        jaccard <- length(intersect(g1, g2)) / length(union(g1, g2))
        
        # Store the similarity value for the pair of GO terms
        sim_matrix[i, j] <- jaccard
        sim_matrix[j, i] <- jaccard
      }
    }
    
    
    # ---------------------------------------------------------------------
    # 3. CREATE NETWORK GRAPH
    # ---------------------------------------------------------------------
    
    # Define the minimum Jaccard similarity required to connect two GO terms
    jaccard_threshold <- 0.2
    
    # Remove connections representing less than 20% gene overlap
    # This simplifies the network by retaining only relatively similar gene sets
    sim_matrix[sim_matrix < jaccard_threshold] <- 0
    
    # Convert the similarity matrix into an undirected weighted network.
    # GO terms become nodes and significant gene-set overlaps become edges.
    g <- graph_from_adjacency_matrix(sim_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
    
    
    # ---------------------------------------------------------------------
    # 4. PREPARE NODE DATA
    # ---------------------------------------------------------------------
    
    # Create a data frame containing the GO identifiers represented by network nodes
    node_data <- data.frame(GOID = V(g)$name, stringsAsFactors = FALSE)
    
    # Identify GO terms present among the top significant terms for Catuai
    node_data$in_Catuai <- node_data$GOID %in% top_Catuai$GOID
    
    # Identify GO terms present among the top significant terms for CR95
    node_data$in_CR95 <- node_data$GOID %in% top_CR95$GOID
    
    # Classify each node according to the cultivar(s) in which it occurs
    node_data$Group <- ifelse(node_data$in_Catuai & node_data$in_CR95, "Both",
                              ifelse(node_data$in_Catuai, "Catuai", "CR95"))
    
    
    # Retrieve CAMERA results for the GO terms represented in the network
    comparison_Catuai <- results_Catuai$res[results_Catuai$res$GOID %in% valid_terms, ]
    comparison_CR95 <- results_CR95$res[results_CR95$res$GOID %in% valid_terms, ]
    
    # Combine the results from both cultivars
    comparison_df <- rbind(comparison_Catuai, comparison_CR95)
    
    # Retain GO identifiers, term descriptions, and the number of genes per term
    term_info <- comparison_df[!duplicated(comparison_df$GOID), c("GOID", "TERM", "NGenes")]
    
    # Add GO-term descriptions and gene counts to the network-node information
    node_data <- merge(node_data, term_info, by = "GOID", sort = FALSE)
    
    # Restore the node order to match the order of vertices in the igraph object
    node_data <- node_data[match(V(g)$name, node_data$GOID), ]
    
    # Store GO-term descriptions as vertex attributes for network visualization
    V(g)$TERM <- node_data$TERM
    
    # Store cultivar classification as a vertex attribute
    V(g)$Group <- node_data$Group
    
    # Store the number of genes associated with each GO term as a vertex attribute
    V(g)$Size <- node_data$NGenes
    
    
    # ---------------------------------------------------------------------
    # 5. PLOT THE ENRICHMENT MAP
    # ---------------------------------------------------------------------
    
    # Set the random seed so that the force-directed network layout
    # is reproducible across different executions of the script
    set.seed(42)
    
    # Construct the enrichment network using ggraph and a Fruchterman-Reingold layout
    # in which connected GO terms are positioned according to their network relationships
    p_emap <- ggraph(g, layout = "fr") +
      
      # Draw edges representing gene-set overlap.
      # Greater edge width corresponds to greater Jaccard similarity
      geom_edge_link(aes(edge_width = weight), alpha = 0.3, color = "darkgray") +
      
      # Draw GO terms as nodes.
      # Node color identifies cultivar association and node size represents
      # the number of genes belonging to the GO term
      geom_node_point(aes(color = Group, size = Size), alpha = 0.8) +
      
      # Add GO-term names to network nodes while attempting to minimize label overlap
      geom_node_text(aes(label = TERM), repel = TRUE, size = 3, max.overlaps = 20) +
      
      # Define the visual range for node sizes
      scale_size_continuous(range = c(3, 10)) +
      
      # Define the visual range for edge widths
      scale_edge_width_continuous(range = c(0.5, 2)) +
      
      # Define colors corresponding to cultivar-specific or shared enrichment
      scale_color_manual(values = c("Catuai" = "firebrick", "CR95" = "steelblue", "Both" = "purple")) +
      
      # Add titles and legends explaining the biological and statistical meaning
      # of the network elements
      labs(title = "Enrichment Network Map",
           subtitle = paste("Edges represent >", jaccard_threshold*100, "% gene overlap (Jaccard Index)"),
           color = "Cultivar", 
           size = "Genes in Term", 
           edge_width = "Overlap") +
      
      # Remove conventional plot axes because the network topology itself
      # represents the relationships between GO terms
      theme_void() +
      
      # Position the legend on the right side of the figure
      theme(
        legend.position = "right",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
    
    
    # Save the cross-cultivar enrichment network as a PNG image
    ggsave("emap_compare_hisat2_camera.png", p_emap, width = 10, height = 8, bg = "white")
    
  } else {
    
    # Report when fewer than two valid GO terms are available for network construction
    message("Not enough significant terms with mapped genes to generate a network map.")
  }
}
