library(edgeR)
library(limma)
library(ggplot2)
library(GO.db)
library(AnnotationDbi)

# =========================================================================
# --- DATA IMPORT & PREPARATION ---
# =========================================================================

counts <- read.table("Matrix_hisat2.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
metadata <- read.table("Coffee_metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

colnames(counts)[colnames(counts) == "Cax7"] <- "Ca7x"
counts <- counts[, match(metadata$sample, colnames(counts))]
counts <- round(counts)

# NOTE: Using the TSV here as in your original script. 
# If you want to match the differential expression script perfectly, 
# consider switching this to readxl::read_excel("GeneAnnotation.xlsx")
annotation <- read.table("Coffee_fullAnnotation.tsv.txt", header = TRUE, sep = "\t",
                         quote = "\"", comment.char = "", fill = TRUE,
                         stringsAsFactors = FALSE, na.strings = c("-", "NA", ""))

has_go <- !is.na(annotation$GOs)
gene2go <- strsplit(annotation$GOs[has_go], ",")
names(gene2go) <- annotation$gene_id[has_go]
term2gene <- split(rep(names(gene2go), lengths(gene2go)), unlist(gene2go))

min_set_size <- 5

# =========================================================================
# --- CORE FUNCTION: CAMERA WORKFLOW ---
# =========================================================================

run_camera_analysis <- function(cultivar_name, counts, metadata, term2gene, min_size) {
  
  # 1. Subset Data
  meta_sub <- metadata[metadata$cultivar == cultivar_name, ]
  counts_sub <- counts[, meta_sub$sample]
  treat <- factor(meta_sub$treatment, levels = c("saline", "xylella"))
  
  # 2. DGEList & Normalization
  dge <- DGEList(counts = counts_sub, group = treat)
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  # 3. Voom Transformation
  design <- model.matrix(~treat)
  v <- voom(dge, design, plot = FALSE)
  
  # 4. Map genes to indices efficiently using limma::ids2indices
  idx <- ids2indices(term2gene, rownames(v))
  idx <- idx[lengths(idx) >= min_size]
  
  # 5. Run CAMERA
  res <- camera(v, index = idx, design = design, contrast = 2)
  res$GOID <- rownames(res)
  
  # 6. Fetch Annotations
  term_info <- suppressMessages(
    select(GO.db, keys = res$GOID, columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
  )
  res <- merge(res, term_info, by = "GOID", all.x = TRUE)
  res$TERM[is.na(res$TERM)] <- res$GOID[is.na(res$TERM)]
  res$ONTOLOGY[is.na(res$ONTOLOGY)] <- "Unknown"
  res <- res[order(res$PValue), ]
  
  # 7. Export Full Results
  write.table(res, file = paste0("CAMERA_", cultivar_name, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  # 8. Filter and Plot Significant Terms
  sig <- res[res$FDR < 0.05, ]
  top <- data.frame() 
  
  # Safety check: Only plot if there are significant results
  if (nrow(sig) > 0) {
    top <- head(sig[order(sig$FDR), ], 20)
    top$TERM <- factor(top$TERM, levels = rev(unique(top$TERM)))
    
    p_bar <- ggplot(top, aes(x = TERM, y = -log10(FDR), fill = Direction)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = c(Up = "firebrick", Down = "steelblue")) +
      labs(title = cultivar_name, x = "GO term", y = "-log10(FDR)") +
      theme_bw() +
      theme(axis.text.y = element_text(size = 6))
    
    ggsave(paste0("barplot_", tolower(cultivar_name), "_hisat2_camera.png"), p_bar, width = 8, height = 6)
  } else {
    message("No significant GO terms (FDR < 0.05) found for ", cultivar_name)
  }
  
  # Attach cultivar name for downstream combining
  res$cultivar <- cultivar_name
  sig$cultivar <- cultivar_name
  if(nrow(top) > 0) top$cultivar <- cultivar_name
  
  return(list(res = res, sig = sig, top = top))
}

# =========================================================================
# --- EXECUTE ANALYSIS ---
# =========================================================================

results_Catuai <- run_camera_analysis("Catuai", counts, metadata, term2gene, min_set_size)
results_CR95 <- run_camera_analysis("CR95", counts, metadata, term2gene, min_set_size)

# =========================================================================
# --- CROSS-CULTIVAR COMPARISONS (BUTTERFLY CHART) ---
# =========================================================================

top_Catuai <- results_Catuai$top
top_CR95 <- results_CR95$top
comparison_terms <- unique(c(top_Catuai$GOID, top_CR95$GOID))

# Proceed only if we have overlapping or top terms to compare
if (length(comparison_terms) > 0) {
  
  comparison_Catuai <- results_Catuai$res[results_Catuai$res$GOID %in% comparison_terms, ]
  comparison_CR95 <- results_CR95$res[results_CR95$res$GOID %in% comparison_terms, ]
  comparison_df <- rbind(comparison_Catuai, comparison_CR95)
  
  # Logically order the Y-axis by minimum FDR across both cultivars
  term_order <- aggregate(FDR ~ TERM, data = comparison_df, FUN = min)
  comparison_df$TERM <- factor(comparison_df$TERM, 
                               levels = term_order$TERM[order(term_order$FDR, decreasing = TRUE)])
  
  # Calculate plotting score: Left (negative) for Catuai, Right (positive) for CR95
  comparison_df$plot_score <- -log10(comparison_df$FDR)
  comparison_df$plot_score[comparison_df$cultivar == "Catuai"] <- -comparison_df$plot_score[comparison_df$cultivar == "Catuai"]
  
  # Calculate symmetric axis limits based on max significance
  max_score <- max(abs(comparison_df$plot_score), na.rm = TRUE)
  
  p_butterfly <- ggplot(comparison_df, aes(x = TERM, y = plot_score, fill = Direction)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
    geom_col(width = 0.7, color = "black", size = 0.2) +
    coord_flip() +
    scale_fill_manual(values = c("Up" = "firebrick", "Down" = "steelblue")) +
    scale_y_continuous(labels = abs, 
                       limits = c(-max_score * 1.05, max_score * 1.05)) +
    labs(title = "Significant GO Terms: Catuai vs. CR95",
         subtitle = "Enriched in Catuai       |       Enriched in CR95",
         x = "GO Term", 
         y = expression("-log"[10]*"(FDR)"),
         fill = "Regulation") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8),
          plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 11),
          legend.position = "right")
  
  ggsave("butterfly_compare_hisat2_camera.png", p_butterfly, width = 10, height = 7)
}