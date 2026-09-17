options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  "ggplot2",
  "VennDiagram",
  "gridExtra",
  "pheatmap",
  "dendsort"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "DESeq2",
  "edgeR",
  "limma",
  "GO.db",
  "AnnotationDbi",
  "WGCNA"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

library(DESeq2)
library(edgeR)
library(limma)
library(ggplot2)
library(VennDiagram)
library(grid)
library(gridExtra)
library(pheatmap)
library(dendsort)
library(GO.db)
library(AnnotationDbi)
library(WGCNA)
library(dplyr)
library(RCy3)

allowWGCNAThreads()


###########################
# 1. WORKING DIRECTORY
############################

setwd(choose.dir())

outpath <- "salido_complete"
dir.create(outpath, showWarnings = FALSE)


############################
# 2. LOAD COUNTS
############################

counts <- read.table(
  "Matrix_hisat2.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

# Clean gene IDs
rownames(counts) <- trimws(rownames(counts))

# Remove duplicated gene IDs
counts <- counts[!duplicated(rownames(counts)), ]

# Remove zero-count genes
counts <- counts[rowSums(counts) > 0, ]


############################
# 3. LOAD METADATA
############################

metadata <- read.table(
  "metadata.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

rownames(metadata) <- metadata$sample

# Clean sample names
colnames(counts) <- trimws(colnames(counts))
metadata$sample <- trimws(metadata$sample)

# Keep only matching samples
common.samples <- intersect(
  colnames(counts),
  metadata$sample
)

counts <- counts[, common.samples, drop = FALSE]

metadata <- metadata[
  match(common.samples, metadata$sample),
  ,
  drop = FALSE
]

rownames(metadata) <- metadata$sample

stopifnot(
  all(colnames(counts) == rownames(metadata))
)


############################
# 4. FACTORS
############################

metadata$cultivar <- factor(
  metadata$cultivar,
  levels = c("Catuai", "CR95")
)

metadata$treatment <- factor(
  metadata$treatment,
  levels = c("saline", "xylella")
)

stopifnot(!any(is.na(metadata$cultivar)))
stopifnot(!any(is.na(metadata$treatment)))


############################
# 5. FIXED BIOLOGICAL SAMPLE ORDER
############################

metadata$group <- interaction(
  metadata$cultivar,
  metadata$treatment,
  sep = "_"
)

desired.order <- c(
  "Catuai_saline",
  "Catuai_xylella",
  "CR95_saline",
  "CR95_xylella"
)

sample.order <- rownames(
  metadata[
    match(
      desired.order,
      metadata$group
    ),
    ,
    drop = FALSE
  ]
)

sample.order <- sample.order[
  !is.na(sample.order)
]

sample.order <- rownames(
  metadata[
    order(
      factor(
        metadata$group,
        levels = desired.order
      )
    ),
    ,
    drop = FALSE
  ]
)

counts <- counts[
  ,
  sample.order,
  drop = FALSE
]

metadata <- metadata[
  sample.order,
  ,
  drop = FALSE
]


############################
# 6. DESEQ2 OBJECT
############################

coldata <- metadata[
  ,
  c("cultivar", "treatment"),
  drop = FALSE
]

dds.full <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ cultivar + treatment + cultivar:treatment
)


############################
# 7. EXPRESSION FILTERING
############################

cpm <- counts(dds.full, normalized = FALSE)

group <- interaction(
  coldata$cultivar,
  coldata$treatment,
  drop = TRUE
)

keep <- rep(FALSE, nrow(cpm))

for (g in levels(group)) {
  idx <- which(group == g)
  cutoff <- ceiling(length(idx) / 2)
  keep <- keep |
    (
      rowSums(
        cpm[, idx, drop = FALSE] >= 10
      ) >= cutoff
    )
}

dds <- dds.full[keep, ]

cat(
  "\nGenes retained after filtering:",
  nrow(dds),
  "\n"
)


############################
# 8. RUN DESEQ2
############################

dds <- DESeq(dds)

cat("\nDESeq2 coefficients:\n")
print(resultsNames(dds))


############################
# 9. VST
############################

vst.obj <- vst(
  dds,
  blind = TRUE
)

vst.mat <- assay(vst.obj)

vst.mat <- vst.mat[
  ,
  sample.order,
  drop = FALSE
]


############################
# 10. SAMPLE ANNOTATION
############################

annotation.samples <- data.frame(
  Cultivar = metadata$cultivar,
  Treatment = metadata$treatment
)

rownames(annotation.samples) <- rownames(metadata)


##############################################################
# 14. ALL DESEQ2 CONTRASTS
##############################################################

res.list <- list()

res.list$Catuai <- results(
  dds,
  contrast = c("treatment", "xylella", "saline"),
  alpha = 0.05
)

interaction.coef <- grep(
  "cultivar.*CR95.*treatment.*xylella|treatment.*xylella.*cultivar.*CR95",
  resultsNames(dds),
  value = TRUE
)

if (length(interaction.coef) != 1) {
  stop("Could not uniquely identify interaction coefficient. Check resultsNames(dds).")
}

res.list$CR95 <- results(
  dds,
  contrast = list(c("treatment_xylella_vs_saline", interaction.coef)),
  alpha = 0.05
)

cultivar.coef <- grep(
  "^cultivar.*CR95.*Catuai",
  resultsNames(dds),
  value = TRUE
)

if (length(cultivar.coef) != 1) {
  stop("Could not uniquely identify cultivar coefficient. Check resultsNames(dds).")
}

res.list$CR95_vs_Catuai_saline <- results(
  dds,
  contrast = c("cultivar", "CR95", "Catuai"),
  alpha = 0.05
)

res.list$CR95_vs_Catuai_xylella <- results(
  dds,
  contrast = list(c(cultivar.coef, interaction.coef)),
  alpha = 0.05
)

res.list$Treatment_main <- results(
  dds,
  contrast = c("treatment", "xylella", "saline"),
  alpha = 0.05
)

res.list$Cultivar_main <- results(
  dds,
  contrast = c("cultivar", "CR95", "Catuai"),
  alpha = 0.05
)

dds.LRT <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ cultivar + treatment + cultivar:treatment
)

dds.LRT <- DESeq(
  dds.LRT,
  test = "LRT",
  reduced = ~ cultivar + treatment
)

res.list$Interaction <- results(
  dds.LRT,
  alpha = 0.05
)


##############################################################
# 15. DEG FUNCTION
##############################################################

alpha <- 0.05
lfc.cutoff <- 1

get.degs <- function(res) {
  x <- as.data.frame(res)
  x$padj[is.na(x$padj)] <- 1
  genes <- rownames(x)[x$padj < alpha & abs(x$log2FoldChange) >= lfc.cutoff]
  genes <- trimws(genes)
  genes <- unique(genes)
  genes
}


############################
# 16. DEG LISTS
############################

DEG.Catuai  <- get.degs(res.list$Catuai)
DEG.CR95    <- get.degs(res.list$CR95)
DEG.saline  <- get.degs(res.list$CR95_vs_Catuai_saline)
DEG.xylella <- get.degs(res.list$CR95_vs_Catuai_xylella)

interaction.genes <- rownames(res.list$Interaction)[
  !is.na(res.list$Interaction$padj) & res.list$Interaction$padj < alpha
]

interaction.genes <- unique(trimws(interaction.genes))


##############################################################
# 40. WGCNA PREPARACIÓN DE DATOS
##############################################################

expr.full <- t(vst.mat)

gene.var.wgcna <- apply(expr.full, 2, var)
top.n.genes <- min(5000, ncol(expr.full))
top.genes <- names(sort(gene.var.wgcna, decreasing = TRUE))[1:top.n.genes]

datExpr <- expr.full[, top.genes, drop = FALSE]

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}


##############################################################
# 42. WGCNA SOFT THRESHOLD
##############################################################

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)

soft.power <- sft$powerEstimate
if (is.na(soft.power)) {
  soft.power <- 6
}

cat("\nSelected soft threshold power:", soft.power, "\n")


##############################################################
# 43. BUILD WGCNA NETWORK
##############################################################

net <- blockwiseModules(
  datExpr,
  power = soft.power,
  TOMType = "signed",
  networkType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  verbose = 0
)

moduleColors <- net$colors
names(moduleColors) <- colnames(datExpr)

print(table(moduleColors))


##############################################################
# 44. WGCNA TRAITS
##############################################################

metadata.wgcna <- metadata[rownames(datExpr), , drop = FALSE]

group.factor <- interaction(metadata.wgcna$cultivar, metadata.wgcna$treatment, sep = "_")
traits <- model.matrix(~ 0 + group.factor)
colnames(traits) <- sub("^group.factor", "", colnames(traits))
rownames(traits) <- rownames(datExpr)

traits <- as.data.frame(traits)
traits$Interaction_Effect <- ifelse(
  metadata.wgcna$cultivar == "CR95" & metadata.wgcna$treatment == "xylella", 1,
  ifelse(metadata.wgcna$cultivar == "CR95" | metadata.wgcna$treatment == "xylella", -0.5, 0)
)

MEs <- orderMEs(net$MEs)


##############################################################
# 46. MODELO DE INTERACCIÓN GLOBAL (ANOVA)
##############################################################

interaction.results <- data.frame(
  ME = colnames(MEs),
  interaction_F = NA,
  interaction_pvalue = NA
)

for (i in seq_len(ncol(MEs))) {
  model <- lm(MEs[, i] ~ cultivar * treatment, data = metadata.wgcna)
  model.anova <- anova(model)
  
  interaction.results$interaction_F[i]      <- model.anova["cultivar:treatment", "F value"]
  interaction.results$interaction_pvalue[i] <- model.anova["cultivar:treatment", "Pr(>F)"]
}

interaction.results$interaction_FDR <- p.adjust(
  interaction.results$interaction_pvalue,
  method = "BH"
)

interaction.results$module <- sub("^ME", "", interaction.results$ME)
interaction.results <- interaction.results[interaction.results$module != "grey", , drop = FALSE]
interaction.results <- interaction.results[order(interaction.results$interaction_pvalue), , drop = FALSE]

write.csv(
  interaction.results,
  file.path(outpath, "WGCNA_module_interaction_results.csv"),
  row.names = FALSE
)


##############################################################
# 46B. SELECCIÓN DEL MÓDULO OBJETIVO
##############################################################

best.module.name <- "purple"
module.genes <- names(moduleColors[moduleColors == best.module.name])

cat(
  "\nMódulo seleccionado para análisis de interacción:", best.module.name,
  "\nTotal de genes en el módulo:", length(module.genes), "\n"
)


##############################################################
# 47. INTEGRACIÓN DE HUB GENES
##############################################################

MM <- cor(datExpr, MEs, use = "p")
MM.best <- MM[module.genes, paste0("ME", best.module.name)]

res.int.df <- as.data.frame(res.list$Interaction)

hub.table <- data.frame(
  gene_id = module.genes,
  module = rep(best.module.name, length(module.genes)),
  moduleMembership = as.numeric(MM.best),
  absMM = abs(as.numeric(MM.best))
)

hub.table$deseq2_log2FC_interaction <- res.int.df[hub.table$gene_id, "log2FoldChange"]
hub.table$deseq2_stat_interaction   <- res.int.df[hub.table$gene_id, "stat"]
hub.table$deseq2_pvalue_interaction <- res.int.df[hub.table$gene_id, "pvalue"]
hub.table$deseq2_padj_interaction   <- res.int.df[hub.table$gene_id, "padj"]

hub.table$is_interaction_DEG <- !is.na(hub.table$deseq2_padj_interaction) & hub.table$deseq2_padj_interaction < 0.05

hub.table <- hub.table[order(-hub.table$absMM), , drop = FALSE]

write.csv(
  hub.table,
  file.path(outpath, "WGCNA_hub_genes.csv"),
  row.names = FALSE
)


##############################################################
# 49. CYTOSCAPE EXPORT PREPARATION
##############################################################

top.hubs <- head(hub.table$gene_id, 150)
module.interaction.degs <- intersect(module.genes, interaction.genes)

network.genes <- unique(c(top.hubs, module.interaction.degs))
max.network.genes <- 200

if (length(network.genes) > max.network.genes) {
  network.genes <- head(
    hub.table$gene_id[hub.table$gene_id %in% network.genes],
    max.network.genes
  )
}

cat("\nGenes exported to Cytoscape:", length(network.genes), "\n")


##############################################################
# 50. TOM SIMILARITY
##############################################################

TOM.network <- TOMsimilarityFromExpr(
  datExpr[, network.genes, drop = FALSE],
  power = soft.power,
  networkType = "signed"
)

dimnames(TOM.network) <- list(network.genes, network.genes)

tom.values <- TOM.network[upper.tri(TOM.network)]
tom.threshold <- as.numeric(quantile(tom.values, probs = 0.95, na.rm = TRUE))


##############################################################
# 51. CYTOSCAPE FILES EXPORT
##############################################################

exportNetworkToCytoscape(
  TOM.network,
  edgeFile = file.path(outpath, paste0("Cytoscape_edges_", best.module.name, ".txt")),
  nodeFile = file.path(outpath, paste0("Cytoscape_nodes_", best.module.name, ".txt")),
  weighted = TRUE,
  threshold = tom.threshold,
  nodeNames = network.genes,
  nodeAttr = rep(best.module.name, length(network.genes))
)


##############################################################
# 52. CYTOSCAPE NODE ATTRIBUTES
##############################################################

cyt.node.attrs <- hub.table[hub.table$gene_id %in% network.genes, , drop = FALSE]
cyt.node.attrs <- cyt.node.attrs[match(network.genes, cyt.node.attrs$gene_id), , drop = FALSE]

cyt.node.attrs$nodeName <- cyt.node.attrs$gene_id

adjacency <- TOM.network > tom.threshold
diag(adjacency) <- FALSE

cyt.node.attrs$network_degree <- rowSums(adjacency)

write.table(
  cyt.node.attrs,
  file = file.path(outpath, paste0("Cytoscape_node_attributes_", best.module.name, ".txt")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


##############################################################
# 53. CONVERSIÓN A CSV E INYECCIÓN A CYTOSCAPE (SIN ESTILOS)
##############################################################

edges_file <- file.path(outpath, "Cytoscape_edges_purple.txt")
edges_df <- read.csv(edges_file, sep = "", header = TRUE, check.names = FALSE)

write.csv(
  edges_df, 
  file = file.path(outpath, "Cytoscape_edges_purple.csv"), 
  row.names = FALSE
)

nodes_file <- file.path(outpath, "Cytoscape_nodes_purple.txt")
if (file.exists(nodes_file)) {
  nodes_df <- read.csv(nodes_file, sep = "", header = TRUE, check.names = FALSE)
  write.csv(
    nodes_df, 
    file = file.path(outpath, "Cytoscape_nodes_purple.csv"), 
    row.names = FALSE
  )
}

# Cargar en Cytoscape
cytoscapePing()

edges_df <- read.csv(file.path(outpath, "Cytoscape_edges_purple.csv"))
nodes_df <- read.csv(file.path(outpath, "Cytoscape_nodes_purple.csv"))

top_edges <- edges_df %>%
  arrange(desc(weight)) %>%
  head(100)

connected_nodes <- unique(c(top_edges$fromNode, top_edges$toNode))
filtered_nodes <- nodes_df %>%
  filter(nodeName %in% connected_nodes)

loadTableData(
  data = filtered_nodes,
  data.key.column = "nodeName",
  table.key.column = "name"
)

cat("\n============================================\n")
cat("PROCESO COMPLETADO: TABLA CARGADA EN CYTOSCAPE\n")
cat("============================================\n")
