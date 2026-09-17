##############################################################
# WGCNA CO-EXPRESSION ANALYSIS
##############################################################

# Load the packages required for WGCNA and Cytoscape connection
library(WGCNA)
library(dplyr)
library(RCy3)

# Allow WGCNA to use multiple processor threads
allowWGCNAThreads()

##############################################################
# 1. PREPARE EXPRESSION DATA FOR WGCNA
##############################################################

# Transpose the VST matrix so samples are rows and genes are columns, as required by WGCNA
expr.full <- t(vst.mat)

# Calculate the expression variance of each gene across all samples
gene.var.wgcna <- apply(expr.full,2,var)

# Select a maximum of 5,000 genes for network construction
top.n.genes <- min(5000,ncol(expr.full))

# Rank genes by variance and retain the 5,000 most variable genes
top.genes <- names(sort(gene.var.wgcna,decreasing=TRUE))[1:top.n.genes]

# Create the expression matrix containing only the selected variable genes
datExpr <- expr.full[,top.genes,drop=FALSE]

# Check for genes or samples with missing values or insufficient information
gsg <- goodSamplesGenes(datExpr,verbose=0)

# Remove problematic samples or genes if WGCNA identifies any
if(!gsg$allOK) datExpr <- datExpr[gsg$goodSamples,gsg$goodGenes,drop=FALSE]

##############################################################
# 2. SELECT SOFT THRESHOLD POWER
##############################################################

# Define the range of candidate soft threshold powers
powers <- 1:20

# Evaluate candidate powers to identify one compatible with approximate scale free topology
sft <- pickSoftThreshold(datExpr,powerVector=powers,verbose=0)

# Store the automatically estimated soft threshold power
soft.power <- sft$powerEstimate

# Use a power of 6 if WGCNA does not automatically identify a suitable value
if(is.na(soft.power)) soft.power <- 6

# Print the selected soft threshold power
cat("\nSelected soft threshold power:",soft.power,"\n")

##############################################################
# 3. CONSTRUCT THE WGCNA NETWORK
##############################################################

# Build a signed weighted co-expression network and identify gene modules
net <- blockwiseModules(datExpr,power=soft.power,TOMType="signed",networkType="signed",minModuleSize=30,reassignThreshold=0,mergeCutHeight=0.25,numericLabels=FALSE,pamRespectsDendro=FALSE,verbose=0)

# Extract the module assigned to each gene
moduleColors <- net$colors

# Assign gene IDs as names of the module color vector
names(moduleColors) <- colnames(datExpr)

# Display the number of genes assigned to each module
print(table(moduleColors))

##############################################################
# 4. PREPARE SAMPLE TRAITS
##############################################################

# Match metadata to the samples present in the WGCNA expression matrix
metadata.wgcna <- metadata[rownames(datExpr),,drop=FALSE]

# Create a factor combining cultivar and treatment for each sample
group.factor <- interaction(metadata.wgcna$cultivar,metadata.wgcna$treatment,sep="_")

# Convert the four cultivar by treatment groups into numeric trait variables
traits <- model.matrix(~0+group.factor)

# Remove the group.factor prefix from trait column names
colnames(traits) <- sub("^group.factor","",colnames(traits))

# Set sample names as row names of the trait matrix
rownames(traits) <- rownames(datExpr)

# Convert the trait matrix to a data frame
traits <- as.data.frame(traits)

# Create an additional numerical variable representing the Cultivar × Treatment interaction pattern
traits$Interaction_Effect <- ifelse(metadata.wgcna$cultivar=="CR95" & metadata.wgcna$treatment=="xylella",1,ifelse(metadata.wgcna$cultivar=="CR95" | metadata.wgcna$treatment=="xylella",-0.5,0))

# Extract and order the module eigengenes
MEs <- orderMEs(net$MEs)

##############################################################
# 5. TEST MODULE ASSOCIATION WITH CULTIVAR × TREATMENT
##############################################################

# Create an empty table to store the interaction statistics for each module
interaction.results <- data.frame(ME=colnames(MEs),interaction_F=NA,interaction_pvalue=NA)

# Fit a linear model for each module eigengene and test the Cultivar × Treatment interaction
for(i in seq_len(ncol(MEs))){model <- lm(MEs[,i]~cultivar*treatment,data=metadata.wgcna);model.anova <- anova(model);interaction.results$interaction_F[i] <- model.anova["cultivar:treatment","F value"];interaction.results$interaction_pvalue[i] <- model.anova["cultivar:treatment","Pr(>F)"]}

# Adjust interaction p values for multiple testing using the Benjamini Hochberg method
interaction.results$interaction_FDR <- p.adjust(interaction.results$interaction_pvalue,method="BH")

# Remove the ME prefix to obtain the module name
interaction.results$module <- sub("^ME","",interaction.results$ME)

# Remove the grey module because it contains genes not assigned to a defined co-expression module
interaction.results <- interaction.results[interaction.results$module!="grey",,drop=FALSE]

# Order modules from lowest to highest interaction p value
interaction.results <- interaction.results[order(interaction.results$interaction_pvalue),,drop=FALSE]

# Export the module interaction results
write.csv(interaction.results,file.path(outpath,"WGCNA_module_interaction_results.csv"),row.names=FALSE)

##############################################################
# 6. SELECT MODULE OF INTEREST
##############################################################

# Select the purple module for downstream analysis
best.module.name <- "purple"

# Obtain all genes assigned to the selected module
module.genes <- names(moduleColors[moduleColors==best.module.name])

# Print the selected module and its number of genes
cat("\nSelected module for interaction analysis:",best.module.name,"\nTotal genes in module:",length(module.genes),"\n")

##############################################################
# 7. CALCULATE MODULE MEMBERSHIP AND INTEGRATE DESEQ2 RESULTS
##############################################################

# Calculate the correlation between each gene and each module eigengene
MM <- cor(datExpr,MEs,use="p")

# Extract module membership values for genes belonging to the purple module
MM.best <- MM[module.genes,paste0("ME",best.module.name)]

# Convert the DESeq2 Cultivar × Treatment interaction results to a data frame
res.int.df <- as.data.frame(res.list$Interaction)

# Create a table containing each gene, its module, module membership and absolute module membership
hub.table <- data.frame(gene_id=module.genes,module=rep(best.module.name,length(module.genes)),moduleMembership=as.numeric(MM.best),absMM=abs(as.numeric(MM.best)))

# Add the DESeq2 interaction log2 fold change for each gene
hub.table$deseq2_log2FC_interaction <- res.int.df[hub.table$gene_id,"log2FoldChange"]

# Add the DESeq2 interaction test statistic for each gene
hub.table$deseq2_stat_interaction <- res.int.df[hub.table$gene_id,"stat"]

# Add the DESeq2 interaction p value for each gene
hub.table$deseq2_pvalue_interaction <- res.int.df[hub.table$gene_id,"pvalue"]

# Add the DESeq2 interaction adjusted p value for each gene
hub.table$deseq2_padj_interaction <- res.int.df[hub.table$gene_id,"padj"]

# Classify genes as interaction DEGs when adjusted p value is below 0.05
hub.table$is_interaction_DEG <- !is.na(hub.table$deseq2_padj_interaction) & hub.table$deseq2_padj_interaction<0.05

# Rank genes from highest to lowest absolute module membership
hub.table <- hub.table[order(-hub.table$absMM),,drop=FALSE]

# Export the hub gene information
write.csv(hub.table,file.path(outpath,"WGCNA_hub_genes.csv"),row.names=FALSE)

##############################################################
# 8. PREPARE GENES FOR CYTOSCAPE
##############################################################

# Select the 150 genes with the highest module membership
top.hubs <- head(hub.table$gene_id,150)

# Identify genes that are both members of the selected module and significant interaction DEGs
module.interaction.degs <- intersect(module.genes,interaction.genes)

# Combine top hub genes and interaction DEGs into a single gene list
network.genes <- unique(c(top.hubs,module.interaction.degs))

# Limit the network to a maximum of 200 genes
max.network.genes <- 200

# If more than 200 genes were selected, retain the first 200 according to hub gene ranking
if(length(network.genes)>max.network.genes) network.genes <- head(hub.table$gene_id[hub.table$gene_id %in% network.genes],max.network.genes)

# Print the number of genes selected for Cytoscape export
cat("\nGenes exported to Cytoscape:",length(network.genes),"\n")

##############################################################
# 9. CALCULATE TOPOLOGICAL OVERLAP FOR THE NETWORK
##############################################################

# Calculate TOM similarity among the selected genes using the same signed network and soft threshold
TOM.network <- TOMsimilarityFromExpr(datExpr[,network.genes,drop=FALSE],power=soft.power,networkType="signed")

# Add gene IDs to the rows and columns of the TOM matrix
dimnames(TOM.network) <- list(network.genes,network.genes)

# Extract unique pairwise TOM similarity values
tom.values <- TOM.network[upper.tri(TOM.network)]

# Use the 95th percentile of TOM similarity as the threshold for retaining network connections
tom.threshold <- as.numeric(quantile(tom.values,probs=0.95,na.rm=TRUE))

##############################################################
# 10. EXPORT NETWORK FOR CYTOSCAPE
##############################################################

# Export the weighted network as edge and node files compatible with Cytoscape
exportNetworkToCytoscape(TOM.network,edgeFile=file.path(outpath,paste0("Cytoscape_edges_",best.module.name,".txt")),nodeFile=file.path(outpath,paste0("Cytoscape_nodes_",best.module.name,".txt")),weighted=TRUE,threshold=tom.threshold,nodeNames=network.genes,nodeAttr=rep(best.module.name,length(network.genes)))

##############################################################
# 11. CREATE CYTOSCAPE NODE ATTRIBUTES
##############################################################

# Keep hub gene information only for genes included in the Cytoscape network
cyt.node.attrs <- hub.table[hub.table$gene_id %in% network.genes,,drop=FALSE]

# Reorder the node information to match the network gene order
cyt.node.attrs <- cyt.node.attrs[match(network.genes,cyt.node.attrs$gene_id),,drop=FALSE]

# Create the nodeName column used by Cytoscape
cyt.node.attrs$nodeName <- cyt.node.attrs$gene_id

# Convert TOM similarities above the threshold into a binary adjacency matrix
adjacency <- TOM.network>tom.threshold

# Remove self connections from the adjacency matrix
diag(adjacency) <- FALSE

# Calculate the number of retained connections for each gene
cyt.node.attrs$network_degree <- rowSums(adjacency)

# Export the complete Cytoscape node attribute table
write.table(cyt.node.attrs,file=file.path(outpath,paste0("Cytoscape_node_attributes_",best.module.name,".txt")),sep="\t",row.names=FALSE,quote=FALSE)

##############################################################
# 12. CONVERT FILES AND SEND DATA TO CYTOSCAPE
##############################################################

# Define the path of the Cytoscape edge file
edges_file <- file.path(outpath,"Cytoscape_edges_purple.txt")

# Read the exported edge file
edges_df <- read.csv(edges_file,sep="",header=TRUE,check.names=FALSE)

# Convert the edge file to CSV format
write.csv(edges_df,file=file.path(outpath,"Cytoscape_edges_purple.csv"),row.names=FALSE)

# Define the path of the Cytoscape node file
nodes_file <- file.path(outpath,"Cytoscape_nodes_purple.txt")

# If the node file exists, read it and convert it to CSV format
if(file.exists(nodes_file)){nodes_df <- read.csv(nodes_file,sep="",header=TRUE,check.names=FALSE);write.csv(nodes_df,file=file.path(outpath,"Cytoscape_nodes_purple.csv"),row.names=FALSE)}

# Check that Cytoscape is open and available through RCy3
cytoscapePing()

# Read the Cytoscape edge information
edges_df <- read.csv(file.path(outpath,"Cytoscape_edges_purple.csv"))

# Read the Cytoscape node information
nodes_df <- read.csv(file.path(outpath,"Cytoscape_nodes_purple.csv"))

# Sort network connections by weight and retain the 100 strongest edges
top_edges <- edges_df %>% arrange(desc(weight)) %>% head(100)

# Identify all genes connected by the 100 strongest edges
connected_nodes <- unique(c(top_edges$fromNode,top_edges$toNode))

# Retain node information only for genes present in the strongest network connections
filtered_nodes <- nodes_df %>% filter(nodeName %in% connected_nodes)

# Add the selected node information to the active Cytoscape network table
loadTableData(data=filtered_nodes,data.key.column="nodeName",table.key.column="name")

# Print a message indicating that the process was completed
cat("\n============================================\nPROCESO COMPLETADO: TABLA CARGADA EN CYTOSCAPE\n============================================\n")
