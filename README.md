Transcriptomic Analysis Workflow
This repository contains the workflow used to process RNA-seq data, identify differentially expressed genes, characterize transcription factors and promoter motifs, perform gene set enrichment analysis, and construct co-expression networks for the Catuaí (CAT) and CR95 coffee cultivars.

1. Data Pre-processing and Quantification
1.1. Prepare the RNA-seq datasets
Start with the raw paired-end RNA-seq datasets for the 15 samples:
Catuaí (CAT): 8 samples
CR95: 7 samples
Perform the initial quality control, read trimming, genome alignment, and read quantification using the Galaxy platform version 26.1.2.dev0 (The Galaxy Community, 2026).

1.2. Assess read quality
Use MultiQC to evaluate the quality and integrity of all raw sequencing datasets (Andrews, n.d.).
Pay particular attention to:
Adapter contamination.
Declining base quality toward the 3′ ends of reads.
Other potential sequencing-quality issues that may affect downstream analyses.

1.3. Trim low-quality and unwanted sequences
Use Trimmomatic to remove adapters and low-quality or repetitive sequences (Bolger et al., 2014).
Configure Trimmomatic with the following parameters:
IlluminaClip — remove Illumina universal adapter sequences.
Trailing: 30 — remove low-quality bases from the 3′ ends of reads.
Headcrop: 14 — remove the first 14 bases from the 5′ ends of reads.
After trimming, verify the resulting read quality and length. The expected output should contain reads with an average length of approximately 135 bp and quality scores of at least Q35.

1.4. Align reads to the reference genome
Use HISAT2 to align the cleaned paired-end reads against reference genome assembly GCF_036785885.1 (Kim et al., 2015).
For each sample, retain:
The HISAT2 alignment summary.
The resulting alignment file (.sam/.bam).
Use the alignment files for downstream gene-level quantification.

1.5. Quantify gene expression
Combine the HISAT2-generated BAM files with the corresponding gene annotation file.
Use featureCounts to quantify the number of reads assigned to genomic features (Liao et al., 2013).
The resulting output should be a raw read-count matrix, which will serve as the primary input for differential expression and downstream transcriptomic analyses.

2. Differential Expression Analysis
Perform the differential expression analysis in RStudio using R version 4.5.3 (R Core Team, 2026).
Prepare the following input files:
Raw read-count matrix.
Sample metadata file.
Custom gene annotation reference file.
Use the following R packages as required:
ggplot2
pheatmap
DESeq2
edgeR
GO.db
AnnotationDbi
ashr
EnhancedVolcano

Use ggplot2 for graphical outputs unless another package is specifically required for a particular visualization.

2.1. Prepare the gene annotation
Create a dictionary that maps the Gene Ontology (GO) identifiers contained in the custom annotation file to their corresponding biological descriptions.
Use this GO dictionary during the downstream functional interpretation of differentially expressed genes.

2.2. Filter lowly expressed genes
Before performing exploratory or differential expression analyses, filter the raw count matrix using the filterByExpr function from edgeR.
Apply a CPM cutoff of 10 to remove genes with negligible expression.
Use the resulting filtered count matrix for subsequent analyses.

2.3. Perform variance stabilization
Apply a variance stabilizing transformation (VST) to the filtered count matrix using DESeq2.
Use the VST-transformed matrix for exploratory analyses and visualization rather than the raw count matrix.

2.4. Evaluate sample clustering using PCA
Perform a principal component analysis (PCA) using the VST-transformed expression matrix.
Use the PCA to evaluate:
Sample clustering.
Similarity among biological replicates.
Potential sample outliers.
Major sources of variation within the dataset.

2.5. Perform differential expression analysis
Use DESeq2 to perform differential expression analysis.
Apply the Wald test to construct pairwise contrasts evaluating the effect of Xylella fastidiosa (XYL) relative to the control treatment (SAL).
Perform the comparisons separately for:
Catuaí: XYL vs. SAL.
CR95: XYL vs. SAL.

2.6. Filter differentially expressed genes
Filter the DESeq2 results according to statistical significance using:
FDR < 0.05
Create two DEG subsets based on absolute log2 fold change:
|log2FC| >= 0
and
|log2FC| >= 1
Use these subsets for subsequent visualization and functional analyses.

2.7. Visualize differential expression results
Generate volcano plots using EnhancedVolcano.
Generate hierarchical clustering heatmaps using pheatmap with the following settings:
Row-scaled expression values.
Euclidean distance.
Complete-linkage clustering.
Use these visualizations to compare expression patterns between the CAT and CR95 cultivars.

2.8. Perform functional annotation
Use the previously generated GO dictionary to associate GO terms with the DEG subsets.
For pathway visualization, generate KEGG Orthology (KO) lists for the DEGs and color-code the genes according to their expression direction.
Export these lists for visualization using KEGG Mapper Color, version 5.

3. Transcription Factor and Motif Analysis
Perform transcription factor (TF) and promoter motif analyses using Miniconda and a dedicated environment named:
genomics-tools
Run the analysis from a Bash terminal using Conda version 25.7.0.
Install the required tools through Bioconda (Grüning et al., 2018):
BEDTools
MEME Suite

3.1. Select genes for promoter analysis
Select the subset of CR95 DEGs that are annotated as protein-coding genes.
Use the reference genome annotation GCF_036785885.1 to identify the protein-coding genes.

3.2. Extract promoter regions
For each selected gene, use BEDTools version 2.31.1 to identify the strand-specific region located 1,000 bp upstream of the gene (Quinlan & Hall, 2010).
Use these genomic coordinates to extract the corresponding promoter sequences from the reference genome.
The resulting FASTA files should contain the promoter sequences used for motif discovery.

3.3. Perform de novo motif discovery
Use MEME version 5.5.9 to identify DNA motifs independently in the CAT and CR95 promoter datasets (Bailey et al., 2015).
Run the analysis with the following settings:
Analyze both DNA orientations.
Enable the reverse-complement option.
Use the zero-or-one occurrence per sequence (ZOOPS) model.
Search for a maximum of 10 motifs.
Set motif widths between 6 and 15 bp.
Set the maximum input size to 500,000 characters.
Perform motif discovery separately for:
CAT promoter sequences
CR95 promoter sequences

3.4. Compare discovered motifs with known TF-binding motifs
Use TOMTOM version 5.5.9 to compare the motifs identified by MEME with known plant transcription factor binding motifs (Gupta et al., 2007).
Use the following databases as references:
JASPAR CORE Plants
ARABIDOPSIS motif database (Rauluseviciute et al., 2024)
Use the resulting matches to identify potential transcription factors associated with the discovered promoter motifs.

3.5. Identify differentially expressed transcription factors
Use the DESeq2 DEG subsets for each cultivar and identify genes annotated as transcription factors.
Use the PFAM descriptions in the custom annotation file to identify TFs and assign them to their corresponding transcription factor families.
Group the identified significant TFs according to:
TF family.
Expression direction.
Use these groups to compare transcription factor responses between CAT and CR95.

3.6. Identify significant motif occurrences
Use FIMO version 5.5.9 to search the promoter sequences for occurrences of the significant motifs identified during the MEME analysis (Grand et al., 2011).
Filter FIMO results using:
q-value < 0.05
Retain only statistically significant motif occurrences.

3.7. Associate motifs with genes
Use BEDTools to match significant FIMO motif occurrences to their corresponding promoter regions.
Use these associations to identify the genes containing each significant motif.
Integrate these genes with:
DESeq2 differential expression results.
Functional annotation.
VST-transformed expression data.
Finally, compare the expression patterns of motif-containing genes between CAT and CR95 using box plots.

4. Gene Set Enrichment Analysis
Perform gene set enrichment analysis independently for:
Catuaí (CAT)
CR95
Use the following R packages:
limma
igraph
ggraph
Use limma for expression preprocessing and the CAMERA competitive gene set test, while igraph and ggraph are used to construct enrichment network visualizations.

4.1. Prepare the expression matrix
Start with the same filtered expression matrix generated using filterByExpr.
Normalize the library sizes using:
calcNormFactors
Then apply the voom transformation to model the relationship between mean expression and variance within the dataset.

4.2. Generate GO gene sets
Generate gene sets from the GO annotation using:
ids2indices
Retain only gene sets containing more than four genes.

4.3. Perform competitive gene set testing
Use the CAMERA function to test whether genes belonging to each GO gene set show coordinated expression patterns relative to genes outside the gene set.
Perform the analysis separately for CAT and CR95.
Consider a gene set significantly enriched when:
FDR < 0.05

4.4. Visualize enriched GO terms
For each cultivar, identify the 20 most significant GO terms.
Generate horizontal bar plots showing:
Significant GO terms.
Enrichment direction.
Relative significance.
Distinguish between upregulated and downregulated enrichment patterns.

4.5. Compare enrichment between cultivars
Use the Jaccard similarity index to quantify the overlap between genes belonging to significant GO terms:
Jaccard Index = intersection / union
Use this measurement to identify GO terms sharing genes between enrichment categories.

4.6. Generate the enrichment network
Use emap to construct an enrichment network map.
Connect GO terms when they share more than 20% of their genes.
Use igraph and ggraph to generate the network visualization.
The resulting network should provide a graphical representation of relationships among significantly enriched biological processes and highlight shared biological functions between CAT and CR95.

5. Co-expression Analysis
Perform co-expression analysis using the WGCNA package in R.
Start with the filtered and VST-transformed count matrix.

5.1. Construct the co-expression network
Calculate the Pearson correlation between all pairwise combinations of genes.
Apply a soft-thresholding power (β) to the correlation matrix to obtain a network with scale-free topology characteristics.
Use the resulting network to identify groups of genes with similar expression patterns.

5.2. Calculate topological overlap
Transform the network into a topological overlap measure (TOM) dissimilarity matrix.
Use the TOM-based dissimilarity matrix to evaluate the similarity between genes based on their shared network connections.

5.3. Identify gene modules
Apply hierarchical clustering using the average linkage method.
Identify gene modules using the Dynamic Tree Cut algorithm.
Each resulting module represents a group of genes with similar co-expression patterns.

5.4. Calculate module membership
Calculate module membership by determining the correlation between each individual gene expression profile and the corresponding module eigengene.
Use module membership values to identify genes that are strongly associated with their respective modules.

5.5. Integrate differential expression results
Integrate the DESeq2 DEG results with the WGCNA co-expression network.
Use this integration to prioritize candidate genes based on their combination of:
Differential expression.
Co-expression relationships.
Network connectivity.

5.6. Export the co-expression network
Use the WGCNA function:
exportNetworkToCytoscape
to export the network for visualization in Cytoscape Web version 1.0.8.
Filter and sort network edges according to their topological connection weight.
Retain the 100 strongest pairwise interactions for visualization.

5.7. Import network attributes into Cytoscape
Import the following information into Cytoscape:
Edge topology metadata.
Node attributes.
DESeq2 DEG status.
Use the RCy3 R package to map these attributes onto the active Cytoscape network.

5.8. Highlight differentially expressed genes
Create a discrete mapping rule based on the DEG_status attribute.
Use the following categories:
DEG_status = TRUE
for genes identified as DEGs, and
DEG_status = FALSE
for genes that are not classified as DEGs.
Use the resulting mapping to visually distinguish interaction-DEGs from non-DEGs within the co-expression network.



