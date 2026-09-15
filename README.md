# FunctionalGenomics_CoffeaArabica

TERMINAL CODE FOR MEME AND FIMO
# Files for the terminal
CR95_induced_genes.txt
Catuai_induced_genes.txt
genomic.gff
GCF_036785885.1_Coffea_Arabica_ET-39_HiFi_genomic.fna
genoma.sizes
# Activate the environment
conda activate genomics-tools

# CR95 PROMOTER ANALYSIS

# Select the CR95 genes that were significantly induced according to DESeq2and keep only genes annotated as protein-coding
grep -Ff CR95_induced_genes.txt genomic.gff | awk -F'\t' '$3=="gene" && $9~/gene_biotype=protein_coding/' > CR95_induced_protein.gff

# Convert the selected CR95 genes from GFF to BED format, the start coordinate is changed to BED format by subtracting 1, the gene name and strand are also kept
awk -F'\t' 'BEGIN{OFS="\t"}{name="";n=split($9,a,";");for(i=1;i<=n;i++){if(a[i]~/^Name=/){name=a[i];sub(/^Name=/,"",name)}}if(name=="")name=$9;print $1,$4-1,$5,name,".",$7}' CR95_induced_protein.gff > CR95_induced_protein.bed

# Obtain the 1,000 bp region upstream of each CR95 protein-coding gene
# -l 1000 takes 1,000 bp upstream
# -r 0 does not include downstream sequence
# -s makes the extraction strand-specific
bedtools flank -i CR95_induced_protein.bed -g genoma.sizes -l 1000 -r 0 -s > CR95_induced_protein_upstream.bed

# Extract the DNA sequence of each CR95 promoter from the Coffea arabica reference genome
# -s keeps the correct strand orientation
# -name keeps the gene information in the FASTA header
bedtools getfasta -fi GCF_036785885.1_Coffea_Arabica_ET-39_HiFi_genomic.fna -bed CR95_induced_protein_upstream.bed -s -name -fo CR95_induced_protein_upstream.fasta

# Count how many CR95 promoter sequences were obtained - This gave 156 protein-coding promoter sequences
grep -c "^>" CR95_induced_protein_upstream.fasta

# Run MEME for de novo motif discovery in the CR95 promoters
# -dna indicates that the input contains DNA sequences
# -revcomp searches both the original sequence and its reverse complement
# -mod zoops allows zero or one occurrence of a motif per sequence
# -nmotifs 10 searches for a maximum of 10 motifs
# -minw 6 and -maxw 15 search for motifs between 6 and 15 bp
# -maxsize 500000 allows the total input sequence size used in the analysis
meme CR95_induced_protein_upstream.fasta -dna -revcomp -mod zoops -nmotifs 10 -minw 6 -maxw 15 -maxsize 500000 -oc meme_CR95_induced


# CATUAI PROMOTER ANALYSIS

# Select the Catuai genes that were significantly induced according to DESeq2 and keep only genes annotated as protein-coding
grep -Ff Catuai_induced_genes.txt genomic.gff | awk -F'\t' '$3=="gene" && $9~/gene_biotype=protein_coding/' > Catuai_induced_protein.gff

# Convert the selected Catuai genes from GFF to BED format, the chromosome, coordinates, gene name and strand are kept
awk -F'\t' 'BEGIN{OFS="\t"}{name="";n=split($9,a,";");for(i=1;i<=n;i++){if(a[i]~/^Name=/){name=a[i];sub(/^Name=/,"",name)}}if(name=="")name=$9;print $1,$4-1,$5,name,".",$7}' Catuai_induced_protein.gff > Catuai_induced_protein.bed

# Obtain the 1,000 bp strand-specific upstream region of each Catuai protein-coding gene
bedtools flank -i Catuai_induced_protein.bed -g genoma.sizes -l 1000 -r 0 -s > Catuai_induced_protein_upstream.bed

# Extract the promoter DNA sequences from the same Coffea arabica reference genome
bedtools getfasta -fi GCF_036785885.1_Coffea_Arabica_ET-39_HiFi_genomic.fna -bed Catuai_induced_protein_upstream.bed -s -name -fo Catuai_induced_protein_upstream.fasta

# Count the number of Catuai protein-coding promoter sequences obtained
grep -c "^>" Catuai_induced_protein_upstream.fasta

# Run MEME independently for the Catuai promoter sequences
# The same parameters used for CR95 were used so both analyses were comparable
meme Catuai_induced_protein_upstream.fasta -dna -revcomp -mod zoops -nmotifs 10 -minw 6 -maxw 15 -maxsize 500000 -oc meme_Catuai_induced


# FIMO - CR95

# Search the CR95 promoter sequences for occurrences of the motifs identified by MEME
# The MEME output file contains all motifs identified in the CR95 analysis
fimo --oc fimo_CR95_motifs meme_CR95_induced/meme.txt CR95_induced_protein_upstream.fasta

# After FIMO, only hits corresponding to CR95 Motif 3 were retained
# and significant occurrences were defined using q-value < 0.05

# These significant Motif 3 occurrences were then matched back to the
# CR95 promoter regions to obtain the corresponding Coffea arabica gene IDs

# The final gene IDs were saved in:
# WRKY7_candidate_targets.txt



TRANSCRIPTION FACTOR PREPARATION FOR FIGURES
# Convert DESeq2 results for each cultivar to data frames
CR95_df <- as.data.frame(res.CR95)
CR95_df$gene <- rownames(CR95_df)
Catuai_df <- as.data.frame(res.Catuai)
Catuai_df$gene <- rownames(Catuai_df)

# Keep one annotation row per gene to avoid duplicated genes during merging
annot_tf <- annotation[!duplicated(annotation$gene_id),c("gene_id","Preferred_name","Description","GOs","PFAMs")]
colnames(annot_tf)[1] <- "gene"

# Combine CR95 and Catuai differential-expression results in a single table
tf_table <- merge(CR95_df[,c("gene","log2FoldChange","padj")],Catuai_df[,c("gene","log2FoldChange","padj")],by="gene",suffixes=c("_CR95","_Catuai"))

# Add functional annotation to each gene
tf_table <- merge(tf_table,annot_tf,by="gene",all.x=TRUE)

# Search annotations for terms/domains commonly associated with plant transcription factors
tf_pattern <- "WRKY|MYB|NAC|NAM|NAP|bZIP|bHLH|AP2|EREBP|ERF|DREB|GATA|MADS|TCP|homeobox|HD-ZIP|C2H2|zinc finger|DOF|GRAS|ARF|SBP|SPL|YABBY|B3 domain|trihelix|heat shock factor|HSF|LBD|LOB domain|RAV|NF-Y|GARP|CAMTA|EIN3|EIL|BES1|BZR|SANT|transcription factor"

# Retain genes whose annotation matched the TF search pattern
TFs <- tf_table[grepl(tf_pattern,paste(tf_table$Preferred_name,tf_table$Description,tf_table$PFAMs),ignore.case=TRUE),]

# Remove duplicated gene IDs
TFs <- TFs[!duplicated(TFs$gene),]

# Identify significantly induced TFs in CR95
CR95_TFs_up <- TFs[!is.na(TFs$padj_CR95) & TFs$padj_CR95<0.05 & TFs$log2FoldChange_CR95>0,]

# Identify significantly induced TFs in Catuai
Catuai_TFs_up <- TFs[!is.na(TFs$padj_Catuai) & TFs$padj_Catuai<0.05 & TFs$log2FoldChange_Catuai>0,]


# FIGURE ___: Transcription factor response to Xylella 

# Convert repressed counts to negative values so they appear to the left of zero
tf_counts$value <- ifelse(tf_counts$direction=="Repressed",-tf_counts$n,tf_counts$n)

# Plot the number of significant TFs per family and direction
p_tf <- ggplot(tf_counts,aes(x=value,y=family,fill=direction))+geom_col(width=0.75)+geom_vline(xintercept=0)+geom_text(aes(label=n),hjust=ifelse(tf_counts$value<0,1.2,-0.2),size=3)+facet_grid(cultivar~.,scales="free_y",space="free_y")+labs(title="Transcription factor response to Xylella",x="Number of significant transcription factors",y="TF family",fill="Direction")+theme_bw()

# Display the figure
p_tf


# FIGURE ___: selected WRKY transcription factors

# Select the WRKY genes highlighted in the final figure
wrky_ids <- c("LOC113715139","LOC113724956","LOC113728125","LOC113731163","LOC113738517","LOC113700452","LOC113743770")

# Assign readable labels to each gene
wrky_labels <- c("LOC113715139"="WRKY7\nLOC113715139","LOC113724956"="WRKY22\nLOC113724956","LOC113728125"="WRKY22\nLOC113728125","LOC113731163"="WRKY25\nLOC113731163","LOC113738517"="WRKY40\nLOC113738517","LOC113700452"="WRKY51\nLOC113700452","LOC113743770"="WRKY34\nLOC113743770")

# Match the VST expression matrix columns with sample metadata
meta_plot <- metadata[match(colnames(vst_mat),metadata$sample),]

# Convert the expression values of the selected genes into long format for ggplot
wrky_long <- do.call(rbind,lapply(wrky_ids,function(g){data.frame(gene=wrky_labels[g],sample=colnames(vst_mat),VST=as.numeric(vst_mat[g,]),cultivar=meta_plot$cultivar,treatment=meta_plot$treatment)}))

# Define the order in which WRKY genes appear in the figure
wrky_long$gene <- factor(wrky_long$gene,levels=wrky_labels[wrky_ids])

# Generate expression plots for each WRKY gene in both cultivars
p_wrky <- ggplot(wrky_long,aes(x=treatment,y=VST))+geom_boxplot(width=0.5,outlier.shape=NA)+geom_jitter(width=0.06,size=1.8)+stat_summary(aes(group=1),fun=mean,geom="line",linewidth=0.8)+stat_summary(fun=mean,geom="point",size=2.5)+facet_grid(rows=vars(gene),cols=vars(cultivar),scales="free_y",switch="y")+labs(title="WRKY transcription factors in response to Xylella",x="Treatment",y="VST expression")+theme_bw()+theme(plot.title=element_text(size=16,hjust=0.5),strip.placement="outside",strip.text.y.left=element_text(angle=0,size=10),strip.text.x=element_text(size=11),axis.text=element_text(size=9),axis.title=element_text(size=11),panel.spacing=unit(0.25,"cm"),plot.margin=margin(15,15,15,150))

# Display the final figure
p_wrky

# Export as PDF
ggsave("Figure6B_WRKY_expression.pdf",p_wrky,width=12,height=12)


# FIGURE __: FIMO genes containing CR95 Motif 3

# Upload the file from FIMO/BEDtools
wrky7_targets <- readLines("WRKY7_candidate_targets.txt")

# Read gene IDs containing a significant CR95 Motif 3 FIMO occurrence
wrky7_targets <- readLines("WRKY7_candidate_targets.txt")

# Extract CR95 DESeq2 results for genes identified by FIMO
target_expr <- CR95_df[CR95_df$gene%in%wrky7_targets,c("gene","log2FoldChange","padj")]

# Prepare one functional annotation row per gene
target_annot <- annotation[!duplicated(annotation$gene_id),c("gene_id","Preferred_name","Description","GOs","PFAMs")]
colnames(target_annot)[1] <- "gene"

# Combine FIMO genes, differential-expression results and functional annotation
wrky7_target_table <- merge(target_expr,target_annot,by="gene",all.x=TRUE)

# Sort candidate genes by DESeq2 adjusted p-value
wrky7_target_table <- wrky7_target_table[order(wrky7_target_table$padj),]

# Search the FIMO genes for annotations related to plant defense/signaling
bio_pattern <- "defense|biotic|pathogen|immune|LRR|NBS|NB-ARC|ethylene|ACC oxidase|ACO|ERF|senescence|NAC|NAM|cell wall|xylan|cellulose|expansin|WRKY"

# Retain genes whose annotation matched these biological terms
wrky7_bio_targets <- wrky7_target_table[grepl(bio_pattern,paste(wrky7_target_table$Preferred_name,wrky7_target_table$Description,wrky7_target_table$GOs,wrky7_target_table$PFAMs),ignore.case=TRUE),]

# Inspect the selected annotated genes
wrky7_bio_targets[,c("gene","Preferred_name","Description","log2FoldChange","padj")]

# Select FIMO-identified genes shown in the final figure
target_ids <- c("LOC113691088","LOC113731163","LOC113724956","LOC113725207","LOC113728125","LOC113738517","LOC113694531","LOC113715139")

# Assign readable labels
target_labels <- c("LOC113691088"="ACO1\nLOC113691088","LOC113731163"="WRKY25\nLOC113731163","LOC113724956"="WRKY22\nLOC113724956","LOC113725207"="ERF5\nLOC113725207","LOC113728125"="WRKY22\nLOC113728125","LOC113738517"="WRKY40\nLOC113738517","LOC113694531"="pgip\nLOC113694531","LOC113715139"="WRKY7\nLOC113715139")

# Match VST matrix samples with metadata
meta_plot <- metadata[match(colnames(vst_mat),metadata$sample),]

# Convert expression values of the eight genes to long format
target_long <- do.call(rbind,lapply(target_ids,function(g){data.frame(gene=target_labels[g],sample=colnames(vst_mat),VST=as.numeric(vst_mat[g,]),cultivar=meta_plot$cultivar,treatment=meta_plot$treatment)}))

# Preserve the selected order of genes
target_long$gene <- factor(target_long$gene,levels=target_labels[target_ids])

# Plot expression in Catuai and CR95
p_targets <- ggplot(target_long,aes(x=treatment,y=VST))+geom_boxplot(width=0.5,outlier.shape=NA)+geom_jitter(width=0.06,size=1.8)+stat_summary(aes(group=1),fun=mean,geom="line",linewidth=0.8)+stat_summary(fun=mean,geom="point",size=2.5)+facet_grid(rows=vars(gene),cols=vars(cultivar),scales="free_y",switch="y")+labs(title="Expression of candidate WRKY7 target genes",subtitle="Expression of CR95 FIMO-identified candidate genes in Catuai and CR95",x="Treatment",y="VST expression")+theme_bw()+theme(plot.title=element_text(size=16,hjust=0.5),plot.subtitle=element_text(size=11,hjust=0.5),strip.placement="outside",strip.text.y.left=element_text(angle=0,size=9),strip.text.x=element_text(size=11),axis.text=element_text(size=9),axis.title=element_text(size=11),panel.spacing=unit(0.2,"cm"),plot.margin=margin(15,15,15,150))

# Display the plot
p_targets

# Export final figure
ggsave("Figure7_WRKY7_candidate_targets.pdf",p_targets,width=12,height=14)


