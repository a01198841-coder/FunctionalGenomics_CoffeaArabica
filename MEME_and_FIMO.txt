BASH CODE FOR MEME AND FIMO
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

# The final gene IDs were saved in:
# WRKY7_candidate_targets.txt
