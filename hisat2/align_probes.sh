#!/bin/bash

# Set variables
TRANSCRIPTOME_FASTA="gencode.v47.transcripts.fa"
PROBE_FASTA="xenium_human_breast_gene_expression_panel_probe_sequences.fasta"
ALIGNER="hisat2"  # or "bowtie2"
INDEX_PREFIX="gencode_transcriptome"
ANNOTATION_GTF="gencode.v47.annotation.gtf"
GENE_BED="chess3.1.1.GRCh38.genes.fmted.bed"

# Create the transcriptome index
# echo "Building index..."
# if [[ "$ALIGNER" == "hisat2" ]]; then
#     hisat2-build $TRANSCRIPTOME_FASTA $INDEX_PREFIX
# else
#     bowtie2-build $TRANSCRIPTOME_FASTA $INDEX_PREFIX
# fi

# Align probes to transcriptome
echo "Aligning probes..."
if [[ "$ALIGNER" == "hisat2" ]]; then
    hisat2 -x $INDEX_PREFIX -f -U $PROBE_FASTA -S probes_transcriptome.sam --score-min L,0,-0.5
else
    bowtie2 -x $INDEX_PREFIX -f -U $PROBE_FASTA -S probes_transcriptome.sam
fi

# Convert SAM to BAM, sort, and index
echo "Processing BAM file..."
samtools view -S -b probes_transcriptome.sam | samtools sort -o probes_transcriptome.bam
samtools index probes_transcriptome.bam

# Remove alignments to the negative strand
# echo "Filtering positive-strand alignments..."
# samtools view -h -F 16 probes_transcriptome.bam -o probes_positive.bam
# samtools index probes_positive.bam

# Convert BAM to BED
echo "Converting BAM to BED..."
bedtools bamtobed -i probes_transcriptome.bam > probes_transcriptome.bed

# Download and extract the GTF annotation file if not already present
# if [[ ! -f "$ANNOTATION_GTF" ]]; then
#     echo "Downloading annotation file..."
#     wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
#     gunzip gencode.v47.annotation.gtf.gz
# fi

# # Extract transcript-to-genome mapping
# echo "Creating transcript-to-genome mapping..."
# awk '$3 == "transcript" {print $12, $1, $4, $5, $7}' $ANNOTATION_GTF | sed 's/"//g' | sed 's/;//g' > transcript_to_genome_map.tsv

# Convert transcriptome-aligned BED to genome coordinates
echo "Mapping transcriptome alignments to genome..."
awk 'NR==FNR {map[$1]=$2"\t"$3"\t"$4"\t"$5; next} {split($1, a, "|"); if (a[1] in map) print map[a[1]], $2, $3, $4}' transcript_to_genome_map.tsv probes_transcriptome.bed > probes_genome.bed

# Find overlaps with gene regions
echo "Finding overlaps with genes..."
bedtools intersect -a probes_genome.bed -b $GENE_BED -wa -wb > overlaps_xenium.bed

# Ensure unique locations in overlaps
# echo "Removing duplicate locations..."
# awk '{
#   key = $2 "\t" $3;
#   if (!(key in seen)) {
#     seen[key] = 1;
#     print;
#   }
# }' overlaps_xenium.bed > overlaps_xenium_unique.bed

# Extract gene names
echo "Extracting gene names..."
awk -F' ' '{split($7, a, "|"); print a[2], $11}' overlaps_xenium.bed > extracted_genes.txt

# Generate gene summary
echo "Generating gene match summary..."
awk '
BEGIN { FS=" "; OFS="\t"; print "Gene", "Matched_Count", "Total_Count", "Mismatches" }
{
  gene = $1;   # First column: the main gene name
  match_gene = $2;  # Second column: the matched gene name

  left[gene]++;  # Count occurrences of the gene from column 1
  if (match_gene == gene) {  
    right[gene]++;  # Count valid matches
  } else {
    mismatch[gene] = (mismatch[gene] ? mismatch[gene] "," match_gene : match_gene);  # Collect mismatches
  }
}
END {
  for (gene in left) {
    print gene, (gene in right ? right[gene] : 0), left[gene], (gene in mismatch ? mismatch[gene] : "None");
  }
}' extracted_genes.txt > gene_summary_transcriptome.txt



echo "Pipeline complete. Results:"
echo "  - Filtered alignments: probes_positive.bam"
echo "  - BED format alignments: probes_positive.bed"
echo "  - Mapped genomic alignments: probes_genome.bed"
echo "  - Gene overlaps: overlaps_xenium_unique.bed"
echo "  - Extracted genes: extracted_genes.txt"
echo "  - Gene match summary: gene_summary_transcriptome.txt"
