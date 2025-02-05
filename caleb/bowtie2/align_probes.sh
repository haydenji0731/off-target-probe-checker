#!/bin/bash


# doesnt actually run, just a script to show how to align probes to the transcriptome

# Align probes to transcriptome
bowtie2 -x gencode_transcriptome -f -U xenium_human_breast_gene_expression_panel_probe_sequences.fasta -S probes_transcriptome.sam --very-sensitive -N 1 -L 10 --score-min L,0,0.5 -a --local


# susbetting the CIGAR string to have middle 30 BP matching
# chatgpt made this
gawk 'BEGIN {
    middle_start = 6; 
    middle_end = 35; 
}
# Pass header lines unchanged
/^@/ { print; next; }
{
    cigar = $6;      # Get the CIGAR string from SAM column 6
    qpos = 1;        # Query (probe) coordinate starts at 1
    valid = 1;       # Assume valid unless proven otherwise
    
    # Process the CIGAR string segment by segment.
    while (match(cigar, /^([0-9]+)([MIDNSHPX=])/, arr)) {
        n = arr[1] + 0;   # Number of bases for this segment
        op = arr[2];      # Operation letter
        
        seg_start = qpos; 
        # Operations that consume query bases (M, I, S, =, X) advance qpos.
        if (op ~ /[MIS=X]/) {
            seg_end = qpos + n - 1;
        } else {
            # D, N, H, or P do not consume query bases.
            seg_end = qpos;
        }
        
        # Determine if this segment overlaps the middle region [6,35].
        overlap_start = (seg_start < middle_start ? middle_start : seg_start);
        overlap_end   = (seg_end   < middle_end   ? seg_end   : middle_end);
        
        if (overlap_start <= overlap_end) {
            # There is overlap with the middle region.
            # For the middle region, only a match is acceptable.
            if (op != "M" && op != "=" && op != "X") {
                valid = 0;
                break;
            }
        }
        
        # Advance the query coordinate for operations that consume query bases.
        if (op ~ /[MIS=X]/)
            qpos += n;
        
        # Chop off the processed segment from the CIGAR string.
        cigar = substr(cigar, RSTART + RLENGTH);
    }
    
    # Finally, only print the SAM record if the middle region was perfectly aligned.
    if (valid == 1)
        print $0;
}' probes_transcriptome.sam > filtered_probes_transcriptome.sam

# can check with this
awk '{print $6}' filtered_probes_transcriptome.sam | sort | uniq -c | sort -nr


# filter out non-protein-coding transcripts
awk '$3 ~ /protein_coding/ { print }' filtered_probes_transcriptome.sam > protein_coding_filtered.sam

# get gene names only
awk '{
    # Split first column using "|" as the delimiter.
    split($1, a, "|");
    # Split third column using "|" as the delimiter.
    split($3, b, "|");
    # a[2] is the gene name from the first column.
    # b[6] is the gene name from the third column.
    print a[2], b[6];
}' protein_coding_filtered.sam > gene_names.txt



# get gene summary
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
}' gene_names.txt > gene_summary_transcriptome.txt







