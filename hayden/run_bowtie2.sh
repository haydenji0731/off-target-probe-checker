#!/usr/bin/env bash

set -x

# static vars that require attention
bowtie="/ccb/salz4-3/hji20/off-target-probe-checker/tools/bowtie2"
index_dir="/ccb/salz4-3/hji20/off-target-probe-checker/indices"

tgt=$1
qry=$2
out_dir=$3
flags="${@4}" # TODO: finish implementing (make index building optional)
p=36

# parse flags
for flag in $flags; do
   if [[ "$flag" == "-i" ]]; then
      if [[ $# -gt 0 ]]; then
            index_dir="$1"
            echo "index_dir set to: $index_dir"
            shift
        else
            echo "Error: No directory provided after -i"
            exit 1
        fi
   fi
done

tgt_prefix=$(basename "$tgt" .fa)
tgt_prefix="${tgt_prefix%.fasta}"
qry_prefix=$(basename "$qry" .fa | basename "$qry" .fasta)

# indexing command (optional)
"${bowtie}/bowtie2-build" $tgt "${index_dir}/${tgt_prefix}" --threads $p

"${bowtie}/bowtie2" -f -x "${index_dir}/${tgt_prefix}" -U $qry \
   --very-sensitive -k 1000 --threads $p | \
   samtools view -b -o "${out_dir}/${qry_prefix}.bam" -@ $p