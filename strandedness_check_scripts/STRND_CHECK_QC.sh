#!/bin/bash

input_directory=$1
type_=$2 #Are exon nts present?


folder_name=$(basename "$input_directory")

output_file="StrandednessCheck_QC_orthologgroups${folder_name}.txt"
output_file_probablt70="Strn_check_fail_files_probablt70_${folder_name}.txt"
echo -e "FILE\tMOTIF\tMean(Motif|Sense)\tStdev(Motif|Sense)\tTotal_sequences" > "$output_file"

for fasta_file in "$input_directory"/*.fa.aln; do
   python strandedness_check_Feb2026.py "$fasta_file" "$output_file" "$output_file_probablt70" "$type_"
   echo "$fasta_file"
done
