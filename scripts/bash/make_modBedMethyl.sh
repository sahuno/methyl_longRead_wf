#!/bin/bash

#create modified bedmethyl file with columns in the following order:
#chr, start, end, strand, methylation level, coverage, 

#original columns:
#c("chrom", "start", "end", "mod_base", "score", "strand", "start2", "end2", "color", "Nvalid_cov", "fraction_modified", "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall")

# Check if input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <bedmethyl_file>"
    exit 1
fi

input_file="$1"
output_file="$2"
awk '{
    # Print columns in the desired order
    print $1, $2, $3, $6, $4, $10, $11, $12, $13, $14, $15, $16, $17, $18
    # print $1, $2, $3, $6, $4, $5, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18

}' "$input_file" > "$output_file"
