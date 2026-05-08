#!/usr/bin/env bash
set -euo pipefail

INPUT="/Users/ka13/My Drive/PhD_AF_mutations/Git_clone_GOF_project/Claude_code_make_AM_PTPN11/Data/AlphaMissense/AlphaMissense_hg38.tsv.gz"
OUTPUT="/Users/ka13/My Drive/PhD_AF_mutations/Git_clone_GOF_project/Claude_code_make_AM_PTPN11/Data/AlphaMissense/AlphaMissense_hg38_PTPN11.tsv"

# Extract the header line (first non-comment line starting with #CHROM), then filter rows for PTPN11 (UniProt Q06124)
{
  gzip -dc "$INPUT" | grep "^#CHROM"
  gzip -dc "$INPUT" | grep -v "^#" | awk -F'\t' '$6 == "Q06124"'
} > "$OUTPUT"

echo "Done. Output: $OUTPUT"
echo "Row count (excluding header): $(tail -n +2 "$OUTPUT" | wc -l | tr -d ' ')"
