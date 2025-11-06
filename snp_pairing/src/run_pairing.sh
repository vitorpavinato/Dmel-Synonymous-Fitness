#!/bin/bash

VCFPATH="../data"
VCFSUFFIX="remade_rooted_lifted_filtered_ann"
OUTPUTPATH="../result"

# Command to run the script that pair synonymous and and short-intron SNPs
uv run python get_short_intron_paired_SNP_allele_counts_codons.py \
  -c +ANY \
  -f \
  -p ZI \
  -n 160 \
  -u ANY \
  -s codons.txt \
  -v $VCFPATH \
  -z $VCFSUFFIX \
  -o $OUTPUTPATH \