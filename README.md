# Dmel-Synonymous-Fitness

### Requirements:
1. Need to have `samtools` installed and available in the command-line.

2. Need to download *Drosophila melanogaster* genome build 6 from NCBI and place the .fna and .fai files in `fasta` folder.
```{bash}
.    
├── fasta
│   ├── GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna
│   ├── GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.fai
```

3. Add the vcfs in data for the respective population (DGRP2 or DPGP3)
```{bash}
.    
├── data
│   ├── dgrp2
│   ├── dpgp3
│   │   ├── vcfs
│   │   │   ├── ZI_Chr2L_remade_rooted_lifted_filtered_ann.vcf
│   │   │   ├── ZI_Chr2R_remade_rooted_lifted_filtered_ann.vcf
│   │   │   ├── ZI_Chr3L_remade_rooted_lifted_filtered_ann.vcf
│   │   │   ├── ZI_Chr3R_remade_rooted_lifted_filtered_ann.vcf
```

4. VCF file name should have `<'pop_code'>_Chr<'chrom'>_<suffix>.vcf`

5. The vcf suffix should be specified in the command-line with the argument `-z` or `--vcfsuffix`.

### How to run the script to pair synonymous sites with short-introns
You can run directly on the command line:
```{bash}
# Set environment variables

# the base folder where the vcf files are
VCFPATH="../data" 

# the vcf filename suffix (after <'pop_code'>_Chr<'chrom'>_' and before .vcf)
VCFSUFFIX="remade_rooted_lifted_filtered_ann" 

# where the paired SFS for each codon change will be saved (it is a basefolder, where each population will have a folder within)
OUTPUTPATH="../result" 

# Command to run the script that pair synonymous and and short-intron SNPs
python get_short_intron_paired_SNP_allele_counts_codons.py \
  -c +ANY \
  -f \
  -p ZI \
  -n 160 \
  -u ANY \
  -s codons.txt \ # File containing all 134 synonymous ordered synonymous codon changes
  -v $VCFPATH \
  -z $VCFSUFFIX \
  -o $OUTPUTPATH \
```

There is a helper bash script with the environment variables and the command already set.
Within the `src/` folder you can run the script:
```{bash}

# Make the script executable, by
chmod +x run_pairing.sh

# Then run it:
./run_pairing.sh
```