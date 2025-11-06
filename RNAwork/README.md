# RNA Secondary Structure Analysis Pipeline

## Overview

This pipeline analyzes the relationship between RNA secondary structure (stem vs. loop regions) and synonymous SNPs in *Drosophila melanogaster* transcripts. The analysis identifies how codon usage patterns differ between structured (stem) and unstructured (loop) regions of mRNA, and examines whether natural selection acts differently on synonymous sites depending on their structural context.

---

## Pipeline Scripts

### 1. `parallel_rna_prediction.py`

**Purpose**: Predicts RNA secondary structure for all transcripts in parallel using ViennaRNA.

**What it does**:
- Takes transcript DNA sequences as input (from pickle file)
- Converts DNA to RNA (T → U)
- Uses ViennaRNA's `RNA.fold()` to predict minimum free energy (MFE) secondary structures
- Classifies every position in each transcript as either:
  - **stem** (paired bases, shown as `(` or `)` in dot-bracket notation)
  - **loop** (unpaired bases, shown as `.` in dot-bracket notation)
- Processes transcripts in parallel using 20 workers for speed
- Saves pre-computed structures to pickle file for fast reuse

**Key features**:
- Multi-threaded processing (~20 workers)
- Convergence-based execution
- Pre-classification of all positions for rapid lookup
- Caching system to avoid re-computation

**Input**:
- `ZI_SynSNPs_transcripts.p` - Pickle file with transcript sequences
- `ZI_SynSNPs_transcript_pos.txt` - Table of SNP positions to determine which transcripts are needed

**Output**:
- `transcript_structures_parallel.p` - Pickle file containing:
  - Secondary structure (dot-bracket notation)
  - Minimum free energy (kcal/mol)
  - Position-by-position stem/loop classifications
  - Success/failure status for each transcript

**Runtime**: ~10-20 minutes for full *D. melanogaster* transcriptome

---

### 2. `analyze_synonymous_snps_stem_loop.py`

**Purpose**: Maps genomic SNP positions to transcript coordinates and determines their RNA structural context.

**What it does**:
- Reads VCF files with SnpEff-annotated synonymous SNPs
- Parses GTF file to build coordinate mapping for each transcript:
  - **Genomic position** → **Transcript position** (accounting for exons and strand)
  - **Transcript position** → **CDS position** (excluding UTRs)
  - **CDS position** → **Codon position** (frame 1/2/3 and absolute codon number)
- Looks up RNA structure classification (stem/loop) for each SNP position
- Handles both forward and reverse strand genes correctly
- Validates coordinate mappings and logs errors

**Key technical aspects**:
- Strand-aware coordinate transformation
- Multi-level position mapping (genomic → transcript → CDS → codon)
- Integration of RNA structure predictions with genomic variants
- Comprehensive error logging for troubleshooting

**Input files**:
- VCF files: `/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/codonpairs/ZI_VCFs/*.vcf`
- GTF annotation: `/mnt/d/genemod/better_dNdS_models/drosophila/Dmel_resources/Drosophila_melanogaster.BDGP6.54.114.chr.gtf`
- RNA structures: `full_transcript_RNA_2nd_structure.pkl` (from script #1)

**Output**:
- `synonymous_snps_stem_loop.tsv` - Table with columns:
  - CHROM, POS, REF, ALT - Variant information
  - CODONs - Codon change (e.g., "agC/agT")
  - CODON_POS - Position within codon (1, 2, or 3)
  - CODON_ABS_POS - Absolute codon number in CDS
  - **STEM_LOOP** - RNA structure classification
  - TRANSCRIPT - Transcript ID
  - GENE - Gene name
- `*_errors.log` - Detailed error log for debugging

**Usage**:
```bash
# Full analysis
python analyze_synonymous_snps_stem_loop.py

# Debug mode (first 100 SNPs only)
python analyze_synonymous_snps_stem_loop.py --debug

# Inspect RNA structure format
python analyze_synonymous_snps_stem_loop.py --inspect-rna-only
```

---

### 3. `analyze_codon_stemloop_frequencies.py`

**Purpose**: Statistical analysis to identify codons that are significantly enriched in stems vs. loops.

**What it does**:
- Reads the SNP classification table from script #2
- Counts how often each codon appears in stem vs. loop regions
- Performs statistical tests to identify significant enrichments:
  - **Binomial test** for stem vs. loop proportions
  - **Chi-square test** (when expected counts ≥ 5)
  - **FDR correction** (Benjamini-Hochberg) for multiple testing
- Calculates fold-change enrichment for each codon
- Identifies codons that are:
  - **Stem-enriched** (appear more often in structured regions)
  - **Loop-enriched** (appear more often in unstructured regions)

**Key outputs**:
- Statistical significance testing (p-values, FDR-corrected)
- Fold-change calculations (observed/expected ratios)
- Lists of significantly enriched codons in each category
- Summary statistics by codon and amino acid

**Input**:
- `synonymous_snps_stem_loop.tsv` (from script #2)
- Command-line parameters:
  - `--input` - Input file path
  - `--output` - Output CSV file
  - `--alpha` - Significance threshold (default: 0.05)
  - `--min-counts` - Minimum total counts to include (default: 10)

**Output**:
- `stem_loop_frequencies_by_codon_full_transcript.txt` - Contains:
  - Observed stem/loop counts for each codon
  - Expected proportions (genome-wide baseline)
  - Fold-change enrichments
  - P-values and FDR-corrected p-values
  - Significance flags

**Statistical approach**:
1. Calculate overall stem/loop proportions across all codons
2. For each codon, test if its stem/loop ratio differs from expected
3. Apply FDR correction for multiple hypothesis testing
4. Report enrichment magnitude and direction

---

## Pipeline Workflow

```
1. parallel_rna_prediction.py
   ↓
   [Predicts RNA structures for all transcripts]
   ↓
   transcript_structures_parallel.p

2. analyze_synonymous_snps_stem_loop.py
   ↓
   [Maps SNPs to transcripts & assigns stem/loop classification]
   ↓
   synonymous_snps_stem_loop.tsv

3. analyze_codon_stemloop_frequencies.py
   ↓
   [Statistical testing for codon enrichments]
   ↓
   stem_loop_frequencies_by_codon_full_transcript.txt
```

---

## Biological Question

**Do synonymous codon choices influence RNA secondary structure, and does natural selection act on this structural variation?**

The pipeline addresses this by:
1. Identifying which RNA regions are structured (stems) vs. unstructured (loops)
2. Determining which codons preferentially occur in each structural context
3. Analyzing site frequency spectra (SFS) to detect selection signatures
4. Testing whether mutations toward stem-enriched or loop-enriched codons show different frequency distributions

---

## Key Findings

From the comprehensive analysis using these scripts:

1. **Position 3 (wobble) dominates** - Synonymous changes at the third codon position have the largest effect on RNA structure
2. **A↔G transitions** have the highest structural impact (0.8875 avg magnitude)
3. **25.4% of synonymous transitions** completely switch stem/loop enrichment preference
4. **CAA→CAG** (Glutamine) has the largest structural effect (magnitude 1.3629)
5. Significant differences in site frequency spectra between stem-enriched and loop-enriched codon contexts

---

## Dependencies

- Python 3.7+
- ViennaRNA Python package (`RNA` module)
- pandas
- numpy
- scipy
- statsmodels (for FDR correction)

**Install ViennaRNA**:
```bash
pip install ViennaRNA
```

---

## Data Files Location

All scripts assume data files are in:
```
/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/codonpairs/RNAwork/
```

---

## Quick Start

```bash
# Navigate to RNAwork directory
cd /mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/codonpairs/RNAwork

# Step 1: Predict RNA structures (10-20 min)
python parallel_rna_prediction.py

# Step 2: Map SNPs to structures (5-10 min)
python analyze_synonymous_snps_stem_loop.py

# Step 3: Statistical analysis (<1 min)
python analyze_codon_stemloop_frequencies.py \
  --input synonymous_snps_stem_loop.tsv \
  --output stem_loop_frequencies_by_codon_full_transcript.txt \
  --alpha 0.05
```

---

## Related Scripts

The RNAwork directory contains additional scripts for:
- Transition analysis (`comprehensive_synonymous_analysis.py`)
- SFS generation (`sum_differential_sfs.py`, `sum_location_direction_sfs.py`)
- Null model simulations (`simulate_random_codon_distributions.py`)
- Full transcript analysis (including UTRs and introns)

See `PIPELINE_EXECUTION_GUIDE.md` for the complete analysis workflow.

---

## Citation

If using these scripts, please cite:
- ViennaRNA Package 2.0 (Lorenz et al., 2011)
- SnpEff (Cingolani et al., 2012)

---

## Author

Jody Hey Lab
Temple University
Generated with assistance from Claude (Anthropic)
Last updated: July 2025
