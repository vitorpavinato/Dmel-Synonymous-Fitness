# Equilibrium Codon Frequency Prediction

## Overview

These scripts predict the **equilibrium frequencies** of synonymous codons under mutation-selection balance. They use population genetic theory to model how mutation rates and selection coefficients combine to determine the long-term frequency distribution of codons for each amino acid.

---

## Biological Context

Synonymous codons (those coding for the same amino acid) are not expected to occur at equal frequencies because:

1. **Mutation bias**: Different nucleotide substitutions occur at different rates
2. **Selection on synonymous sites**: "Silent" mutations can affect:
   - Translation efficiency
   - mRNA stability
   - RNA secondary structure
   - Splicing regulation

The equilibrium frequency of each codon depends on the balance between:
- **Mutation** (pushing toward codons with higher mutation rates)
- **Selection** (favoring codons with higher fitness)

---

## Scripts

### 1. `calculate_equilibrium_frequencies.py`

**Purpose**: Calculate equilibrium codon frequencies using **mutation rates only** (neutral model).

**Approach**:
- Builds a mutation rate matrix **Q** for synonymous codons of each amino acid
- Uses two independent methods to find equilibrium:
  1. **Linear algebra**: Finds the eigenvector corresponding to eigenvalue = 0
  2. **Simulation**: Iteratively applies the transition matrix until convergence
- Compares both methods to verify accuracy

**Key equation**:
```
dπ/dt = π · Q = 0  (at equilibrium)
```

Where:
- **π** = vector of codon frequencies
- **Q** = mutation rate matrix (Q[i,j] = rate from codon i to codon j)

**Mutation rates used** (from published literature):
```python
Transition rates:
  A→G: 0.091836735
  G→A: 0.144132653
  C→T: 0.144132653
  T→C: 0.091836735

Transversion rates:
  A→C: 0.051020408
  A→T: 0.076530612
  C→G: 0.052295918
  G→T: 0.084183673
  (and reverse complements)
```

**Output**:
- `equilibrium_codon_frequencies.json` - Complete results for all amino acids
- Console output showing both methods agree within tolerance

**Usage**:
```bash
python calculate_equilibrium_frequencies.py
```

**Requirements**:
- numpy
- json
- `analyze_synonymous_mutations.py` (for genetic code definitions)

---

### 2. `Estimating_equilibrium_codon_frequencies_2.py`

**Purpose**: Calculate equilibrium codon frequencies using **both mutation rates AND selection coefficients** (selection-mutation balance).

**Approach**:
- Reads estimated selection coefficients (2Ns values) from model fitting results
- Uses mutation rates from *Drosophila* mutation accumulation studies
- Calculates equilibrium frequencies using the formula:

**Key equation**:
```
πᵢ = Z⁻¹ · exp(2·Nsᵢ) · ∏ⱼ≠ᵢ (uⱼᵢ/uᵢⱼ)
```

Where:
- **πᵢ** = equilibrium frequency of codon i
- **Nsᵢ** = population-scaled selection coefficient for codon i
- **uᵢⱼ** = mutation rate from codon i to codon j
- **Z** = normalization constant (sum of all πᵢ)

This formula combines:
1. **Selection term**: `exp(2·Nsᵢ)` - favors codons with higher fitness
2. **Mutation bias term**: Product of mutation rate ratios - favors codons that are easier to mutate TO

**Mutation rates used** (updated 7/24/2025 from Assaf et al. 2017):
```python
A→C: 0.03    C→A: 0.0925
A→G: 0.055   G→A: 0.2
A→T: 0.0575  T→A: 0.0575
C→G: 0.065   G→C: 0.065
C→T: 0.2     T→C: 0.055
G→T: 0.0925  T→G: 0.03
```

**Input**:
- Selection coefficient file from model fitting (e.g., `ZI_single_codon_pair_fitting_m5.out`)
- Contains estimated 2Ns values for each codon
- Format:
  ```
  Amino acid: L
  CTA    0.12345
  CTC   -0.23456
  CTG    0.34567
  ...
  ```

**Output**:
- Text file with predicted equilibrium frequencies for each codon
- Format: `Amino_Acid  Codon  Frequency`
- Example output file: `ZI_single_codon_pair_fitting_m5_shuffle_R_mm_codon_equilibrium_prediction_NEW_7_24_2025.out`

**Usage**:
```bash
# Edit script to set input/output file paths, then run:
python Estimating_equilibrium_codon_frequencies_2.py
```

**Note**: As commented in the script, "does not work very well" - the predictions may not match observed frequencies perfectly, indicating additional biological factors not captured by the simple model.

---

## Comparison of the Two Scripts

| Feature | calculate_equilibrium_frequencies.py | Estimating_equilibrium_codon_frequencies_2.py |
|---------|--------------------------------------|-----------------------------------------------|
| **Model** | Neutral (mutation only) | Selection-mutation balance |
| **Input** | Genetic code + mutation rates | Selection coefficients + mutation rates |
| **Selection** | None (assumes all codons equally fit) | Includes fitness differences (2Ns values) |
| **Methods** | Linear algebra + simulation | Analytical formula |
| **Use case** | Null expectation under no selection | Predicted equilibrium given estimated selection |
| **Output format** | JSON | Tab-delimited text |

---

## Theoretical Background

### Mutation-Selection Balance

In a population of effective size **N**, the equilibrium frequency of a codon depends on:

1. **Mutation pressure**: Rate at which mutations create and destroy that codon
2. **Selection pressure**: Fitness advantage/disadvantage (selection coefficient **s**)
3. **Genetic drift**: Random sampling effects (strength depends on **N**)

The key parameter is **2Ns** (or **4Ns** in diploids):
- If |2Ns| << 1: Mutation dominates (nearly neutral)
- If |2Ns| >> 1: Selection dominates (drift negligible)
- If |2Ns| ≈ 1: Both mutation and selection matter

### Multi-allele Wright-Fisher Model

For k synonymous codons, the equilibrium distribution follows from the stationary distribution of the mutation-selection-drift process. Under certain assumptions (weak selection, reversible mutation), this can be approximated analytically.

---

## Data Sources

**Mutation rates**:
- Schrider et al. 2013. *Genetics* 194:937-954
- Keightley et al. 2009. *Genome Research* 19:1195-1201
- **Updated rates**: Assaf et al. 2017. *Genome Research* 27:1988-2000

**Selection coefficients**:
- Estimated from site frequency spectrum (SFS) fitting
- Model fits differential selection on synonymous codons
- Input files generated by codon pair analysis pipeline

---

## Limitations and Caveats

1. **Simplified mutation model**: Assumes single-nucleotide mutations only
2. **Constant selection**: Assumes selection coefficients don't vary with context
3. **Homogeneous population**: Doesn't account for population structure
4. **No epistasis**: Treats each codon independently
5. **Equilibrium assumption**: Real populations may not be at equilibrium
6. **Missing factors**: Doesn't include:
   - tRNA abundance effects
   - mRNA structure constraints
   - Translational selection context
   - Gene expression level

As noted in the script: "does not work very well" - observed codon frequencies often deviate from predictions, suggesting these additional factors are important.

---

## Applications

These scripts are useful for:

1. **Null model generation**: What frequencies would we expect under neutral evolution?
2. **Testing selection hypotheses**: Do observed frequencies match predictions under estimated selection?
3. **Model validation**: Does the selection-mutation balance model explain codon usage?
4. **Comparative analysis**: How do predicted vs. observed frequencies differ across genes/species?

---

## Future Improvements

To better match observed frequencies, the model could be extended to include:
- Gene expression level (highly expressed genes show stronger codon bias)
- tRNA gene copy numbers (proxy for tRNA abundance)
- Local sequence context (neighboring codons affect mutation rates)
- mRNA structure (stem/loop regions as shown in RNA structure analysis)
- Recent demographic changes (populations not at equilibrium)

---

## Dependencies

**calculate_equilibrium_frequencies.py**:
```
numpy
json
analyze_synonymous_mutations (custom module)
```

**Estimating_equilibrium_codon_frequencies_2.py**:
```
numpy
scipy (for null_space, though not currently used)
matplotlib (imported but not used)
```

---

## Example Workflow

```bash
# 1. Calculate neutral expectation
python calculate_equilibrium_frequencies.py
# Output: equilibrium_codon_frequencies.json

# 2. Fit selection model to SFS data (separate pipeline)
# Generates: ZI_single_codon_pair_fitting_m5.out

# 3. Predict equilibrium with selection
python Estimating_equilibrium_codon_frequencies_2.py
# Output: *_codon_equilibrium_prediction.out

# 4. Compare predicted vs. observed frequencies
# (Analysis done in R or separate Python script)
```

---

## Related Work

These scripts are part of a larger project analyzing:
- Site frequency spectra (SFS) for synonymous variants
- Selection on codon pairs
- RNA secondary structure effects on codon usage
- Mutation-selection balance in *Drosophila melanogaster*

See also:
- `/codonpairs/manuscript/archive/RNAwork/` - RNA structure analysis
- `/codonpairs/single_codon_pair_work/` - Selection coefficient estimation

---

## Author

Jody Hey Lab
Temple University
Scripts written with assistance from Claude (Anthropic)
Date: July 24, 2025
