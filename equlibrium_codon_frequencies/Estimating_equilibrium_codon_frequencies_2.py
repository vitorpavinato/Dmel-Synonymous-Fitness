""" written by claude
    based on twocodon model in modelling_the_frequency_of_favored_and_disfavored_codons.nb
    I asked claude the following
        "I'm pasting here some mathematica cells that outline the equilibrium frequency of codons 
        for a two-fold amino acid as a function of relative mutation rates between them and the 
        selection coefficient.    I need to extend this model now for the equilibrium frequencies 
        of codons for a k-fold amino acid (k= 2,3,4 or 6).  In the two-fold case we do not need 
        the absolute selection coefficients as we end up just using the difference between the 
        absolute fitnesses,  which is s in one direction and -s in the other direction.  
        I suspect with k > 2  we might need to specify k absolute fitnesses,  
        or at least pick one to be 1, and the other k-1 to be parameterized. "
    then I asked about an algabraic solution and it gave this,  which looks sort of sensible 
    π₁ = Z^(-1) · e^(4Ns₁) · (u₂₁/u₁₂)·(u₃₁/u₁₃)·(u₄₁/u₁₄)
    π₂ = Z^(-1) · e^(4Ns₂) · (u₁₂/u₂₁)·(u₃₂/u₂₃)·(u₄₂/u₂₄)
    π₃ = Z^(-1) · e^(4Ns₃) · (u₁₃/u₃₁)·(u₂₃/u₃₂)·(u₄₃/u₃₄)
    π₄ = Z^(-1) · e^(4Ns₄) · (u₁₄/u₄₁)·(u₂₄/u₄₂)·(u₃₄/u₄₃)

    where Z = the sum   (i.e. normalizing )

    reads a file of estiamtes of 2Ns values and then given mutation rate values in 
        Two papers that estimate 
	        Schrider DR, Houle D, Lynch M, Hahn MW. 2013. Rates and genomic consequences of spontaneous mutational events in Drosophila melanogaster. Genetics 194:937-954.
	        Keightley PD, Trivedi U, Thomson M, Oliver F, Kumar S, Blaxter ML. 2009. Analysis of the genome sequences of three Drosophila melanogaster spontaneous mutation accumulation lines. Genome Research 19:1195-1201.
    generates predicted equilibrium codon frequencies

    does not work very well 
"""
import numpy as np
from scipy.linalg import null_space
import matplotlib.pyplot as plt


def get_codon_mutation_rate(from_codon, to_codon):
    """
    Returns the relative mutation rate from one codon to another.
    If codons differ by more than one base, returns 0.
    
    Parameters:
    from_codon (str): The source codon (3 nucleotides)
    to_codon (str): The target codon (3 nucleotides)
    
    Returns:
    float: The relative mutation rate, or 0 if codons are identical
           or differ by more than one base
    
    Raises:
    ValueError: If an invalid codon is provided
    """
    # Mutation rate matrix for individual nucleotides, these are relative values 
    mutation_rates = {
        'A': {'A': 0, 'C': 0.051020408, 'G': 0.091836735, 'T': 0.076530612},
        'C': {'A': 0.084183673, 'C': 0, 'G': 0.052295918, 'T': 0.144132653},
        'G': {'A': 0.144132653, 'C': 0.052295918, 'G': 0, 'T': 0.084183673},
        'T': {'A': 0.076530612, 'C': 0.091836735, 'G': 0.051020408, 'T': 0}
    }
#updated 7/24/2025 with rates from the summary in Assaf ZJ, Tilk S, Park J, Siegal ML, Petrov DA. 2017. Deep sequencing of natural and experimental populations of Drosophila melanogaster reveals biases in the spectrum of new mutations. Genome Research 27:1988-2000.
    mutation_rates = {
        'A': {'A': 0, 'C': 0.03, 'G': 0.055, 'T': 0.0575},
        'C': {'A': 0.0925, 'C': 0, 'G': 0.065, 'T': 0.2},
        'G': {'A': 0.2, 'C': 0.065, 'G': 0, 'T': 0.0925},
        'T': {'A': 0.0575, 'C': 0.055, 'G': 0.03, 'T': 0}
    }
    
    # Convert input to uppercase to handle lowercase inputs
    from_codon = from_codon.upper()
    to_codon = to_codon.upper()
    
    # Validate inputs
    if len(from_codon) != 3 or len(to_codon) != 3:
        raise ValueError("Both codons must be exactly 3 nucleotides long")
    
    valid_bases = ['A', 'C', 'G', 'T']
    for base in from_codon + to_codon:
        if base not in valid_bases:
            raise ValueError(f"Invalid nucleotide '{base}'. Must be one of: A, C, G, or T")
    
    # If codons are identical, return 0
    if from_codon == to_codon:
        return 0
    
    # Count the number of differences between codons
    diff_positions = []
    for i in range(3):
        if from_codon[i] != to_codon[i]:
            diff_positions.append(i)
    
    # If codons differ by more than one base, return 0
    if len(diff_positions) > 1:
        return 1e-10 # keeps code workable for 6 fold amino acids 
        # return 0
    
    # Return the mutation rate for the single base that differs
    position = diff_positions[0]
    from_base = from_codon[position]
    to_base = to_codon[position]
    
    return mutation_rates[from_base][to_base]

# Example usage
# print(get_codon_mutation_rate('ATG', 'ACG'))  # Should return the A->C rate: 0.051020408
# print(get_codon_mutation_rate('ATG', 'AAG'))  # Should return the T->A rate: 0.076530612
# print(get_codon_mutation_rate('ATG', 'ACC'))  # Should return 0 (differs by more than one base)
def build_codon_mutation_matrix(codons):
    """
    Builds a mutation rate matrix for a list of codons.
    
    Parameters:
    codons (list): List of codons
    
    Returns:
    numpy.ndarray: A k×k matrix of mutation rates between codons
    """
    import numpy as np
    
    k = len(codons)
    mutation_matrix = np.zeros((k, k))
    
    # Fill the matrix with mutation rates
    for i in range(k):
        for j in range(k):
            if i != j:  # Skip diagonal elements
                from_codon = codons[i]
                to_codon = codons[j]
                mutation_matrix[i, j] = get_codon_mutation_rate(from_codon, to_codon)
    
    return mutation_matrix

def get_codons_fitnesses(infilename):
    
    f = open(infilename,'r')
    lines = f.readlines()

    AAL = []
    codonLL = []
    sLL = []
    i = 0
    ai = 0
    while len(codonLL) < 18:
        if "Amino acid:" in lines[i]:
            AAL.append(lines[i].split()[2])
            codonLL.append([])
            sLL.append([])
            i += 1
            while len(lines[i]) > 2:
                if "---" not in lines[i]:
                    ls = lines[i].strip().split()
                    codonLL[ai].append(ls[0])
                    sLL[ai].append(float(ls[1]))
                i += 1
            ai += 1
        else:
            i += 1
    return codonLL,sLL,AAL


def main():
    # infilename = "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/codonpairs/ZI_single_codon_pair_fitting_m5.out"
    # outfilename = "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/codonpairs/ZI_single_codon_pair_fitting_m5_codon_equilibrium_prediction.out"

    infilename = "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/codonpairs/single_codon_pair_work/ZI/ZI_single_codon_pair_fitting_m5_shuffle_R_mm.out"
    outfilename = "/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/codonpairs/single_codon_pair_work/ZI/ZI_single_codon_pair_fitting_m5_shuffle_R_mm_codon_equilibrium_prediction_NEW_7_24_2025.out"
    of = open(outfilename,'w')
    of.write("generated by Estimating_equilibrium_codon_frequencies_2.py\n")
    of.write("input file: {}\n".format(infilename))
    of.write("\nAmino_Acid\tCodon\tFrequency\n")
    codonlists, slists,AAlist  = get_codons_fitnesses(infilename)

    for codons, slist, AA in zip(codonlists, slists, AAlist):
        k = len(codons)
        selection_coeffs = np.array(slist)
        mutation_rates = build_codon_mutation_matrix(codons)

        equi = []
        for ci in range(len(codons)):
            mprod = 1
            for j in range(len(codons)):
                if j != ci :
                    mprod *= mutation_rates[ci][j] / mutation_rates[j][ci]
                    # mprod *= mutation_rates[j][ci] / mutation_rates[ci][j]
            equi.append(np.exp(2 * slist[ci]) * mprod)
        Z = sum(equi)
        equi = equi/Z 
        print(AA, equi)
        for ci in range(len(codons)):
            of.write("{}\t{}\t{:.5f}\n".format(AA,codons[ci],equi[ci]))

        # exact, approx = compare_methods(selection_coeffs, mutation_rates, N)
        pass
    of.close()
    
        
    # # Example usage for a 4-fold amino acid
    # k = 4
    
    # # Selection coefficients (first is reference with s=0)
    # selection_coeffs = np.array([0.0, -0.001, -0.002, -0.003])
    # selection_coeffs = np.array([1.0, 0.999, 0.998, 0.997])
    
    # # Mutation rates matrix (from i to j)
    # # For simplicity, we'll use equal mutation rates in all directions
    # mutation_rates = np.ones((k, k)) * 1e-8
    # np.fill_diagonal(mutation_rates, 0)  # Zero diagonal
    
    # # Introduce some bias in mutation rates
    # mutation_rates[0, 1] = 2e-8  # Higher mutation rate from codon 1 to 2
    # mutation_rates[2, 3] = 3e-8  # Higher mutation rate from codon 3 to 4
    
    # mutation_rates = np.ones((k, k)) * 1e-5
    # np.fill_diagonal(mutation_rates, 0)  # Zero diagonal
    
    # # Introduce some bias in mutation rates
    # mutation_rates[0, 1] = 2e-5  # Higher mutation rate from codon 1 to 2
    # mutation_rates[2, 3] = 3e-5  # Higher mutation rate from codon 3 to 4

    # # Effective population size
    # N = 1000
    
    # # Calculate and compare equilibrium frequencies
    # exact, approx = compare_methods(selection_coeffs, mutation_rates, N)
    
    # return exact, approx

if __name__ == "__main__":
    main()