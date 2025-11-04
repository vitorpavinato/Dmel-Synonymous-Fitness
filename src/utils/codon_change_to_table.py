import os
import re

def read_codon_pairs(filename):
    """Read the codon pairs from the input file."""
    with open(filename, 'r') as f:
        return [line.strip() for line in f]

def extract_metrics(filepath):
    """Extract the metrics from the results file."""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
            
            # Define patterns to match the metrics
            patterns = {
                'AIC': r'AIC\s+(-?\d+\.?\d*)',
                'likelihood': r'likelihood\s+(-?\d+\.?\d*)',
                'thetaratio': r'thetaratio\s+(-?\d+\.?\d*)\s+\((-?\d+\.?\d*)\s+-\s+(-?\d+\.?\d*)\)',
                '2Ns': r'2Ns\s+(-?\d+\.?\d*)\s+\((-?\d+\.?\d*)\s+-\s+(-?\d+\.?\d*)\)'
            }
            
            results = {}
            for metric, pattern in patterns.items():
                match = re.search(pattern, content)
                if match:
                    if metric in ['thetaratio', '2Ns']:
                        # Store the main value and confidence interval
                        results[metric] = (match.group(1), 
                                         f"({match.group(2)} - {match.group(3)})")
                    else:
                        results[metric] = match.group(1)
                else:
                    results[metric] = "NA"
            
            return results
    except Exception as e:
        print(f"Error processing {filepath}: {str(e)}")
        return None

def process_results(base_dir, codon_pairs_file):
    """Process all results files for given codon pairs."""
    # Read codon pairs
    codon_pairs = read_codon_pairs(codon_pairs_file)
    
    # Print header
    print("codon_pair\tAIC\tlikelihood\tthetaratio\tthetaratio_CI\t2Ns\t2Ns_CI")
    
    # Process each codon pair
    for codon_pair in codon_pairs:
        # Construct the expected filepath
        filepath = os.path.join(
            base_dir, 
            "results/snp_pairing", 
            "dpgp3_imputed_vcf_n197/codon_specific_results/", 
            codon_pair,
            f"ZI_{codon_pair}_197_+ANY-IMPUTED-SNPs_SYN_Qratio_fixed2Ns_nc196_estimates.out"
        )
        
        if os.path.exists(filepath):
            metrics = extract_metrics(filepath)
            if metrics:
                print(f"{codon_pair}\t"
                      f"{metrics['AIC']}\t"
                      f"{metrics['likelihood']}\t"
                      f"{metrics['thetaratio'][0]}\t"
                      f"{metrics['thetaratio'][1]}\t"
                      f"{metrics['2Ns'][0]}\t"
                      f"{metrics['2Ns'][1]}")
        else:
            print(f"{codon_pair}\tNA\tNA\tNA\tNA\tNA\tNA")

# Example usage:
if __name__ == "__main__":
    # These paths should be adjusted to match your actual directory structure
    base_dir = "/Users/vitorpavinato/WorkDir/SF_Ratios_syn"
    codon_pairs_file = "/Users/vitorpavinato/WorkDir/SF_Ratios_syn/results/snp_pairing/dpgp3_imputed_vcf_n197/codon_changes.txt"
    
    process_results(base_dir, codon_pairs_file)
