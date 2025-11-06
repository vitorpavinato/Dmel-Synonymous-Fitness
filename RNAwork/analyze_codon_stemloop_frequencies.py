#!/usr/bin/env python3
"""
Analysis of stem and loop frequencies for each of the 61 sense codons.
"""

import pickle
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact
from collections import defaultdict, Counter
import argparse

# Define the 61 sense codons (all codons except stop codons)
STOP_CODONS = {'TAA', 'TAG', 'TGA'}
ALL_CODONS = {
    ''.join([n1, n2, n3]) 
    for n1 in 'ATCG' 
    for n2 in 'ATCG' 
    for n3 in 'ATCG'
}
SENSE_CODONS = sorted(list(ALL_CODONS - STOP_CODONS))

def load_pickle_file(filename):
    """Load a pickle file and return the data."""
    try:
        with open(filename, 'rb') as f:
            data = pickle.load(f)
        print(f"Loaded {filename}: {len(data)} entries")
        return data
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None

def examine_data_structure(transcripts_data, structures_data):
    """Examine the structure of the loaded data."""
    print("\n=== Data Structure Examination ===")
    
    # Check transcripts data
    if transcripts_data:
        sample_key = list(transcripts_data.keys())[0]
        sample_transcript = transcripts_data[sample_key]
        print(f"Sample transcript key: {sample_key}")
        print(f"Sample transcript type: {type(sample_transcript)}")
        if isinstance(sample_transcript, str):
            print(f"Sample transcript length: {len(sample_transcript)}")
            print(f"Sample transcript start: {sample_transcript[:60]}...")
        else:
            print(f"Sample transcript: {sample_transcript}")
    
    # Check structures data  
    if structures_data:
        sample_key = list(structures_data.keys())[0]
        sample_structure = structures_data[sample_key]
        print(f"Sample structure key: {sample_key}")
        print(f"Sample structure type: {type(sample_structure)}")
        print(f"Sample structure: {sample_structure}")
        
        # If it's a dict, examine its contents
        if isinstance(sample_structure, dict):
            print("Structure dict keys:", list(sample_structure.keys()))
            for key, value in sample_structure.items():
                print(f"  {key}: {type(value)} - {str(value)[:100]}...")

def count_codon_stem_loop(transcripts_data, structures_data):
    """
    Count stem and loop occurrences for each codon.
    """
    codon_counts = {codon: {'stem': 0, 'loop': 0, 'other': 0} for codon in SENSE_CODONS}
    total_codons_processed = 0
    transcripts_processed = 0
    
    # Find common transcript IDs
    common_transcripts = set(transcripts_data.keys()) & set(structures_data.keys())
    print(f"Found {len(common_transcripts)} common transcripts")
    
    for transcript_id in common_transcripts:
        try:
            sequence = transcripts_data[transcript_id]
            structure_info = structures_data[transcript_id]
            
            # Check if structure prediction was successful
            if not isinstance(structure_info, dict) or not structure_info.get('success', False):
                continue
                
            # Use position_classifications if available, otherwise fall back to dot-bracket
            position_classifications = structure_info.get('position_classifications')
            
            if position_classifications:
                # Use the direct position classifications
                # Process codons (reading frame starts at position 0)
                for i in range(0, len(sequence) - 2, 3):
                    codon = sequence[i:i+3].upper()
                    
                    # Skip if not a sense codon
                    if codon not in SENSE_CODONS:
                        continue
                    
                    # Get classifications for the three positions of this codon
                    # Position numbers in the dict are 1-based, so add 1
                    pos1_class = position_classifications.get(i+1, 'other')
                    pos2_class = position_classifications.get(i+2, 'other') 
                    pos3_class = position_classifications.get(i+3, 'other')
                    
                    # Count stem and loop positions in this codon
                    stem_positions = sum(1 for cls in [pos1_class, pos2_class, pos3_class] if cls == 'stem')
                    loop_positions = sum(1 for cls in [pos1_class, pos2_class, pos3_class] if cls == 'loop')
                    
                    # Classify the codon based on majority structure
                    if stem_positions > loop_positions:
                        codon_counts[codon]['stem'] += 1
                    elif loop_positions > stem_positions:
                        codon_counts[codon]['loop'] += 1
                    else:
                        # Tie or other - count as other
                        codon_counts[codon]['other'] += 1
                    
                    total_codons_processed += 1
            
            else:
                # Fall back to dot-bracket notation
                structure = structure_info.get('structure', '')
                if not structure or len(structure) != len(sequence):
                    continue
                    
                # Process codons (reading frame starts at position 0)
                for i in range(0, len(sequence) - 2, 3):
                    codon = sequence[i:i+3].upper()
                    
                    # Skip if not a sense codon
                    if codon not in SENSE_CODONS:
                        continue
                    
                    # Get structure for the three positions of this codon
                    codon_structure = structure[i:i+3]
                    
                    # Classify based on structure
                    # Count stem ('(' or ')') and loop ('.' or other)
                    stem_count = codon_structure.count('(') + codon_structure.count(')')
                    loop_count = codon_structure.count('.')
                    
                    # Classify the codon based on majority structure
                    if stem_count > loop_count:
                        codon_counts[codon]['stem'] += 1
                    elif loop_count > stem_count:
                        codon_counts[codon]['loop'] += 1
                    else:
                        # Tie or other - count as other
                        codon_counts[codon]['other'] += 1
                    
                    total_codons_processed += 1
            
            transcripts_processed += 1
            
            if transcripts_processed % 1000 == 0:
                print(f"Processed {transcripts_processed} transcripts...")
                
        except Exception as e:
            print(f"Error processing transcript {transcript_id}: {e}")
            continue
    
    print(f"Processed {transcripts_processed} transcripts")
    print(f"Processed {total_codons_processed} total codons")
    
    return codon_counts

def perform_statistical_analysis(codon_counts):
    """
    Perform statistical analysis to test for differences in stem/loop frequencies.
    """
    results = []
    
    # Calculate overall totals
    total_stem = sum(counts['stem'] for counts in codon_counts.values())
    total_loop = sum(counts['loop'] for counts in codon_counts.values())
    total_other = sum(counts['other'] for counts in codon_counts.values())
    grand_total = total_stem + total_loop + total_other
    
    # Overall proportions (excluding 'other')
    total_stem_loop = total_stem + total_loop
    if total_stem_loop == 0:
        print("Error: No stem/loop data found!")
        return pd.DataFrame()
    
    overall_stem_prop = total_stem / total_stem_loop
    overall_loop_prop = total_loop / total_stem_loop
    
    print(f"\nOverall statistics:")
    print(f"Total stem: {total_stem:,}")
    print(f"Total loop: {total_loop:,}")
    print(f"Total other: {total_other:,}")
    print(f"Overall stem proportion: {overall_stem_prop:.4f}")
    print(f"Overall loop proportion: {overall_loop_prop:.4f}")
    
    # Test each codon
    for codon in SENSE_CODONS:
        counts = codon_counts[codon]
        stem_count = counts['stem']
        loop_count = counts['loop']
        other_count = counts['other']
        total_count = stem_count + loop_count + other_count
        stem_loop_total = stem_count + loop_count
        
        if stem_loop_total == 0:
            continue
        
        # Observed proportions
        obs_stem_prop = stem_count / stem_loop_total
        obs_loop_prop = loop_count / stem_loop_total
        
        # Expected counts based on overall proportions
        exp_stem_count = stem_loop_total * overall_stem_prop
        exp_loop_count = stem_loop_total * overall_loop_prop
        
        # Chi-square test
        chi2_pvalue = np.nan
        chi2_stat = np.nan
        try:
            if exp_stem_count >= 5 and exp_loop_count >= 5:
                observed = np.array([stem_count, loop_count])
                expected = np.array([exp_stem_count, exp_loop_count])
                chi2_stat, chi2_pvalue = chi2_contingency([observed, expected - observed])[:2]
        except:
            pass
        
        # Fisher's exact test (2x2 contingency table)
        fisher_pvalue = np.nan
        try:
            # Create 2x2 table: [this_codon_stem, this_codon_loop], [other_codons_stem, other_codons_loop]
            other_stem = total_stem - stem_count
            other_loop = total_loop - loop_count
            
            if other_stem >= 0 and other_loop >= 0:
                contingency_table = [[stem_count, loop_count], [other_stem, other_loop]]
                _, fisher_pvalue = fisher_exact(contingency_table)
        except:
            pass
        
        # Calculate fold changes
        stem_fold_change = (obs_stem_prop / overall_stem_prop) if overall_stem_prop > 0 else np.inf
        loop_fold_change = (obs_loop_prop / overall_loop_prop) if overall_loop_prop > 0 else np.inf
        
        results.append({
            'codon': codon,
            'stem_count': stem_count,
            'loop_count': loop_count,
            'other_count': other_count,
            'total_count': total_count,
            'stem_loop_total': stem_loop_total,
            'obs_stem_prop': obs_stem_prop,
            'obs_loop_prop': obs_loop_prop,
            'exp_stem_prop': overall_stem_prop,
            'exp_loop_prop': overall_loop_prop,
            'stem_fold_change': stem_fold_change,
            'loop_fold_change': loop_fold_change,
            'chi2_stat': chi2_stat,
            'chi2_pvalue': chi2_pvalue,
            'fisher_pvalue': fisher_pvalue
        })
    
    return pd.DataFrame(results)

def apply_multiple_testing_correction(results_df, alpha=0.05):
    """Apply FDR correction to p-values."""
    try:
        from statsmodels.stats.multitest import multipletests
        
        # Use Fisher's exact test p-values for correction
        valid_pvals = results_df['fisher_pvalue'].dropna()
        
        if len(valid_pvals) > 0:
            reject, pvals_corrected, _, _ = multipletests(valid_pvals, alpha=alpha, method='fdr_bh')
            
            # Map corrected p-values back
            results_df['fisher_pvalue_fdr'] = np.nan
            results_df.loc[results_df['fisher_pvalue'].notna(), 'fisher_pvalue_fdr'] = pvals_corrected
            results_df['significant_fdr'] = results_df['fisher_pvalue_fdr'] < alpha
        else:
            results_df['fisher_pvalue_fdr'] = np.nan
            results_df['significant_fdr'] = False
            
    except ImportError:
        print("Warning: statsmodels not available, skipping FDR correction")
        results_df['fisher_pvalue_fdr'] = results_df['fisher_pvalue']
        results_df['significant_fdr'] = results_df['fisher_pvalue'] < alpha
    
    return results_df

def main():
    parser = argparse.ArgumentParser(description='Analyze stem/loop frequencies by codon')
    parser.add_argument('--transcripts', 
                       default='/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/codonpairs/RNAwork/ZI_SynSNPs_transcripts.p',
                       help='Pickle file with transcript sequences')
    parser.add_argument('--structures',
                       default='/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/codonpairs/RNAwork/transcript_structures_parallel.p',
                       help='Pickle file with RNA structure predictions')
    parser.add_argument('-o', '--output',
                       default='/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/codonpairs/RNAwork/stem_loop_frequencies_by_codon.txt',
                       help='Output file')
    parser.add_argument('-a', '--alpha', type=float, default=0.05,
                       help='Significance level for FDR (default: 0.05)')
    parser.add_argument('--examine-only', action='store_true',
                       help='Only examine data structure, do not perform analysis')
    
    args = parser.parse_args()
    
    # Load data
    print("Loading transcript sequences...")
    transcripts_data = load_pickle_file(args.transcripts)
    
    print("Loading RNA structure predictions...")
    structures_data = load_pickle_file(args.structures)
    
    if not transcripts_data or not structures_data:
        print("Failed to load required data files!")
        return
    
    # Examine data structure
    examine_data_structure(transcripts_data, structures_data)
    
    if args.examine_only:
        return
    
    # Count codon occurrences
    print("\nCounting codon stem/loop occurrences...")
    codon_counts = count_codon_stem_loop(transcripts_data, structures_data)
    
    # Perform statistical analysis
    print("Performing statistical analysis...")
    results_df = perform_statistical_analysis(codon_counts)
    
    if len(results_df) == 0:
        print("No results generated!")
        return
    
    # Apply multiple testing correction
    results_df = apply_multiple_testing_correction(results_df, args.alpha)
    
    # Sort by significance
    results_df = results_df.sort_values('fisher_pvalue')
    
    # Save results
    results_df.to_csv(args.output, sep='\t', index=False)
    print(f"Results saved to {args.output}")
    
    # Print summary
    significant = results_df['significant_fdr'].sum()
    total = len(results_df)
    
    print(f"\nSummary (FDR < {args.alpha}):")
    print(f"Significant codons: {significant}/{total} ({100*significant/total:.1f}%)")
    
    if significant > 0:
        print(f"\nTop 10 most significant results:")
        display_cols = ['codon', 'stem_count', 'loop_count', 'obs_stem_prop', 'stem_fold_change', 
                       'fisher_pvalue', 'fisher_pvalue_fdr', 'significant_fdr']
        print(results_df[display_cols].head(10).to_string(index=False))
        
        # Summary of enrichments
        stem_enriched = results_df[(results_df['significant_fdr']) & (results_df['stem_fold_change'] > 1)]
        loop_enriched = results_df[(results_df['significant_fdr']) & (results_df['loop_fold_change'] > 1)]
        
        print(f"\nSignificantly stem-enriched codons ({len(stem_enriched)}):")
        for _, row in stem_enriched.head(20).iterrows():
            print(f"  {row['codon']}: {row['stem_fold_change']:.2f}x enriched (p={row['fisher_pvalue']:.2e})")
        
        print(f"\nSignificantly loop-enriched codons ({len(loop_enriched)}):")
        for _, row in loop_enriched.head(20).iterrows():
            print(f"  {row['codon']}: {row['loop_fold_change']:.2f}x enriched (p={row['fisher_pvalue']:.2e})")

if __name__ == "__main__":
    main()
