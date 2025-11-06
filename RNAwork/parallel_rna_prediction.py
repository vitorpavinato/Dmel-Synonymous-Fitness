#!/usr/bin/env python3
"""
Parallel RNA secondary structure prediction script.

Uses multiprocessing to predict structures for transcripts in parallel
and saves results to a pickle file for later use.
"""

import pickle
import RNA
import multiprocessing as mp
from pathlib import Path
import pandas as pd
from tqdm import tqdm
import time
import sys

def predict_structure_worker(transcript_data):
    """
    Worker function to predict structure for a single transcript.
    
    Args:
        transcript_data: tuple of (transcript_id, sequence)
    
    Returns:
        tuple of (transcript_id, structure_info_dict)
    """
    transcript_id, sequence = transcript_data
    
    try:
        # Convert DNA to RNA (T -> U)
        rna_sequence = sequence.replace('T', 'U')
        
        # Predict structure
        structure, mfe = RNA.fold(rna_sequence)
        
        return transcript_id, {
            'structure': structure,
            'mfe': mfe,
            'sequence_length': len(sequence),
            'success': True
        }
    except Exception as e:
        return transcript_id, {
            'structure': None,
            'mfe': None,
            'sequence_length': len(sequence),
            'success': False,
            'error': str(e)
        }

def classify_all_positions(structure):
    """
    Pre-classify all positions in a structure for faster lookup.
    
    Args:
        structure: dot-bracket structure string
        
    Returns:
        dict mapping 1-based position to 'stem'/'loop'/'unknown'
    """
    position_classifications = {}
    
    for i, char in enumerate(structure):
        pos_1based = i + 1
        if char in '()':
            position_classifications[pos_1based] = 'stem'
        elif char == '.':
            position_classifications[pos_1based] = 'loop'
        else:
            position_classifications[pos_1based] = 'unknown'
    
    return position_classifications

def process_transcript_chunk(chunk_data):
    """
    Process a chunk of transcripts in a single worker process.
    
    Args:
        chunk_data: tuple of (worker_id, transcript_list)
        
    Returns:
        dict of transcript_id -> structure_info
    """
    worker_id, transcript_list = chunk_data
    results = {}
    
    print(f"Worker {worker_id}: Processing {len(transcript_list)} transcripts")
    
    for transcript_id, sequence in transcript_list:
        transcript_id, structure_info = predict_structure_worker((transcript_id, sequence))
        
        # If successful, pre-classify all positions
        if structure_info['success'] and structure_info['structure']:
            structure_info['position_classifications'] = classify_all_positions(structure_info['structure'])
        else:
            structure_info['position_classifications'] = {}
        
        results[transcript_id] = structure_info
    
    print(f"Worker {worker_id}: Completed {len(results)} transcripts")
    return results

def main():
    print("Starting parallel RNA secondary structure prediction...")
    
    # File paths
    base_path = Path("/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/codonpairs/RNAwork")
    
    pickle_file = base_path / "ZI_SynSNPs_transcripts.p"
    input_table = base_path / "ZI_SynSNPs_transcript_pos.txt"
    output_structures_file = base_path / "transcript_structures_parallel.p"
    
    # Load transcript sequences
    print(f"Loading transcript sequences from: {pickle_file}")
    with open(pickle_file, 'rb') as f:
        transcript_sequences = pickle.load(f)
    print(f"Loaded {len(transcript_sequences)} transcript sequences")
    
    # Load SNP table to get unique transcripts needed
    print(f"Loading SNP table from: {input_table}")
    df = pd.read_csv(input_table, sep='\t')
    unique_transcripts_needed = df['transcriptID'].unique()
    print(f"Found {len(unique_transcripts_needed)} unique transcripts needed")
    
    # Filter to only transcripts we actually need and have sequences for
    transcripts_to_process = []
    missing_sequences = 0
    
    for transcript_id in unique_transcripts_needed:
        if transcript_id in transcript_sequences:
            transcripts_to_process.append((transcript_id, transcript_sequences[transcript_id]))
        else:
            missing_sequences += 1
    
    print(f"Will process {len(transcripts_to_process)} transcripts")
    if missing_sequences > 0:
        print(f"Warning: {missing_sequences} transcripts missing sequence data")
    
    # Split transcripts into chunks for parallel processing
    num_workers = 20
    chunk_size = len(transcripts_to_process) // num_workers + 1
    transcript_chunks = []
    
    for i in range(num_workers):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, len(transcripts_to_process))
        if start_idx < len(transcripts_to_process):
            chunk = transcripts_to_process[start_idx:end_idx]
            transcript_chunks.append((i, chunk))
    
    print(f"Split transcripts into {len(transcript_chunks)} chunks for {num_workers} workers")
    for i, (worker_id, chunk) in enumerate(transcript_chunks):
        print(f"  Worker {worker_id}: {len(chunk)} transcripts")
    
    # Process chunks in parallel
    start_time = time.time()
    print(f"Starting parallel processing with {num_workers} workers...")
    
    with mp.Pool(processes=num_workers) as pool:
        chunk_results = pool.map(process_transcript_chunk, transcript_chunks)
    
    # Combine results from all workers
    all_structures = {}
    for chunk_result in chunk_results:
        all_structures.update(chunk_result)
    
    end_time = time.time()
    processing_time = end_time - start_time
    
    print(f"Parallel processing completed in {processing_time:.1f} seconds")
    print(f"Processed {len(all_structures)} transcripts")
    
    # Count success/failure rates
    successful = sum(1 for info in all_structures.values() if info['success'])
    failed = len(all_structures) - successful
    success_rate = (successful / len(all_structures)) * 100
    
    print(f"Success rate: {successful}/{len(all_structures)} ({success_rate:.1f}%)")
    if failed > 0:
        print(f"Failed predictions: {failed}")
    
    # Calculate some statistics
    if successful > 0:
        mfe_values = [info['mfe'] for info in all_structures.values() if info['success'] and info['mfe'] is not None]
        if mfe_values:
            avg_mfe = sum(mfe_values) / len(mfe_values)
            min_mfe = min(mfe_values)
            max_mfe = max(mfe_values)
            print(f"MFE statistics: avg={avg_mfe:.2f}, min={min_mfe:.2f}, max={max_mfe:.2f} kcal/mol")
    
    # Save results to pickle file
    print(f"Saving structure predictions to: {output_structures_file}")
    with open(output_structures_file, 'wb') as f:
        pickle.dump(all_structures, f)
    
    print("Structure prediction complete!")
    print(f"Results saved to: {output_structures_file}")
    print(f"Average processing time per transcript: {processing_time/len(all_structures):.3f} seconds")

if __name__ == "__main__":
    # Check if RNA module is available
    try:
        import RNA
        print(f"ViennaRNA version: {RNA.__version__ if hasattr(RNA, '__version__') else 'unknown'}")
    except ImportError:
        print("Error: ViennaRNA Python package not found. Please install it first.")
        print("Try: pip install ViennaRNA")
        sys.exit(1)
    
    main()
