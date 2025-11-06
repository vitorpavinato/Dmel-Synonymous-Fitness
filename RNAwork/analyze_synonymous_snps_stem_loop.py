#!/usr/bin/env python3
"""
Analyze synonymous SNPs for stem/loop classification based on RNA secondary structure predictions.
Maps VCF positions to transcript positions accounting for exons and CDS regions.
"""

import pickle
import argparse
import sys
from collections import defaultdict
import re
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class TranscriptMapper:
    """Handles coordinate mapping between genomic and transcript positions"""
    
    def __init__(self, transcript_id, exons, cds_regions, strand):
        self.transcript_id = transcript_id
        self.exons = sorted(exons, key=lambda x: x[0])  # Sort by start position
        self.cds_regions = sorted(cds_regions, key=lambda x: x[0])
        self.strand = strand
        
        # Build mapping tables
        self._build_genomic_to_transcript_map()
        self._build_transcript_to_cds_map()
        
    def _build_genomic_to_transcript_map(self):
        """Build mapping from genomic position to transcript position"""
        self.genomic_to_transcript = {}
        transcript_pos = 0
        
        for exon_start, exon_end in self.exons:
            for genomic_pos in range(exon_start, exon_end + 1):
                if self.strand == '+':
                    self.genomic_to_transcript[genomic_pos] = transcript_pos
                    transcript_pos += 1
                else:  # Reverse strand
                    # For reverse strand, we still map sequentially but will need to reverse later
                    self.genomic_to_transcript[genomic_pos] = transcript_pos
                    transcript_pos += 1
        
        # For reverse strand, reverse the transcript positions
        if self.strand == '-':
            max_pos = transcript_pos - 1
            for genomic_pos in self.genomic_to_transcript:
                self.genomic_to_transcript[genomic_pos] = max_pos - self.genomic_to_transcript[genomic_pos]
    
    def _build_transcript_to_cds_map(self):
        """Build mapping from transcript position to CDS position"""
        self.transcript_to_cds = {}
        
        # First, map genomic CDS positions to transcript positions
        cds_transcript_positions = []
        for cds_start, cds_end in self.cds_regions:
            for genomic_pos in range(cds_start, cds_end + 1):
                if genomic_pos in self.genomic_to_transcript:
                    cds_transcript_positions.append(self.genomic_to_transcript[genomic_pos])
        
        # Sort CDS positions in transcript coordinates
        cds_transcript_positions.sort()
        
        # Map transcript positions to CDS positions (0-based for CDS)
        for cds_pos, transcript_pos in enumerate(cds_transcript_positions):
            self.transcript_to_cds[transcript_pos] = cds_pos
    
    def get_transcript_position(self, genomic_pos):
        """Get transcript position from genomic position"""
        return self.genomic_to_transcript.get(genomic_pos)
    
    def get_cds_position(self, transcript_pos):
        """Get CDS position from transcript position"""
        return self.transcript_to_cds.get(transcript_pos)
    
    def validate_mapping(self):
        """Perform consistency checks on the mapping"""
        errors = []
        
        # Check that exons don't overlap
        for i in range(len(self.exons) - 1):
            if self.exons[i][1] >= self.exons[i+1][0]:
                errors.append(f"Overlapping exons: {self.exons[i]} and {self.exons[i+1]}")
        
        # Check that CDS regions are within exons
        for cds_start, cds_end in self.cds_regions:
            cds_in_exon = False
            for exon_start, exon_end in self.exons:
                if cds_start >= exon_start and cds_end <= exon_end:
                    cds_in_exon = True
                    break
            if not cds_in_exon:
                # Check if CDS spans multiple exons
                cds_positions = set(range(cds_start, cds_end + 1))
                exon_positions = set()
                for exon_start, exon_end in self.exons:
                    exon_positions.update(range(exon_start, exon_end + 1))
                if not cds_positions.issubset(exon_positions):
                    errors.append(f"CDS region {cds_start}-{cds_end} not fully within exons")
        
        # Check transcript length matches exon lengths
        expected_length = sum(end - start + 1 for start, end in self.exons)
        actual_length = len(self.genomic_to_transcript)
        if expected_length != actual_length:
            errors.append(f"Transcript length mismatch: expected {expected_length}, got {actual_length}")
        
        return errors


def parse_gtf_file(gtf_file):
    """Parse GTF file to extract CDS and exon coordinates for each transcript"""
    transcript_data = defaultdict(lambda: {'exons': [], 'cds': [], 'strand': None, 'gene_id': None})
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom, source, feature, start, end, score, strand, frame, attributes = fields
            
            if feature not in ['exon', 'CDS']:
                continue
            
            # Parse attributes
            attr_dict = {}
            for attr in attributes.split('; '):
                if ' ' in attr:
                    key, value = attr.split(' ', 1)
                    attr_dict[key] = value.strip('"')
            
            if 'transcript_id' not in attr_dict:
                continue
            
            transcript_id = attr_dict['transcript_id']
            start = int(start)
            end = int(end)
            
            # Store data
            transcript_data[transcript_id]['strand'] = strand
            transcript_data[transcript_id]['gene_id'] = attr_dict.get('gene_id', 'unknown')
            
            if feature == 'exon':
                transcript_data[transcript_id]['exons'].append((start, end))
            elif feature == 'CDS':
                transcript_data[transcript_id]['cds'].append((start, end))
    
    # Convert defaultdict to regular dict
    return dict(transcript_data)


def load_rna_structure_data(pickle_file):
    """Load RNA secondary structure predictions from pickle file"""
    with open(pickle_file, 'rb') as f:
        return pickle.load(f)


def inspect_rna_structure_format(rna_structure_data, num_examples=3):
    """Inspect the format of RNA structure data"""
    logger.info("\nInspecting RNA structure data format:")
    
    # Get a few example transcripts
    transcript_ids = list(rna_structure_data.keys())[:num_examples]
    
    for transcript_id in transcript_ids:
        data = rna_structure_data[transcript_id]
        logger.info(f"\nTranscript {transcript_id}:")
        logger.info(f"  Type: {type(data)}")
        
        if isinstance(data, dict):
            logger.info(f"  Keys: {list(data.keys())}")
            
            # Special handling for position_classifications
            if 'position_classifications' in data:
                pos_class = data['position_classifications']
                logger.info(f"  position_classifications: list of length {len(pos_class) if isinstance(pos_class, list) else 'N/A'}")
                if isinstance(pos_class, list) and len(pos_class) > 0:
                    logger.info(f"    First 20 elements: {pos_class[:20]}")
                    unique_values = set(pos_class)
                    logger.info(f"    Unique values: {unique_values}")
                    logger.info(f"    Counts: stem={pos_class.count('stem')}, loop={pos_class.count('loop')}")
            
            # Show other fields
            for key, value in data.items():
                if key != 'position_classifications':  # Already handled above
                    if isinstance(value, list):
                        logger.info(f"  {key}: list of length {len(value)}")
                        if len(value) > 0 and key in ['structure', 'sequence']:
                            logger.info(f"    First few elements: {value[:5] if isinstance(value[0], str) else str(value[:5][:20])}")
                    elif isinstance(value, str) and len(value) > 50:
                        logger.info(f"  {key}: string of length {len(value)}")
                        logger.info(f"    First 50 chars: {value[:50]}")
                    else:
                        logger.info(f"  {key}: {type(value)}")
        elif isinstance(data, list):
            logger.info(f"  List length: {len(data)}")
            if len(data) > 0:
                logger.info(f"  First few elements: {data[:10]}")
                unique_values = set(data)
                logger.info(f"  Unique values in list: {unique_values}")


def parse_eff_annotation(eff_string):
    """Parse SnpEff EFF annotation string"""
    # EFF=SYNONYMOUS_CODING(LOW|SILENT|agC/agT|S36|665|CG12581|protein_coding|CODING|FBtr0078961|2|T)
    match = re.match(r'EFF=(\w+)\(([^)]+)\)', eff_string)
    if not match:
        return None
    
    effect_type = match.group(1)
    if effect_type != 'SYNONYMOUS_CODING':
        return None
    
    fields = match.group(2).split('|')
    if len(fields) < 10:
        return None
    
    return {
        'impact': fields[0],
        'functional_class': fields[1],
        'codon_change': fields[2],
        'aa_change': fields[3],
        'codon_number': int(fields[4]) if fields[4].isdigit() else None,
        'gene_name': fields[5],
        'biotype': fields[6],
        'coding_status': fields[7],
        'transcript_id': fields[8],
        'exon_number': int(fields[9]) if fields[9].isdigit() else None,
        'alt_allele': fields[10] if len(fields) > 10 else None
    }


def process_vcf_line(line, transcript_mappers, rna_structure_data, error_log):
    """Process a single VCF line and return SNP data if it's a synonymous SNP"""
    if line.startswith('#'):
        return None
    
    fields = line.strip().split('\t')
    if len(fields) < 8:
        return None
    
    chrom = fields[0]
    pos = int(fields[1])
    ref = fields[3]
    alt = fields[4]
    info = fields[7]
    
    # Check if biallelic
    if ',' in alt:
        return None
    
    # Find EFF annotation
    eff_match = re.search(r'EFF=[^;]+', info)
    if not eff_match:
        return None
    
    # Parse EFF annotation
    eff_data = parse_eff_annotation(eff_match.group(0))
    if not eff_data:
        return None
    
    transcript_id = eff_data['transcript_id']
    
    # Check if we have mapping data for this transcript
    if transcript_id not in transcript_mappers:
        error_log.append(f"No mapping data for transcript {transcript_id} at {chrom}:{pos}")
        return None
    
    mapper = transcript_mappers[transcript_id]
    
    # Get transcript position
    transcript_pos = mapper.get_transcript_position(pos)
    if transcript_pos is None:
        error_log.append(f"Position {pos} not found in transcript {transcript_id} exons")
        return None
    
    # Get CDS position
    cds_pos = mapper.get_cds_position(transcript_pos)
    if cds_pos is None:
        error_log.append(f"Transcript position {transcript_pos} not in CDS for {transcript_id}")
        return None
    
    # Get codon frame position (1, 2, or 3)
    codon_frame_pos = (cds_pos % 3) + 1
    
    # Get absolute codon position (1-based)
    absolute_codon_pos = (cds_pos // 3) + 1
    
    # Get stem/loop classification
    stem_loop = 'unknown'
    if transcript_id in rna_structure_data:
        # The RNA structure data should contain a dict with 'position_classifications' field
        # containing a list of 'stem'/'loop' strings, one for each base position in the full transcript
        rna_data = rna_structure_data[transcript_id]
        
        if isinstance(rna_data, dict) and 'position_classifications' in rna_data:
            stem_loop_list = rna_data['position_classifications']
            
            if isinstance(stem_loop_list, list) and transcript_pos < len(stem_loop_list):
                # RNA structure list is 0-based, transcript_pos is 0-based
                stem_loop = stem_loop_list[transcript_pos]
                # Validate that it's actually 'stem' or 'loop'
                if stem_loop not in ['stem', 'loop']:
                    error_log.append(f"Invalid stem/loop value '{stem_loop}' for {transcript_id} at position {transcript_pos}")
                    stem_loop = 'unknown'
            else:
                error_log.append(f"Transcript position {transcript_pos} out of range for RNA structure data (length {len(stem_loop_list) if isinstance(stem_loop_list, list) else 'N/A'}) for {transcript_id}")
        else:
            error_log.append(f"No 'position_classifications' field in RNA data for transcript {transcript_id}")
    else:
        error_log.append(f"No RNA structure data for transcript {transcript_id}")
    
    return {
        'CHROM': chrom,
        'POS': pos,
        'REF': ref,
        'ALT': alt,
        'CODONs': eff_data['codon_change'],
        'CODON_POS': codon_frame_pos,
        'CODON_ABS_POS': absolute_codon_pos,
        'STEM_LOOP': stem_loop,
        'TRANSCRIPT': transcript_id,
        'GENE': eff_data['gene_name']
    }


def main(args):
    """Main processing function"""
    # Load RNA structure data
    logger.info(f"Loading RNA structure data from {args.rna_pickle}")
    rna_structure_data = load_rna_structure_data(args.rna_pickle)
    logger.info(f"Loaded RNA structure data for {len(rna_structure_data)} transcripts")
    
    # In debug mode or inspect-only mode, inspect the RNA structure format
    if args.debug or args.inspect_rna_only:
        inspect_rna_structure_format(rna_structure_data)
    
    # If inspect-only mode, exit here
    if args.inspect_rna_only:
        return
    
    # Parse GTF file
    logger.info(f"Parsing GTF file {args.gtf_file}")
    transcript_data = parse_gtf_file(args.gtf_file)
    logger.info(f"Parsed GTF data for {len(transcript_data)} transcripts")
    
    # Build transcript mappers
    logger.info("Building transcript coordinate mappers")
    transcript_mappers = {}
    validation_errors = []
    
    for transcript_id, data in transcript_data.items():
        if not data['exons'] or not data['cds']:
            continue
        
        mapper = TranscriptMapper(
            transcript_id,
            data['exons'],
            data['cds'],
            data['strand']
        )
        
        # Validate mapping
        errors = mapper.validate_mapping()
        if errors:
            for error in errors:
                validation_errors.append(f"{transcript_id}: {error}")
        
        transcript_mappers[transcript_id] = mapper
    
    logger.info(f"Built mappers for {len(transcript_mappers)} transcripts")
    if validation_errors:
        logger.warning(f"Found {len(validation_errors)} validation errors")
    
    # Process VCF files
    results = []
    error_log = validation_errors.copy()
    snp_count = 0
    
    vcf_dir = Path(args.vcf_dir)
    vcf_files = list(vcf_dir.glob("*.vcf"))
    
    for vcf_file in vcf_files:
        logger.info(f"Processing {vcf_file}")
        
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                result = process_vcf_line(line, transcript_mappers, rna_structure_data, error_log)
                if result:
                    results.append(result)
                    snp_count += 1
                    
                    if args.debug and snp_count >= 100:
                        logger.info("Debug mode: Stopping after 100 SNPs")
                        break
        
        if args.debug and snp_count >= 100:
            break
    
    # Write results
    logger.info(f"Writing results to {args.output_file}")
    with open(args.output_file, 'w') as f:
        # Write header
        headers = ['CHROM', 'POS', 'REF', 'ALT', 'CODONs', 'CODON_POS', 'CODON_ABS_POS', 'STEM_LOOP', 'TRANSCRIPT', 'GENE']
        f.write('\t'.join(headers) + '\n')
        
        # Write data
        for result in results:
            row = [str(result[h]) for h in headers]
            f.write('\t'.join(row) + '\n')
    
    # Write error log
    if error_log:
        error_file = args.output_file.replace('.tsv', '_errors.log')
        logger.info(f"Writing error log to {error_file}")
        with open(error_file, 'w') as f:
            for error in error_log:
                f.write(error + '\n')
    
    logger.info(f"Processed {snp_count} synonymous SNPs")
    logger.info(f"Found {len(error_log)} errors/warnings")
    logger.info("Done!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze synonymous SNPs for stem/loop classification")
    parser.add_argument("--rna-pickle", default="/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/codonpairs/RNAwork/full_transcript_RNA_2nd_structure.pkl",
                        help="Path to RNA structure pickle file")
    parser.add_argument("--gtf-file", default="/mnt/d/genemod/better_dNdS_models/drosophila/Dmel_resources/Drosophila_melanogaster.BDGP6.54.114.chr.gtf",
                        help="Path to GTF file")
    parser.add_argument("--vcf-dir", default="/mnt/d/genemod/better_dNdS_models/popgen/Drosophila_SFS_and_SFRatios/codonpairs/ZI_VCFs",
                        help="Directory containing VCF files")
    parser.add_argument("--output-file", default="synonymous_snps_stem_loop.tsv",
                        help="Output file name")
    parser.add_argument("--debug", action="store_true",
                        help="Debug mode: process only first 100 synonymous SNPs")
    parser.add_argument("--inspect-rna-only", action="store_true",
                        help="Only inspect RNA structure format and exit")
    
    args = parser.parse_args()
    main(args)
