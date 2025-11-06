#!/bin/bash

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Configuration - consider moving to config file
readonly PYTHON_PATH="/Users/vitorpavinato/Library/Caches/pypoetry/virtualenvs/sfratios-kGgQ8yf2-py3.12/bin/python"
readonly POPULATIONS=("NC" "ZI")
readonly POPULATION_CODES=("dgrp2" "dpgp3")
readonly DOWNSAMPLE=160
readonly MODE="+ANY-IMPUTED-SNPs"
readonly WORKDIR="/Users/vitorpavinato/WorkDir/SF_Ratios_syn/results/snp_pairing"
readonly SFRATIOS="/Users/vitorpavinato/WorkDir/SFRatios/SFRatios.py"

# Command line arguments with better defaults
DEBUG_MODE="${DEBUG_MODE:-false}"
POPNAME="${POPNAME:-NC}"
PARALLEL_JOBS="${PARALLEL_JOBS:-1}"  # Allow parallel processing

# Function to display usage
usage() {
    echo "Usage: $0"
    echo "Environment variables:"
    echo "  POPNAME=<population>     Population name (default: NC)"
    echo "  DEBUG_MODE=<true/false>  Enable debug mode (default: false)"
    echo "  PARALLEL_JOBS=<number>   Number of parallel jobs (default: 1)"
    echo ""
    echo "Valid populations: ${POPULATIONS[*]}"
    exit 1
}

# Function to log messages with timestamps
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" >&2
}

# Function to validate population and get code
get_population_code() {
    local popname="$1"
    local popcode=""
    
    for i in "${!POPULATIONS[@]}"; do
        if [[ "${POPULATIONS[$i]}" == "$popname" ]]; then
            popcode="${POPULATION_CODES[$i]}"
            break
        fi
    done
    
    if [[ -z "$popcode" ]]; then
        log "ERROR: Invalid population name '$popname'. Valid options: ${POPULATIONS[*]}"
        exit 1
    fi
    
    echo "$popcode"
}

# Function to validate required files and directories
validate_environment() {
    local popcode="$1"
    local codon_changes_file="$2"
    
    # Check Python interpreter
    if [[ ! -x "$PYTHON_PATH" ]]; then
        log "ERROR: Python interpreter not found or not executable: $PYTHON_PATH"
        exit 1
    fi
    
    # Check SFRatios script
    if [[ ! -f "$SFRATIOS" ]]; then
        log "ERROR: SFRatios script not found: $SFRATIOS"
        exit 1
    fi
    
    # Check codon changes file
    if [[ ! -f "$codon_changes_file" ]]; then
        log "ERROR: Codon changes file not found: $codon_changes_file"
        exit 1
    fi
    
    # Check if file is empty
    if [[ ! -s "$codon_changes_file" ]]; then
        log "ERROR: Codon changes file is empty: $codon_changes_file"
        exit 1
    fi
}

# Function to process a single codon change
process_codon_change() {
    local codon_change="$1"
    local popname="$2"
    local inputdir="$3"
    local outputdir="$4"
    local index="$5"
    
    local codon_inputdir="${inputdir}/${codon_change}"
    local codon_outputdir="${outputdir}/${codon_change}"
    local input_file="${codon_inputdir}/${popname}_SI_and_SYN_${codon_change}_nc${DOWNSAMPLE}_${MODE}_SFSs.txt"
    
    log "Processing [$index] $codon_change"
    
    # Check if input file exists
    if [[ ! -f "$input_file" ]]; then
        log "WARNING: Input file not found, skipping: $input_file"
        return 1
    fi
    
    # Create output directory
    mkdir -p "$codon_outputdir"
    
    if [[ "$DEBUG_MODE" == "true" ]]; then
        echo "=== DEBUG MODE - Codon $codon_change ==="
        echo "Input file: $input_file"
        echo "Output dir: $codon_outputdir"
        echo "Command: $PYTHON_PATH $SFRATIOS -a \"$input_file\" -d fixed2Ns -f isfolded -g -i 5 -u -p \"${popname}_${codon_change}_${DOWNSAMPLE}_${MODE}\" -r \"$codon_outputdir\""
        echo "=== END DEBUG INFO ==="
        return 0
    fi
    
    # Execute SFRatios with error handling
    if ! "$PYTHON_PATH" "$SFRATIOS" \
        -a "$input_file" \
        -d fixed2Ns \
        -f isfolded \
        -g -i 5 -u \
        -p "${popname}_${codon_change}_${DOWNSAMPLE}_${MODE}" \
        -r "$codon_outputdir"; then
        log "ERROR: SFRatios failed for codon change: $codon_change"
        return 1
    fi
    
    log "Completed [$index] $codon_change"
    return 0
}

# Main execution
main() {
    log "Starting SFRatios processing script"
    log "Population: $POPNAME, Debug: $DEBUG_MODE, Parallel jobs: $PARALLEL_JOBS"
    
    # Get population code
    local popcode
    popcode=$(get_population_code "$POPNAME")
    log "Population code: $popcode"
    
    # Set up paths
    local inputdir="${WORKDIR}/${popcode}_imputed_vcf/codon_specific_results"
    local outputdir="${WORKDIR}/${popcode}_imputed_vcf/codon_specific_results"
    local codon_changes_file="${WORKDIR}/${popcode}_imputed_vcf/codon_changes.txt"
    
    # Validate environment
    validate_environment "$popcode" "$codon_changes_file"
    
    # Create output directory
    mkdir -p "$outputdir"
    
    # Read codon changes into array (more efficient)
    local codon_changes=()
    while IFS=',' read -r codon_change _; do
        [[ -n "$codon_change" ]] && codon_changes+=("$codon_change")
    done < "$codon_changes_file"
    
    local total_pairs=${#codon_changes[@]}
    log "Total codon pairs to process: $total_pairs"
    
    if [[ $total_pairs -eq 0 ]]; then
        log "ERROR: No codon changes found in file"
        exit 1
    fi
    
    # Process codon changes
    local success_count=0
    local failure_count=0
    
    if [[ "$PARALLEL_JOBS" -gt 1 ]]; then
        # Parallel processing using GNU parallel or xargs
        log "Processing in parallel with $PARALLEL_JOBS jobs"
        export -f process_codon_change log
        export PYTHON_PATH SFRATIOS DOWNSAMPLE MODE DEBUG_MODE
        
        printf '%s\n' "${codon_changes[@]}" | \
        nl -nln | \
        xargs -n2 -P"$PARALLEL_JOBS" -I{} bash -c '
            index=$(echo {} | cut -d" " -f1)
            codon_change=$(echo {} | cut -d" " -f2-)
            process_codon_change "$codon_change" "'"$POPNAME"'" "'"$inputdir"'" "'"$outputdir"'" "$index"
        '
    else
        # Sequential processing
        for i in "${!codon_changes[@]}"; do
            local index=$((i + 1))
            local codon_change="${codon_changes[$i]}"
            
            if process_codon_change "$codon_change" "$POPNAME" "$inputdir" "$outputdir" "$index"; then
                ((success_count++))
            else
                ((failure_count++))
            fi
        done
    fi
    
    # Summary
    log "Processing complete!"
    log "Successful: $success_count"
    log "Failed: $failure_count"
    log "Total: $total_pairs"
    
    if [[ $failure_count -gt 0 ]]; then
        log "WARNING: Some codon changes failed to process"
        exit 1
    fi
}

# Handle help flag
if [[ "${1:-}" == "-h" ]] || [[ "${1:-}" == "--help" ]]; then
    usage
fi

# Run main function
main "$@"
