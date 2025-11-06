"""
Wrapper script to process multiple SFS files at once.
"""

import os
import glob
import subprocess
import argparse


def parse_arguments():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Process multiple SFS files using reformat_sfs_file.py'
    )

    parser.add_argument(
        '-d', '--directory',
        default='.',
        help='Root directory to search for files (default: current directory)'
    )
    parser.add_argument(
        '-p', '--pattern',
        default='**/NC_SI_and_SYN_*_nc160_rule_+ANY_SFSs.txt',
        help='Glob pattern to match files (default: **/NC_SI_and_SYN_*_nc160_rule_+ANY_SFSs.txt)'
    )
    parser.add_argument(
        '--backup',
        action='store_true',
        help='Create backup files (default: False)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print detailed progress (default: False)'
    )

    return parser.parse_args()


def process_sfs_files(root_dir: str, pattern: str, backup: bool = False, verbose: bool = False) -> None:
    """
    Find all files matching pattern and process them using reformat_sfs_file.py
    
    Args:
        root_dir: Root directory to search for files
        pattern: Glob pattern to match files
        backup: Whether to create backup files
        verbose: Whether to print detailed progress
    """
    # Find all matching files
    matching_files = glob.glob(os.path.join(root_dir, pattern), recursive=True)
    
    if verbose:
        print(f"Found {len(matching_files)} files to process")
    
    # Process each file
    for i, filepath in enumerate(matching_files, 1):
        if verbose:
            print(f"\nProcessing file {i}/{len(matching_files)}: {filepath}")
        
        cmd = ['python', 'reformat_sfs_file.py', '-f', filepath]
        if backup:
            cmd.append('--backup')
            
        # Run the command
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            if verbose:
                print(result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error processing {filepath}: {e}")
            if verbose:
                print(f"Error output: {e.stderr}")
            continue
            
    print(f"\nTotal files processed: {len(matching_files)}")

def main():
    """
    Main function
    """
    # Parse command line arguments
    args = parse_arguments()

    # Process files
    process_sfs_files(
        root_dir=args.directory,
        pattern=args.pattern,
        backup=args.backup,
        verbose=args.verbose
    )

if __name__ == "__main__":
    main()
