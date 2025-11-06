"""
Script to format output SFS files from SFRatios/drosophila_work/get_short_intron_paired_SNP_allele_counts.py
to make it ready for SFRatios/SFRatios.py.

The output produced by get_short_intron_paired_SNP_allele_counts.py has:
    line 1: Neutral SFS header
    line 2: Neutral SFS with value != 0 for the 0-bin
    line 3: Selected SFS header
    line 4: Selected SFS with value != 0 for the 0-bin

SFRatios.py require a file containing two SFSs as follow:
data text file format: 
    line 1:  arbitrary text (usually a header)
    line 2:  neutral SFS beginning with 0 bin (which is ignored)
    line 3:  blank
    line 4:  selected SFS beginning with 0 bin (which is ignored)
"""

import argparse
import os
import shutil


def parse_arguments():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(description='Filter SNPeff-consistency output based on mode and effects')

    parser.add_argument('-f', '--inputfile', help='Input table from SNPeff-consistency')
    parser.add_argument('--backup', action='store_true', help='Make a backup of the inputfile (default: False)')
    args = parser.parse_args()

    return args


def make_backup_filename(inputfile: str) -> str:
    """
    Take the inputfile and attach 'bckup' at the end of the filename.
    With shutil.copy, make a backup of inputfile in the same path as the
    original file.
    """

    # Process original filepath and filename
    filepath, ffile = os.path.split(inputfile)
    filename, fextension = ffile.split(".")

    # Make filename for the backup file
    if len(filepath) > 0:
        inputfilebckp = f"{filepath}/{filename}_bckp.{fextension}"
    else:
        inputfilebckp = f"{filename}_bckp.{fextension}"

    return inputfilebckp


def reformat_sfs_file(inputfile: str) -> None:
    """
    Reads a SFS file, create a new header and save the re-formatted file in
    the same inputfile.
    """

    # Read the input file
    with open(inputfile, 'r', encoding="utf-8") as f:
        lines = f.readlines()

    # Make new header
    nsfs_header = lines[0].strip().split("_")
    ssfs_header = lines[2].strip().split("_")

    new_header_prefix = "_".join([nsfs_header[0], ssfs_header[0]])
    new_header_suffix = "_".join(nsfs_header[1:])
    new_header = "_".join([new_header_prefix, new_header_suffix])

    # Process the lines
    new_lines = []
    new_lines.append(f"{new_header}\n")

    # Process first data line (replace first number with 0)
    neutral_data = lines[1].strip().split()
    neutral_data[0] = "0"
    new_lines.append(" ".join(neutral_data) + "\n\n")

    # Process second data line (replace first number with 0)
    selection_data = lines[3].strip().split()
    selection_data[0] = "0"
    new_lines.append(" ".join(selection_data) + "\n")

    # Write back to the original file
    with open(inputfile, 'w', encoding="utf-8") as f:
        f.writelines(new_lines)

    print(f"Processed: {inputfile}")


def main():
    """
    Main function
    """

    args = parse_arguments()
    inputfile = args.inputfile

    if args.backup:
        inputfilebckp = make_backup_filename(inputfile)
        shutil.copy(inputfile, inputfilebckp)

    reformat_sfs_file(inputfile)

if __name__ == "__main__":
    main()
