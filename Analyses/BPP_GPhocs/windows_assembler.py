#!/usr/bin/env python3

"""
windows_assembler.py

Author: Manuel Riveros Escalona
Created: 2026-03-12
Contact: manuelriveros72@hotmail.com

Description
-----------
Combines window FASTA files from multiple samples into locus-based alignments.

Each input FASTA file contains sequences corresponding to genomic windows
(e.g. generated with BEDTools getfasta). The script:

1. Combines homologous windows across samples
2. Filters loci by maximum missing data threshold
3. Writes locus alignments in PHYLIP format
4. Generates a combined G-PhoCS input file
5. Organizes PHYLIP loci into subdirectories

Usage
-----
python windows_assembler.py -i window_list.txt [-n 2000] [-m 0.2]

Arguments
---------
-i / --input_file_list
    File containing paths to window FASTA files (one per line)

-n / --files_per_subdir
    Maximum number of loci per output directory (default: 2000)

-m / --max_missing
    Maximum fraction of missing bases allowed per locus (default: 0.2)

Input
-----
FASTA files containing window sequences for each sample.

Output
------
combined_loci/
    gphocs_input.txt     Combined input file for G-PhoCS
    block_1/
        locus1.phy
        locus2.phy
        ...
    block_2/
        ...

Dependencies
------------
Python >= 3.8
Biopython

License
-------
Free for academic use.
"""

import os
import argparse
import re
from Bio import SeqIO
from collections import defaultdict


def combine_sequences(input_files):

    combined_sequences = defaultdict(dict)
    locus_lengths = {}
    expected_loci = None

    for input_file in input_files:

        sample = os.path.basename(input_file).split(".")[0]

        sample_loci = []

        for record in SeqIO.parse(input_file, "fasta"):

            locus = re.sub(r"[^A-Za-z0-9_.]","_", record.id).lstrip("_")
            seq = str(record.seq).upper()

            sample_loci.append(locus)

            if locus not in locus_lengths:
                locus_lengths[locus] = len(seq)
            else:
                if len(seq) != locus_lengths[locus]:
                    raise ValueError(
                        f"Sequence length mismatch at {locus} for {sample}"
                    )

            combined_sequences[locus][sample] = seq

        if expected_loci is None:
            expected_loci = set(sample_loci)
        else:
            if set(sample_loci) != expected_loci:
                raise ValueError(
                    f"Locus mismatch detected in sample {sample}"
                )

    return combined_sequences, locus_lengths


def write_outputs(combined_sequences, locus_lengths, output_directory,
                  files_per_subdir, max_missing):

    os.makedirs(output_directory, exist_ok=True)

    gphocs_file = open(os.path.join(output_directory, "gphocs_input.txt"), "w")

    subdir_count = 1
    subdir_index = 0
    locus_index = 0

    subdir_path = os.path.join(output_directory, f"block_{subdir_count}")
    os.makedirs(subdir_path, exist_ok=True)

    for locus, species_sequences in combined_sequences.items():

        seqs = species_sequences
        samples = sorted(seqs.keys())
        length = locus_lengths[locus]

        # missing data filtering
        total_missing = sum(seqs[s].count("N") for s in samples)
        missing_fraction = total_missing / (length * len(samples))

        if missing_fraction > max_missing:
            continue

        # write G-PhoCS
        gphocs_file.write(f"{locus} {len(samples)} {length}\n")
        for s in samples:
            gphocs_file.write(f"{s} {seqs[s]}\n")
        gphocs_file.write("\n")

        # directory splitting
        if subdir_index == files_per_subdir:
            subdir_count += 1
            subdir_index = 0
            subdir_path = os.path.join(output_directory, f"block_{subdir_count}")
            os.makedirs(subdir_path, exist_ok=True)

        phylip_file = os.path.join(subdir_path, f"{locus}.phy")

        with open(phylip_file, "w") as out:

            out.write(f"{len(samples)} {length}\n")

            for s in samples:
                out.write(f"^{s} {seqs[s]}\n")

        subdir_index += 1
        locus_index += 1

        if locus_index % 1000 == 0:
            print(f"{locus_index} loci processed")

    gphocs_file.close()


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input_file_list",
        required=True,
        help="file containing list of window fasta files"
    )

    parser.add_argument(
        "-n",
        "--files_per_subdir",
        type=int,
        default=2000,
        help="max loci per folder"
    )

    parser.add_argument(
        "-m",
        "--max_missing",
        type=float,
        default=0.2,
        help="maximum fraction of missing bases allowed"
    )

    args = parser.parse_args()

    with open(args.input_file_list) as f:
        input_files = [x.strip() for x in f]

    combined_sequences, locus_lengths = combine_sequences(input_files)

    write_outputs(
        combined_sequences,
        locus_lengths,
        "combined_loci",
        args.files_per_subdir,
        args.max_missing
    )


if __name__ == "__main__":
    main()
