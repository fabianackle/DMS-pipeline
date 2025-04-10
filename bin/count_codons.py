#!/usr/bin/env python3
import argparse
from collections import Counter
import concurrent.futures
import itertools
import subprocess

import dnaio
import polars as pl
import pysam


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id")
    parser.add_argument("--bam")
    parser.add_argument("--reference")
    parser.add_argument("--positions", nargs="+", type=int)
    return parser.parse_args()


def read_refrence(fasta_file):
    """Returns the name and the sequence of a fasta file"""
    with dnaio.open(fasta_file) as file:
        for record in file:
            name = record.name
            sequence = record.sequence
            break
    return name, sequence


def count_codons_for_pos(bam, ref_name, position):
    samfile = pysam.AlignmentFile(bam, "rb")
    start = position
    end = position + 3
    codons = []

    for read in samfile.fetch(ref_name, start, end):
        if read.is_unmapped:
            continue

        codon_bases = []
        codon_qualities = []
        for ref_pos in range(start, end):
            if ref_pos in read.get_reference_positions():
                read_pos = read.get_reference_positions(full_length=True).index(ref_pos)
            else:
                read_pos = -1

            if read_pos != -1:
                codon_bases.append(read.query_sequence[read_pos])
                codon_qualities.append(int(read.query_qualities[read_pos]))
            else:
                codon_bases.append("N")

        if ("N" not in codon_bases) and all(q >= 35 for q in codon_qualities):
            codon = "".join(codon_bases)
            codons.append(codon)

    samfile.close()

    counted_codons = {"start": start, "end": end}
    counted = Counter(codons)
    counted_codons.update(counted)

    return counted_codons


def main():
    args = parse_arguments()

    sample_id = args.sample_id
    bam = args.bam
    reference = args.reference
    positions = args.positions

    ref_name, _ = read_refrence(reference)

    subprocess.run(["samtools", "index", bam], check=True)

    with concurrent.futures.ProcessPoolExecutor() as executor:
        data = list(executor.map(
            count_codons_for_pos,
            itertools.repeat(bam),
            itertools.repeat(ref_name),
            positions
        ))

    pl.DataFrame(data).write_csv(f"{sample_id}.csv")


if __name__ == "__main__":
    main()
