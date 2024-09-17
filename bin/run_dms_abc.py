#!/usr/bin/env python3
"""Running Gianni's adapted DMS_ABC pipeline"""
import argparse
import json
import os

import dnaio

import codon_truncation_bam_overlap
import counter
#import create_count_file
import create_count_file_stand_alone_onlymutcod_wt_spike_efref_NNK


def parse_arguments():
    """Parse the command line arguments"""
    parser = argparse.ArgumentParser(description="Settings for DMS_ABC script.")
    parser.add_argument("--bam")
    parser.add_argument("--reference")  # reference fasta file
    parser.add_argument("--positions", nargs='+', type=int)  # list codons to be analyzed
    parser.add_argument("--nnk_positions", nargs='+', type=int)
    parser.add_argument("--wt_ref_position")
    parser.add_argument("--wt_codon")
    parser.add_argument("--wt_count", default=10000, type=int)
    parser.add_argument("--readingframes", action='store_true')  # bool multiple reading frames
    parser.add_argument("--frameshift_position", default=0, type=int)
    parser.add_argument("--frameshift_offset", default=0, type=int)
    return parser.parse_args()


def get_reading_frames(reference_sequence, frameshift_position, frameshift_offset):
    """Returns the two reading frames"""
    frame1 = reference_sequence[:frameshift_position]
    frame2 = reference_sequence[frameshift_position - frameshift_offset:]
    return frame1, frame2


def read_refrence(fasta_file):
    """Returns the name and the sequence of a fasta file"""
    with dnaio.open(fasta_file) as file:
        for record in file:
            name = record.name
            sequence = record.sequence
            break
    return name, sequence


def DMS_processing(parameters):
    """Running the DMS_ABC scripts"""
    input_file = parameters["bam"]
    codontruncated_file = input_file[:-4] + "_codontruncated.bam"
    triplet_count_file = input_file[:-4] + "_triplet_count.txt"

    frameshift_position = parameters["frameshift_position"]
    frameshift_offset = parameters["frameshift_offset"]
    reference_name = parameters["reference_name"]
    reference_sequence = parameters["reference_sequence"]
    positions = parameters["positions"]
    nnk_positions = parameters["nnk_positions"]
    wt_ref_position = parameters["wt_ref_position"]
    wt_codon = parameters["wt_codon"]
    wt_count = parameters["wt_count"]

    codon_truncation_bam_overlap.codon_truncation(input_file, codontruncated_file, frameshift_position, frameshift_offset, reference_name)

    counter.count_mutants(codontruncated_file, triplet_count_file, positions, reference_name, reference_sequence)

    #create_count_file.make_HDF5(triplet_count_file, reference_sequence, frameshift_position, frameshift_offset)

    create_count_file_stand_alone_onlymutcod_wt_spike_efref_NNK.make_HDF5(triplet_count_file, reference_sequence, frameshift_position, frameshift_offset, nnk_positions, wt_ref_position, wt_codon, wt_count)


def main():
    args = parse_arguments()

    reference_name, reference_sequence = read_refrence(args.reference)

    parameters = {
        'bam': args.bam,
        'reference_name': reference_name,
        'reference_sequence': reference_sequence,
        'positions': args.positions,
        'nnk_positions': args.nnk_positions,
        'wt_ref_position': args.wt_ref_position,
        'wt_codon': args.wt_codon,
        'wt_count': args.wt_count,
        'readingframes': args.readingframes,
        'frameshift_position': args.frameshift_position,
        'frameshift_offset': args.frameshift_offset
    }

    if args.readingframes:
        frame1, frame2 = get_reading_frames(reference_sequence, args.frameshift_position, args.frameshift_offset)
        parameters['frame1'] = frame1
        parameters['frame2'] = frame2

    DMS_processing(parameters)


if __name__ == "__main__":
    main()
