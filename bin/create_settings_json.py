#!/usr/bin/env python3
"""
Creates settings json file for the DMS_ABC script.
"""
import argparse
import json


def save_settings(filepath, settings):
    """Save settings to a JSON file."""
    with open(filepath, 'w') as file:
        json.dump(settings, file, indent=4)


def get_readingframes(reference_sequence, frameshift_position, frameshift_offset):
    """Returns the two reading frames"""
    frame1 = reference_sequence[:frameshift_position]
    frame2 = reference_sequence[frameshift_position - frameshift_offset:]
    return frame1, frame2


def main():
    parser = argparse.ArgumentParser(description="Create settings for DMS_ABC script.")
    parser.add_argument("--jsonpath", default="settings.json")
    parser.add_argument("--data_files")
    parser.add_argument("--job_name")
    parser.add_argument("--input_dir")
    parser.add_argument("--output_dir")
    parser.add_argument("--position_list", nargs='+', type=int)  # list codons to be analyzed
    parser.add_argument("--reference_sequence")
    parser.add_argument("--readingframes", action='store_true')  # bool multiple reading frames
    parser.add_argument("--frameshift_position", default=0, type=int)
    parser.add_argument("--frameshift_offset", default=0, type=int)
    args = parser.parse_args()

    settings = {
        'data_files': [args.data_files],
        'job_name': args.job_name,
        'input_dir': args.input_dir,
        'output_dir': args.output_dir,
        'position_list': args.position_list,
        'reference_sequence': args.reference_sequence,
        'readingframes': args.readingframes,
        'frameshift_position': args.frameshift_position,
        'frameshift_offset': args.frameshift_offset
    }

    if bool(args.readingframes) == True:
        position = int(args.frameshift_position)
        offset = int(args.frameshift_offset)
        frame1, frame2 = get_readingframes(args.reference_sequence, position, offset)
        settings['frameshift_position'] = position
        settings['frameshift_offset'] = offset
        settings['frame1'] = frame1
        settings['frame2'] = frame2

    save_settings(args.jsonpath, settings)


if __name__ == "__main__":
    main()
