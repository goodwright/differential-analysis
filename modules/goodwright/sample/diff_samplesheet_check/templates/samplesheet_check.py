#!/usr/bin/env python

import os
import sys
import errno
import argparse
import platform

MIN_COLS = 2
HEADER = ["sample_id", "condition"]


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def dump_versions(process_name):
    with open("versions.yml", "w") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")


def check_samplesheet(process_name, samplesheet, counts, output):
    """
    This function checks that the samplesheet follows the following structure:

    sample_id,condition,factor(optional)
    sample_1,A,factor_1
    sample_2,A,factor_2
    sample_3,B,factor_1
    sample_4,B,factor_2
    sample_5,B,factor_1
    """
    # Dump version file
    dump_versions(process_name)

    # Init
    sample_dict = {}
    count_dict = {}

    with open(samplesheet, "r") as fin:

        # Check header
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        ACTUAL_HEADER_LEN = len(header)

        if ACTUAL_HEADER_LEN < MIN_COLS:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        if header[: MIN_COLS] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        # Check sample entries
        line_no = 1
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check if its just a blank line so we dont error
            if line.strip() == "":
                continue

            # Check valid number of columns per row
            if len(lspl) != ACTUAL_HEADER_LEN:
                print_error(
                    "Invalid number of columns (found {} should be {})! - line no. {}".format(
                        len(lspl), len(HEADER), line_no
                    ),
                    "Line",
                    line,
                )

            # Check valid number of populated columns per row
            num_cols = len([x for x in lspl if x])
            if num_cols < ACTUAL_HEADER_LEN:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(ACTUAL_HEADER_LEN),
                    "Line",
                    line,
                )

            # Check sample name entries
            factor = None
            if ACTUAL_HEADER_LEN > MIN_COLS:
                sample, condition, factor = lspl[: (MIN_COLS + 1)]
            else:
                sample, condition = lspl[: MIN_COLS]

            # Check for spaces
            if sample:
                if sample.find(" ") != -1:
                    print_error("Sample ID contains spaces!", "Line", line)
            if condition:
                if condition.find(" ") != -1:
                    print_error("Condition contains spaces!", "Line", line)
            if factor:
                if factor.find(" ") != -1:
                    print_error("Factor contains spaces!", "Line", line)

            # Check for dots
            if sample:
                if sample.find(".") != -1:
                    print_error("Sample ID contains dots!", "Line", line)
            if condition:
                if condition.find(".") != -1:
                    print_error("Condition contains dots!", "Line", line)
            if factor:
                if factor.find(".") != -1:
                    print_error("Factor contains dots!", "Line", line)

            # Set base sample info
            sample_info = [sample, condition]

            # Collect additional sample data
            extra_data = lspl[len(HEADER) :]
            sample_info = sample_info + extra_data

            # Create sample mapping dictionary
            if sample not in sample_dict:
                sample_dict[sample] = sample_info
            else:
                print_error("Samplesheet contains duplicate rows!", "Line", line)

            line_no = line_no + 1

    # Load sample ids from counts file
    with open(counts, "r") as fin:
        # Load header
        header = [x.strip('"') for x in fin.readline().strip().split(",")]

        # Strip first column
        header = header[1:]

        # Collect unique sample ids
        for sample in header:
            count_dict[sample] = True

    # Check that all samples in samplesheet are present in counts file
    for sample in sample_dict.keys():
        if sample not in count_dict.keys():
            print_error("Sample {} not found in counts file!".format(sample))

    # Check that all samples in counts file are present in samplesheet
    for sample in count_dict.keys():
        if sample not in sample_dict.keys():
            print_error("Sample {} in counts file not found in samplesheet!".format(sample))

    # Calculate output header
    ouput_header = ["sample_id", "condition"]
    extra_header = header[len(HEADER) :]
    ouput_header = ouput_header + extra_header

    # Write validated samplesheet with appropriate columns
    if len(sample_dict) > 0:
        out_dir = os.path.dirname(output)
        make_dir(out_dir)
        with open(output, "w") as fout:
            fout.write(",".join(ouput_header) + "\n")

            for sample in sorted(sample_dict.keys()):
                    fout.write(",".join(sample_dict[sample]) + "\n")


if __name__ == "__main__":

    # Allows switching between nextflow templating and standalone python running using arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    parser.add_argument("--samplesheet", default="!{samplesheet}")
    parser.add_argument("--counts", default="!{counts}")
    parser.add_argument("--output", default="!{output}")
    args = parser.parse_args()

    check_samplesheet(args.process_name, args.samplesheet, args.counts, args.output)
