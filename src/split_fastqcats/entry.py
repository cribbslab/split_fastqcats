#!/usr/bin/env python3

'''
Split-fastqcats.py - pipeline to split concatenated long_read cDNA reads
=========================================================================

:Tags: ONT, long-read

To use a specific workflow, type::

    split-fastqcats <workflow> [workflow options] [workflow arguments]

For this message and a list of available keywords type::

    split-fastqcats --help

To get help for a specific workflow, type::

    split-fastqcats <workflow> --help
'''

import os
import sys
import re
import glob
import importlib.util
import split_fastqcats


def print_list_in_columns(l, ncolumns):
    if not l:
        return ""

    max_width = max(len(x) for x in l) + 3
    n = len(l) // ncolumns + (len(l) % ncolumns > 0)
    columns = [l[i * n:(i + 1) * n] for i in range(ncolumns)]

    for i in range(len(columns)):
        while len(columns[i]) < n:
            columns[i].append("")

    rows = zip(*columns)
    pattern = " ".join([f"%-{max_width}s"] * ncolumns)
    return "\n".join([pattern % row for row in rows])


def dynamic_import_and_run(filepath, modulename, args):
    spec = importlib.util.spec_from_file_location(modulename, filepath)
    if spec is None or spec.loader is None:
        print(f"Could not load module from {filepath}")
        sys.exit(1)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    if hasattr(module, 'main'):
        module.main(args)
    else:
        print(f"Module {modulename} does not have a main() function.")
        sys.exit(1)


def main(argv=None):
    if argv is None:
        argv = sys.argv

    if len(argv) == 1 or argv[1] in ("--help", "-h"):
        # Show help and available pipelines
        package_path = os.path.dirname(os.path.abspath(split_fastqcats.__file__))
        pipelines = glob.glob(os.path.join(package_path, "pipeline_*.py"))

        print("Usage: split-fastqcats <command>\n")
        print("Available pipeline commands:\n")
        print(print_list_in_columns(
            sorted([os.path.basename(p)[len("pipeline_"):-len(".py")] for p in pipelines]), 2
        ))
        print("\nSpecial helper commands:\n")
        print(print_list_in_columns(["primer_pair_split", "barcode_split"], 2))
        return

    command = re.sub("-", "_", argv[1])
    args = argv[2:]

    # Base path of the installed package
    base_path = os.path.dirname(os.path.abspath(split_fastqcats.__file__))

    # Handle python fastq-splitters directly
    if command == "primer_pair_split":
        helper_path = os.path.join(base_path, "python", "fastq_splitter_fuzzy.py")
        dynamic_import_and_run(helper_path, "fastq_splitter_fuzzy", args)
        return

    elif command == "barcode_split":
        helper_path = os.path.join(base_path, "python", "fastq_splitter_index_fuzzy.py")
        dynamic_import_and_run(helper_path, "fastq_splitter_index_fuzzy", args)
        return

    # Handle pipeline commands
    pipeline_file = f"pipeline_{command}.py"
    pipeline_path = os.path.join(base_path, pipeline_file)

    if os.path.exists(pipeline_path):
        dynamic_import_and_run(pipeline_path, f"pipeline_{command}", args)
    else:
        print(f"Error: command '{command}' not found.")
        print("Use '--help' to list available commands.")
        sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
