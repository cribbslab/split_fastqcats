#!/usr/bin/env python3

'''
split-fastqcats - pipeline to split concatenated long_read cDNA reads
======================================================================

Usage:
    split-fastqcats <workflow> [workflow options] [workflow arguments]

Help:
    split-fastqcats --help
    split-fastqcats <workflow> --help
'''

import os
import sys
import re
import glob
import importlib.util
import runpy
import split_fastqcats
from split_fastqcats.version import __version__


# Map commands to helper script filenames
HELPER_COMMANDS = {
    "primer_pair_split": "python/fastq_splitter_fuzzy.py",
    "barcode_split": "python/fastq_splitter_index_fuzzy.py"
}


def print_list_in_columns(l, ncolumns=2):
    if not l:
        return
    max_width = max(len(x) for x in l) + 3
    rows = (len(l) + ncolumns - 1) // ncolumns
    columns = [l[i * rows:(i + 1) * rows] for i in range(ncolumns)]
    for col in columns:
        col.extend([""] * (rows - len(col)))
    return "\n".join(
        " ".join(f"{col[i]:<{max_width}}" for col in columns)
        for i in range(rows)
    )


def main(argv=None):
    if argv is None:
        argv = sys.argv
    
    if len(argv) >= 1 and argv[1] in ("--version"):
        print(f"split-fastqcats version {__version__}")
        return

    if len(argv) == 1 or argv[1] in ("--help", "-h"):
        print(__doc__)
        package_path = os.path.dirname(os.path.abspath(split_fastqcats.__file__))

        pipelines = sorted([
            os.path.basename(p)[len("pipeline_"):-len(".py")]
            for p in glob.glob(os.path.join(package_path, "pipeline_*.py"))
        ])

        print("Available pipeline commands:\n")
        print(print_list_in_columns(pipelines))
        print("\nSpecial helper commands:\n")
        print(print_list_in_columns(list(HELPER_COMMANDS.keys())))
        return

    command = re.sub("-", "_", argv[1])
    args = argv[2:]

    package_path = os.path.dirname(os.path.abspath(split_fastqcats.__file__))
    pipeline_path = os.path.join(package_path, f"pipeline_{command}.py")

    # Run CGAT-style pipeline
    if os.path.exists(pipeline_path):
        sys.argv = [f"pipeline_{command}"] + args
        runpy.run_path(pipeline_path, run_name="__main__")
        return

    # Run helper script via mapped filename
    elif command in HELPER_COMMANDS:
        helper_file = HELPER_COMMANDS[command]
        helper_path = os.path.join(package_path, helper_file)
        if not os.path.exists(helper_path):
            print(f"Error: helper script '{helper_file}' not found.")
            sys.exit(1)
        sys.argv = [helper_path] + args
        runpy.run_path(helper_path, run_name="__main__")
        return

    # Command not found
    print(f"Error: command '{command}' not recognized.")
    sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())