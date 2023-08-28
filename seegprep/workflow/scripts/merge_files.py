from pathlib import Path
import sys
import logging
import os
import re


def rename_files(rule, in_files, out_files):
    reg = re.compile(rule)
    in_files_tmp = list(filter(reg.search, in_files))
    out_files_tmp = list(filter(reg.search, out_files))
    if len(out_files_tmp) > 0:
        # Rename the first file
        os.rename(in_files_tmp[0], out_files_tmp[0])
        # Delete the rest
        for file in in_files_tmp[1:]:
            os.remove(file)


def main():
    in_files = snakemake.input.in_files
    out_files = snakemake.output.out_files
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        # Convert to list if str is provided (not sure if needed)
        if type(out_files) == str:
            out_files = [out_files]
        else:  # It's a snakemake.io.Namedlist
            out_files = list(out_files)
        print(out_files)
        # Manage possible scenarios by separating tsv files into 2
        # First case of region ID
        print(out_files)
        print(type(out_files))
        rename_files("regionID", in_files, out_files)
        # Case of reref
        rename_files("reref", in_files, out_files)
    except:
        logging.exception("Got exception on main handler")
        raise


if __name__ == "__main__":
    main()
