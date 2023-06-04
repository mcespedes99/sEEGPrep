from pathlib import Path
import sys
import logging
import os
import re

# Adding path to import cleanSEEG
path = str(Path(Path(__file__).parent.absolute()).parent.parent.parent.absolute())
# print(path)
sys.path.append(path)

def rename_files(rule, tsv_files, out_tsv):
    reg = re.compile(rule)
    tsv_files_tmp = list(filter(reg.search, tsv_files))
    out_tsv_tmp = list(filter(reg.search, out_tsv))
    if len(out_tsv_tmp)>0:
        # Rename the first file
        os.rename(tsv_files_tmp[0], out_tsv_tmp[0])
        # Delete the rest
        for file in tsv_files_tmp[1:]:
            os.remove(file)

def main():
    tsv_files  = snakemake.input.tsv_files
    out_tsv = snakemake.output.out_files
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        # Convert to list if str is provided (not sure if needed)
        if type(out_tsv) == str:
            out_tsv = [out_tsv]
        else: # It's a snakemake.io.Namedlist
            out_tsv = list(out_tsv)
        print(out_tsv)
        # Manage possible scenarios by separating tsv files into 2
        # First case of region ID
        print(out_tsv)
        print(type(out_tsv))
        rename_files('regionID', tsv_files, out_tsv)
        # Case of reref
        rename_files('reref', tsv_files, out_tsv)
    except:
        logging.exception('Got exception on main handler')
        raise    

if __name__=="__main__":
    main()