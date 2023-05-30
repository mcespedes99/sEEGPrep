from pathlib import Path
import sys
import logging
import os

# Adding path to import cleanSEEG
path = str(Path(Path(__file__).parent.absolute()).parent.parent.parent.absolute())
# print(path)
sys.path.append(path)

def main():
    out_edf = str(snakemake.output.out_edf)
    tmp_edf = str(Path(out_edf).with_suffix(''))+'_tmp.edf'
    # Rename 
    os.rename(tmp_edf, out_edf)

if __name__=="__main__":
    main()