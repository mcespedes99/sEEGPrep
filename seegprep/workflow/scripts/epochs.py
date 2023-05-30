from pathlib import Path
import sys
import logging
import os

# Adding path to import cleanSEEG
path = str(Path(Path(__file__).parent.absolute()).parent.parent.parent.absolute())
# print(path)
sys.path.append(path)

# Import cleanSEEG
from clean_seeg import cleanSEEG

def main():
    edf = snakemake.input.edf
    processes = int(snakemake.threads)
    tmp_file = str(snakemake.output.tmp_file)
    tmpdir = snakemake.resources.tmpdir
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        # Create dummy tmp file
        Path(tmp_file).touch()
        # Update out_dir
        path = Path(tmp_file)
        out_dir = path.parent.absolute()

        # Call class
        seegTF = cleanSEEG(edf,
                           processes = processes)

        # Extract epochs
        seegTF.extract_epochs(event_label='awake trigger',
                              out_root=out_dir,
                              tmpdir = tmpdir,
                              snakemake = True)
    
    except:
        logging.exception('Got exception on main handler')
        raise

if __name__=="__main__":
    main()