from pathlib import Path
import sys
import logging

# Adding path to import cleanSEEG
path = str(Path(Path(__file__).parent.absolute()).parent.parent.parent.absolute())
# print(path)
sys.path.append(path)

# Import cleanSEEG
from clean_seeg import cleanSEEG

def main():
    edf_path = snakemake.input.edf
    chn_tsv_path = snakemake.input.tsv
    processes = snakemake.config['processes']
    out_edf = snakemake.output.out_edf
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
      # Call class
      seegTF = cleanSEEG(edf_path, 
                    chn_tsv_path,
                    methodPLI = snakemake.config['methodPLI'], 
                    lineFreq = snakemake.config['lineFreq'],
                    bandwidth = snakemake.config['bandwidth'],
                    n_harmonics = snakemake.config['n_harmonics'],
                    processes = processes)

      # Apply filters
      clean_signal = seegTF.reject_PLI(write_edf=True,
                                       out_edf_path=out_edf)
    except:
        logging.exception('Got exception on main handler')
        raise

if __name__=="__main__":
    main()