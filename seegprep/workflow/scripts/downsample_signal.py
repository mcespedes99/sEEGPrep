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
    processes = snakemake.config['processes']
    out_edf = snakemake.output.out_edf
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        # Call class
        seegTF = cleanSEEG(edf_path,
                        processes = processes)

        # Apply downsampling
        target_srate = snakemake.config['target_srate']
        signal_dsG, downsampledSrate = seegTF.downsample(target_srate=target_srate,
                                                        write_edf = True,
                                                        out_edf_path = out_edf)
    except:
        logging.exception('Got exception on main handler')
        raise

if __name__=="__main__":
    main()