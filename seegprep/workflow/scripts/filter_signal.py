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
    out_tsv = snakemake.output.out_tsv
    noncon_to_con_tf_path = snakemake.input.tf
    subject = snakemake.params.freesurf_patient
    subjects_dir = snakemake.params.freesurf_dir
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
      # Call class
      seegTF = cleanSEEG(edf_path, 
                    chn_tsv_path,
                    RmTrendMethod = snakemake.config['detrend_method'],
                    methodPLI = snakemake.config['methodPLI'], 
                    lineFreq = snakemake.config['lineFreq'],
                    bandwidth = snakemake.config['bandwidth'],
                    n_harmonics = snakemake.config['n_harmonics'],
                    noiseDetect = snakemake.config['noiseDetect'],
                    highpass = snakemake.config['highpass'], 
                    maxFlatlineDuration = snakemake.config['maxFlatlineDuration'], 
                    trsfPath = noncon_to_con_tf_path, 
                    epoch_autoreject = snakemake.config['epoch_length'],
                    processes = processes)

      # Apply filters
      clean, df_epochs = seegTF.clean_epochs(subject = subject, 
                                              subjects_dir = subjects_dir,
                                              return_interpolated=False, 
                                              write_edf_clean = True,
                                              out_edf_path_clean = out_edf,
                                              write_tsv = True,
                                              out_tsv_path = out_tsv,
                                              verbose = False
                                            )
    except:
        logging.exception('Got exception on main handler')
        raise

if __name__=="__main__":
    main()