from pathlib import Path
import sys

# Adding path to import cleanSEEG
path = str(Path(Path(__file__).parent.absolute()).parent.parent.parent.absolute())
# print(path)
sys.path.append(path)

# Import cleanSEEG
from clean_seeg import cleanSEEG

def main():
    edf_path = snakemake.input.edf
    chn_tsv_path = snakemake.input.tsv
    processes = int(snakemake.config['processes'])
    out_edf = snakemake.output.out_edf
    out_tsv = snakemake.output.out_tsv
    noncon_to_con_tf_path = snakemake.input.tf
    subject = snakemake.params.freesurf_patient
    subjects_dir = snakemake.params.freesurf_dir

    # Call class
    seegTF = cleanSEEG(edf_path, # Using downsampled edf
                   chn_tsv_path, 
                   subject, 
                   subjects_dir, 
                   cleanPLI = True, 
                   RmTrendMethod = 'LinearDetrend',
                   methodPLI = 'NotchFilter', 
                   lineFreq = 60,
                   bandwidth = 4,
                   n_harmonics = 1, # Only removing fundamental freq
                   noiseDetect = True,
                   # highpass = None, 
                   maxFlatlineDuration = 5, 
                   trsfPath=noncon_to_con_tf_path, # This is the only one I'm changing from default 
                   epoch_length=5,
                   processes = processes)

    # Apply filters
    clean, df_epochs = seegTF.clean_epochs(return_interpolated=False, 
                                                            write_edf_clean = True,
                                                            out_edf_path_clean = out_edf,
                                                            write_tsv = True,
                                                            out_tsv_path = out_tsv,
                                                            verbose = False
                                                        )

if __name__=="__main__":
    main()