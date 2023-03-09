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
    parc_path = str(snakemake.input.parc)
    noncon_to_con_tf_path = snakemake.input.tf
    processes = int(snakemake.config['processes'])
    out_edf = snakemake.output.out_edf
    out_tsv = snakemake.output.out_tsv

    # Call class
    seegTF = cleanSEEG(edf_path, # Using downsampled edf
                   chn_tsv_path,
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

    # Identify regions
    df = seegTF.identify_regions(parc_path,
                             use_reref = False,
                             write_tsv = True,
                             out_tsv_path = out_tsv,
                             df_cols = None, ## Using default as it was written in previous step
                             use_clean = False,
                             discard_wm_un = True, # TODO: put as argument
                             write_edf = True,
                             out_edf_path = out_edf)

if __name__=="__main__":
    main()