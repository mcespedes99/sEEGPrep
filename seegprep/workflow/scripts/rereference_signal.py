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
                   epoch_length=5,
                   processes = processes)

    # Apply rereferencing
    df_cols = { # TODO: change to parameter!
            'type': 'type',
            'label': 'label',
            'x': 'x',
            'y': 'y',
            'z': 'z',
            'group': 'orig_group'
        }
    seegTF.rereference(out_edf, write_tsv = True, out_tsv_path = out_tsv, df_cols = df_cols)

if __name__=="__main__":
    main()