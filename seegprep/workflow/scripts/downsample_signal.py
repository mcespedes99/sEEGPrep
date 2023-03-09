from pathlib import Path
import sys

# Adding path to import cleanSEEG
path = str(Path(Path(__file__).parent.absolute()).parent.parent.parent.absolute())
# print(path)
sys.path.append(path)

# Import cleanSEEG
from clean_seeg import cleanSEEG

def main():
    edf_path = str(snakemake.input.edf)
    processes = int(snakemake.config['processes'])
    out_edf = str(snakemake.output.out_edf)

    # print(edf_path)
    # print(chn_tsv_path)
    # print(subject)
    # print(subjects_dir)
    # print(noncon_to_con_tf_path)
    # print(out_edf)
    # Call class
    seegTF = cleanSEEG(edf_path, 
                    cleanPLI = True, 
                    RmTrendMethod = 'LinearDetrend',
                    methodPLI = 'Zapline', 
                    lineFreq = 60,
                    bandwidth = 4,
                    noiseDetect = True,
                    highpass = [0.5, 1.5], #I set it to [0.5, 1.5] to improve comp cost
                    maxFlatlineDuration = 5, 
                    epoch_length=5,
                    processes = processes
                    )

    # Apply downsampling
    signal_dsG, downsampledSrate = seegTF.downsample(target_srate=200, write_edf = True, out_edf_path = out_edf)

if __name__=="__main__":
    main()