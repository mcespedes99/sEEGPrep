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
    processes = snakemake.config['processes']
    out_edf = snakemake.output.out_edf

    # print(edf_path)
    # print(chn_tsv_path)
    # print(subject)
    # print(subjects_dir)
    # print(noncon_to_con_tf_path)
    # print(out_edf)
    # Call class
    seegTF = cleanSEEG(edf_path,
                       processes = processes)

    # Apply downsampling
    target_srate = snakemake.config['target_srate']
    signal_dsG, downsampledSrate = seegTF.downsample(target_srate=target_srate,
                                                     write_edf = True,
                                                     out_edf_path = out_edf)

if __name__=="__main__":
    main()