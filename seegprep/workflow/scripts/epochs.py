from pathlib import Path
import sys

# Adding path to import cleanSEEG
path = str(Path(Path(__file__).parent.absolute()).parent.parent.parent.absolute())
# print(path)
sys.path.append(path)

# Import cleanSEEG
from clean_seeg import cleanSEEG

def main():
    edf = snakemake.input.edf
    processes = int(snakemake.threads)
    out_edf = snakemake.output.out_edf
    tmpdir = snakemake.params.tmpdir

    # Call class
    seegTF = cleanSEEG(edf,
                       processes = processes)

    # Extract epochs
    seegTF.extract_epochs(event_label='awake trigger',
                          out_edf_path=out_edf,
                          tmpdir = tmpdir)

if __name__=="__main__":
    main()