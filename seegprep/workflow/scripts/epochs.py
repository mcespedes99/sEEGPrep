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
    event = snakemake.config['event']
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        # Create dummy tmp file
        Path(tmp_file).touch()
        # Update out_dir
        path = Path(tmp_file)
        out_dir = path.parent.absolute()

        # Default
        n_samples = None
        duration = None
        # Get event information
        event_start, type_length, length = event
        assert type_length in ['n_samples', 'duration']
        if type_length == 'n_samples':
            n_samples = length
        else:
            duration = length

        # Call class
        seegTF = cleanSEEG(edf,
                           processes = processes)

        # Extract epochs
        # Extract epoch
        seegTF.extract_epochs(
            n_samples = n_samples,
            event_dur = duration,
            event_start=event_start,
            out_root=out_dir,
            tmpdir=tmpdir,
        )
    
    except:
        logging.exception('Got exception on main handler')
        raise

if __name__=="__main__":
    main()