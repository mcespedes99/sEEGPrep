from pathlib import Path
import sys
import logging

print("init")

# Adding path to import cleanSEEG
path = str(Path(Path(__file__).parent.absolute()).parent.parent.parent.absolute())
# print(path)
sys.path.append(path)

# Import cleanSEEG
from clean_seeg import cleanSEEG


def main():
    edf_path = snakemake.input.edf
    chn_tsv_path = snakemake.input.tsv
    processes = snakemake.config["processes"]
    out_edf = snakemake.output.out_edf
    out_tsv = snakemake.output.out_tsv
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        # Call class
        seegTF = cleanSEEG(
            edf_path, chn_tsv_path, processes=processes  # Using downsampled edf
        )

        # Define names of columns in tsv file
        dict_keys = ["type", "label", "x", "y", "z", "group"]
        dict_vals = snakemake.config["tsv_cols"]
        df_cols = dict(zip(dict_keys, dict_vals))
        # df_cols = { # TODO: change to parameter!
        #         'type': 'type',
        #         'label': 'label',
        #         'x': 'x',
        #         'y': 'y',
        #         'z': 'z',
        #         'group': 'orig_group'
        #     }
        # Apply rereferencing
        seegTF.rereference(
            out_edf, write_tsv=True, out_tsv_path=out_tsv, df_cols=df_cols
        )

    except:
        logging.exception("Got exception on main handler")
        raise


if __name__ == "__main__":
    main()
