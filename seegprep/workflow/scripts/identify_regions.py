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
    edf_path, chn_tsv_path = snakemake.input.edf_tsv
    parc_list = snakemake.input.parc_list
    colortable_file = snakemake.params.colortable
    tfm = snakemake.input.tfm
    processes = snakemake.config["processes"]
    out_edf = snakemake.output.out_edf
    out_tsv = snakemake.output.out_tsv
    masks_out = snakemake.output.out_masks
    colormask_out = snakemake.output.out_colormask
    out_json = snakemake.output.out_json
    reref_run = snakemake.params.reref_run
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        # Call class
        seegTF = cleanSEEG(
            edf_path,  # Using downsampled edf
            chn_tsv_path,
            tfm=[(tfm, False)],
            processes=processes,
        )
        if not reref_run:
            print("reref not run before regions_id")
            # Define names of columns in tsv file
            dict_keys = ["type", "label", "x", "y", "z", "group"]
            dict_vals = snakemake.config["tsv_cols"]
            df_cols = dict(zip(dict_keys, dict_vals))
        else:
            df_cols = None
        # Identify regions
        df = seegTF.identify_regions(
            parc_list,
            colortable_file,
            use_reref=False,
            write_tsv=True,
            out_tsv_path=out_tsv,
            df_cols=df_cols,  ## Using default as it was written in previous step
            use_clean=False,
            discard_wm_un=snakemake.config["discard_wm_un"],
            write_edf=True,
            out_edf_path=out_edf,
            vol_version=True,
            masks_out=masks_out,
            colormask_out=colormask_out,
            json_out=out_json,
        )
    except:
        logging.exception("Got exception on main handler")
        raise


if __name__ == "__main__":
    main()
