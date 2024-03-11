from pathlib import Path
import json
import logging

# Import cleanSEEG
from clean_seeg import cleanSEEG


def main():
    edf_path = snakemake.input.edf
    electrodes_tsv = snakemake.input.electrodes_tsv
    parc_path = snakemake.input.parc
    colortable_file = snakemake.input.colortable
    processes = snakemake.threads
    reference = snakemake.params.reference_edf
    out_edf = snakemake.output.out_edf
    out_tsv = snakemake.output.out_tsv
    out_json = snakemake.params.out_json
    out_mask = snakemake.params.out_mask
    report_df = snakemake.output.report_df
    report_json = snakemake.output.report_json
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        # Call class
        seegTF = cleanSEEG(
            edf_path,  # Using downsampled edf
            processes=processes,
        )

        # Identify regions
        _, df_report_filt, report = seegTF.identify_regions(
            electrodes_tsv,
            parc_path,
            colortable_file,
            reference,
            use_reref = False,
            write_tsv = True,
            out_tsv_path = out_tsv,
            discard_wm_un = snakemake.config["discard_wm_un"],
            write_edf = True,
            out_edf_path = out_edf,
            vol_version = snakemake.config["vol_version"],
            json_out = out_json,
            mask_out = out_mask,
            return_report=True
        )

        # Save
        # Save
        df_report_filt.to_csv(report_df, index=False, sep="\t")
        with open(report_json, 'w') as json_file:
            json.dump(report, json_file)
    except:
        logging.exception("Got exception on main handler")
        raise


if __name__ == "__main__":
    main()
