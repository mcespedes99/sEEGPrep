from pathlib import Path
import sys
import logging
import json

# Import cleanSEEG
from clean_seeg import cleanSEEG


def main():
    edf_path = snakemake.input.edf
    chn_tsv_path = snakemake.input.chn_tsv
    # t1_path = snakemake.input.t1
    processes = snakemake.threads
    out_edf = snakemake.output.out_edf
    report_df = snakemake.output.report_df
    report_json = snakemake.output.report_json
    # subject = snakemake.params.freesurf_patient
    # subjects_dir = snakemake.params.freesurf_dir
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        # Call class
        seegTF = cleanSEEG(
            edf_path,
            RmTrendMethod=snakemake.config["detrend_method"],
            highpass=snakemake.config["highpass"],
            processes=processes,
        )

        # Apply filters
        _, df_report_filt, report_filt = seegTF.drift_correction(
            chn_tsv_path,
            out_edf_path_clean=out_edf,
        )
        # Save
        df_report_filt.to_csv(report_df, index=False, sep="\t")
        with open(report_json, 'w') as json_file:
            json.dump(report_filt, json_file)
    except:
        logging.exception("Got exception on main handler")
        raise


if __name__ == "__main__":
    main()
