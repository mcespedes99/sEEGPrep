from pathlib import Path
import json
import logging


# Import cleanSEEG
from clean_seeg import cleanSEEG

def main():
    edf_path = snakemake.input.edf
    electrodes_tsv = snakemake.input.electrodes_tsv
    processes = snakemake.threads
    out_edf = snakemake.output.out_edf
    out_tsv = snakemake.output.out_tsv
    report_df = snakemake.output.report_df
    report_json = snakemake.output.report_json
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        # Call class
        seegTF = cleanSEEG(edf_path, # Using downsampled edf
                        processes = processes)

        # Apply rereferencing
        df_reref, report_reref = seegTF.rereference(electrodes_tsv,
                                                    out_edf_path = out_edf,
                                                    out_tsv_path = out_tsv, 
                                                    return_report=True)
        # Save
        df_reref.to_csv(report_df, index=False, sep="\t")
        with open(report_json, 'w') as json_file:
            json.dump(report_reref, json_file)
    except:
        logging.exception('Got exception on main handler')
        raise

if __name__=="__main__":
    main()