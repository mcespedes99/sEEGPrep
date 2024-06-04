from pathlib import Path
import sys
import logging

# Import cleanSEEG
from clean_seeg import cleanSEEG

def main():
    edf_path = snakemake.input.edf
    chn_csv_path = snakemake.input.chn_csv
    processes = snakemake.threads
    out_edf = snakemake.output.out_edf
    report = snakemake.output.report_file
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        # Call class
        seegTF = cleanSEEG(edf_path,
                        processes = processes)

        # Apply downsampling
        target_srate = snakemake.config['target_srate']
        _, dn_report = seegTF.downsample(chn_csv_path,
                                                  target_srate=target_srate,
                                                  out_edf_path = out_edf,
                                                  return_report=True)
        dn_report.to_csv(report, index=False, sep="\t")
    except:
        logging.exception('Got exception on main handler')
        raise

if __name__=="__main__":
    main()