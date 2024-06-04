from pathlib import Path
import json
import logging

# Import cleanSEEG
from clean_seeg import cleanSEEG

def main():
    edf_path, chn_tsv_path = snakemake.input.edf_tsv
    processes = snakemake.threads
    out_edf = snakemake.output.out_edf
    report_json = snakemake.output.report_json
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
      # Call class
      seegTF = cleanSEEG(edf_path, 
                    methodPLI = snakemake.config['methodPLI'], 
                    lineFreq = snakemake.config['lineFreq'],
                    bandwidth = snakemake.config['bandwidth'],
                    n_harmonics = snakemake.config['n_harmonics'],
                    processes = processes)

      # Apply filters
      _, report_PLI = seegTF.reject_PLI(chn_tsv_path,
                                       out_edf_path = out_edf,
                                       return_report = True)
      with open(report_json, 'w') as json_file:
         json.dump(report_PLI, json_file)
    except:
        logging.exception('Got exception on main handler')
        raise

if __name__=="__main__":
    main()