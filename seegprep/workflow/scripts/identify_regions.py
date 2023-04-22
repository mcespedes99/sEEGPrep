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
    parc_path = snakemake.input.parc
    noncon_to_con_tf_path = snakemake.input.tf
    processes = snakemake.config['processes']
    out_edf = snakemake.output.out_edf
    out_tsv = snakemake.output.out_tsv
    reref_run = snakemake.params.reref_run
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        # Call class
        seegTF = cleanSEEG(edf_path, # Using downsampled edf
                        chn_tsv_path,
                        trsfPath = noncon_to_con_tf_path,
                        processes = processes)
        if not reref_run:
            print('reref not run before regions_id')
            # Define names of columns in tsv file
            dict_keys = ['type','label','x','y','z','group']
            dict_vals = snakemake.config['tsv_cols']
            df_cols = dict(zip(dict_keys, dict_vals))
            # df_cols = { # TODO: change to parameter in config file!
            #         'type': 'type',
            #         'label': 'label',
            #         'x': 'x',
            #         'y': 'y',
            #         'z': 'z',
            #         'group': 'orig_group'
            #     }
        else:
            df_cols = None
        # Identify regions
        df = seegTF.identify_regions(parc_path,
                                use_reref = False,
                                write_tsv = True,
                                out_tsv_path = out_tsv,
                                df_cols = df_cols, ## Using default as it was written in previous step
                                use_clean = False,
                                discard_wm_un = snakemake.config['discard_wm_un'],
                                write_edf = True,
                                out_edf_path = out_edf)
    except:
        logging.exception('Got exception on main handler')
        raise

if __name__=="__main__":
    main()