from pathlib import Path
import sys
import logging
import nibabel as nb

# Adding path to import cleanSEEG
path = str(Path(Path(__file__).parent.absolute()).parent.parent.parent.absolute())
# print(path)
sys.path.append(path)

# Import cleanSEEG
from clean_seeg import cleanSEEG

def main():
    edf_path = snakemake.input.edf
    chn_tsv_path = snakemake.input.tsv
    # t1_path = snakemake.input.t1
    processes = snakemake.config['processes']
    out_edf = snakemake.output.out_edf
    out_tsv = snakemake.params.out_tsv
    noncon_to_con_tf_path = snakemake.params.tf
    subject = snakemake.params.freesurf_patient
    subjects_dir = snakemake.params.freesurf_dir
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
      tfm_list=[]
      if snakemake.config['noiseDetect']:
        # Manage transforms from mri ras to Freesurfer RAS
        # Read T1 transforms
        t1 = nb.load(t1_path)
        Torig = t1.header.get_vox2ras_tkr()
        tfm_list = [(t1.affine, True), (Torig, False)]
        # Check if there's a tf from non-contrast to contrast
        if noncon_to_con_tf_path != None:
          tfm_list = [(noncon_to_con_tf_path, False)] + tfm_list
      # Call class
      seegTF = cleanSEEG(edf_path, 
                    chn_tsv_path,
                    RmTrendMethod = snakemake.config['detrend_method'],
                    methodPLI = snakemake.config['methodPLI'], 
                    lineFreq = snakemake.config['lineFreq'],
                    bandwidth = snakemake.config['bandwidth'],
                    n_harmonics = snakemake.config['n_harmonics'],
                    noiseDetect = snakemake.config['noiseDetect'],
                    highpass = snakemake.config['highpass'], 
                    maxFlatlineDuration = snakemake.config['maxFlatlineDuration'], 
                    tfm = tfm_list, # Only has to have the tfm from the space in the tsv to the one in the 
                    epoch_autoreject = snakemake.config['epoch_length'],
                    processes = processes)

      # Apply filters
      if snakemake.config['noiseDetect']:
        clean, df_epochs = seegTF.clean_epochs(subject = subject, 
                                                subjects_dir = subjects_dir,
                                                return_interpolated=False, 
                                                write_edf_clean = True,
                                                out_edf_path_clean = out_edf,
                                                write_tsv = True,
                                                out_tsv_path = out_tsv,
                                                verbose = False
                                              )
      else:
        clean = seegTF.clean_epochs(subject = subject, 
                                    subjects_dir = subjects_dir,
                                    return_interpolated=False, 
                                    write_edf_clean = True,
                                    out_edf_path_clean = out_edf,
                                    write_tsv = True,
                                    out_tsv_path = out_tsv,
                                    verbose = False
                                  )
    except:
        logging.exception('Got exception on main handler')
        raise

if __name__=="__main__":
    main()