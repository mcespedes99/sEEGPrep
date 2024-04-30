import mne
import nibabel as nb
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyedflib


def readRegMatrix(trsfPath):
    with open(trsfPath) as (f):
        return np.loadtxt(f.readlines())
    
def get_chn_labels(chn_csv_path, electrodes_edf):
    """Gets label for each electrode.
    Parameters
    ----------
    ch_csv_path : str
        Path to csv containing electrodes positions.

    Returns : list with labels.
    -------
        
    """
    elec_info = pd.read_csv(chn_csv_path, sep='\t')
    labels_csv = elec_info['label'].values.tolist()
    # Filter to get only labels that are also in edf file
    labels = [label for label in labels_csv if label in electrodes_edf]
    return labels

def get_chn_positions(chn_csv_path, electrodes_edf, tfm_list=None):
    """Creates dictionary with the position of each electrode.
    Parameters
    ----------
    ch_csv_path : str
        Path to csv containing electrodes positions.

    Returns : dictionary with position (x,y,z) of each electrode.
    -------
        
    """
    elec_pos = pd.read_csv(chn_csv_path, sep='\t')
    chn_pos = {}
    for i in np.arange(len(elec_pos)):
        label = elec_pos.loc[[i], ['label']].values[0][0]
        if label in electrodes_edf:
            pos = elec_pos.loc[[i], ['x','y','z']].values[0]/1000
            for tfm_path, inv_bool in tfm_list:
                tfm = readRegMatrix(tfm_path)
                if inv_bool:
                    tfm = np.linalg.inv(tfm)
                pos = mne.transforms.apply_trans(tfm, pos)
            pos = pos.tolist()
            chn_pos[label] = pos
    return chn_pos


# Open edf file
edf_in = pyedflib.EdfReader('/home/mcesped/projects/ctb-akhanf/cfmm-bids/Khan/epi_iEEG/ieeg/bids/sub-076/ses-001/ieeg/sub-076_ses-001_task-full_run-01_ieeg.edf')
# Extract labels
labels = edf_in.getSignalLabels()
edf_in.close()

csv_path = '/home/mcesped/projects/ctb-akhanf/cfmm-bids/Khan/clinical_imaging/epi_iEEG/derivatives/seega_coordinates/sub-P076/sub-P076_space-native_SEEGA.tsv'

elec_edf = labels
chn_labels = get_chn_labels(csv_path, elec_edf)

tfm = '/home/mcesped/projects/ctb-akhanf/cfmm-bids/Khan/clinical_imaging/epi_iEEG/derivatives/atlasreg/sub-P076/sub-P076_acq-noncontrast_desc-rigid_from-noncontrast_to-contrast_type-ras_xfm.txt'

chn_pos = get_chn_positions(csv_path, elec_edf, [(tfm, False)])

subject = 'sub-P076'
subjects_dir = '/home/mcesped/projects/ctb-akhanf/cfmm-bids/Khan/clinical_imaging/epi_iEEG/derivatives/fastsurfer/'

subj_trans = mne.coreg.estimate_head_mri_t(subject, subjects_dir)
mri_to_head_trans = mne.transforms.invert_transform(subj_trans)
print('Start transforming mri to head')
print(mri_to_head_trans)

montage_mri = mne.channels.make_dig_montage(chn_pos, coord_frame='mri')
montage = montage_mri.copy()
montage.add_estimated_fiducials(subject, subjects_dir)
montage.apply_trans(mri_to_head_trans)

info = mne.create_info(ch_names=chn_labels,
                        ch_types=['seeg'] * len(chn_labels),
                        sfreq=200)

# compute the transform to head for plotting
trans = mne.channels.compute_native_head_t(montage)
# note that this is the same as:
# ``mne.transforms.invert_transform(
#      mne.transforms.combine_transforms(head_mri_t, mri_mni_t))``

view_kwargs = dict(azimuth=105, elevation=100, focalpoint=(0, 0, -15))
brain = mne.viz.Brain(
    subject,
    subjects_dir=subjects_dir,
    cortex="low_contrast",
    alpha=0.25,
    background="white"
)
brain.add_sensors(info, trans=trans)
brain.add_head(alpha=0.25, color="tan")
brain.show_view(distance=400, **view_kwargs)