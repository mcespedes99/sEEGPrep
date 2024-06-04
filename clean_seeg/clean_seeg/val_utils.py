import pandas as pd 


def get_chn_labels(chn_csv_path, electrodes_edf):
    """Gets label for each electrode.
    Parameters
    ----------
    ch_csv_path : str
        Path to csv containing electrodes positions.

    Returns : list with labels.
    -------

    """
    elec_info = pd.read_csv(chn_csv_path, sep="\t")
    labels_csv = elec_info["name"].values.tolist()
    # Filter to get only labels that are also in edf file
    labels = [label for label in electrodes_edf if label in labels_csv]
    discarded_labels = [label for label in electrodes_edf if label not in labels_csv]
    return labels, discarded_labels