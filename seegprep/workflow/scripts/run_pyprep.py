import pyedflib
import mne
import numpy as np
from pyprep.find_noisy_channels import NoisyChannels
import json
import logging


def get_bad_channels(edf, channels_tsv, max_dur_epoch=20):  # max dur in minutes
    import os
    import pandas as pd

    # Get size of file in GB and number of samples
    size_file = os.path.getsize(edf) / 10**9
    # Get info from channels.tsv (only SEEG, ECOG or EEG)
    tsv_obj = pd.read_csv(channels_tsv, sep="\t")
    channels_tsv = tsv_obj['name'].to_list()
    # Get the first raw data
    edf_obj = pyedflib.EdfReader(edf)
    n_samples = edf_obj.getNSamples().astype(int)[0]
    # Get indexes of relevant channels
    chns_edf = edf_obj.getSignalLabels()
    ch_idx = []
    chn_list = []
    for chn in channels_tsv:
        if chn in chns_edf:
            ch_idx.append(chns_edf.index(chn))
            chn_list.append(chn)
    srate = edf_obj.getSampleFrequencies()[0]

    # Get number of chunks based on available memory or max_dur
    criteria = np.argmax(
        [
            int(np.ceil(n_samples / (max_dur_epoch * 60 * srate))),
            1,
        ]
    )
    if criteria == 0:
        n_extract = int((max_dur_epoch * 60 * srate))
        n_chunks = int(np.ceil(n_samples / n_extract))
    else:
        n_chunks = 1
        n_extract = n_samples

    results = []
    # Analize per chunk
    for n in range(0, n_chunks):
        n_get = min(n_extract, n_samples - n * n_extract)
        signals_process = np.zeros((len(ch_idx), n_get))
        for idx, chn in enumerate(ch_idx):
            signals_process[idx, :] = edf_obj.readSignal(chn, n * n_extract, n_get)
        # Build new raw obj
        # Read EDF
        raw = mne.io.read_raw_edf(edf)
        # Rebuild info
        new_type = [
            type
            for idx, type in enumerate(raw.info.get_channel_types())
            if raw.info.ch_names[idx] in chn_list
        ]
        print(len(new_type), len(chn_list))
        new_info = mne.create_info(
            chn_list, ch_types=new_type, sfreq=raw.info["sfreq"]
        )
        new_info.set_montage(raw.info.get_montage())
        new_raw = mne.io.RawArray(signals_process, new_info)
        del signals_process
        # Run pyprep
        try:
            nd = NoisyChannels(new_raw, do_detrend=True, random_state=1337)
            del new_raw
            nd.find_all_bads(ransac=False)
            key = f"segment_{n+1}"
            results.append(
                {
                    "bad_by_nan": nd.bad_by_nan,
                    "bad_by_flat": nd.bad_by_flat,
                    "bad_by_deviation": nd.bad_by_deviation,
                    "bad_by_hf_noise": nd.bad_by_hf_noise,
                    "bad_by_correlation": nd.bad_by_correlation,
                    "bad_by_SNR": nd.bad_by_SNR,
                    "bad_by_dropout": nd.bad_by_dropout,
                }
            )
            del nd
        except ValueError as ex:
            if "No appropriate channels found for the given picks" in str(ex):
                del new_raw
                results.append(
                    {
                        "bad_by_nan": [],
                        "bad_by_flat": chn_list,
                        "bad_by_deviation": [],
                        "bad_by_hf_noise": [],
                        "bad_by_correlation": [],
                        "bad_by_SNR": [],
                        "bad_by_dropout": [],
                    }
                )
            else:
                raise ex
    # Close file
    edf_obj.close()
    return results


def main():
    edf = snakemake.input.edf
    channels_tsv = snakemake.input.channels_tsv
    # modality = snakemake.params.modality
    out_json = snakemake.output.out_json
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        results = get_bad_channels(edf, channels_tsv)
        # Save
        with open(out_json, "w") as outfile:
            json.dump(results, outfile)
    except:
        logging.exception("Got exception on main handler")
        raise


if __name__ == "__main__":
    main()
