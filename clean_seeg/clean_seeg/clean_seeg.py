from .utils import (
    plot_filter,
    downsampling_mne,
    remove_trend,
)
from .clean_PLI import removePLI_chns, zapline, cleanline, notch_filt
from .data_manager import (
    create_EDF,
    create_bipolars,
    extract_location,
    apply_bipolar_criteria,
    get_chn_info,
    create_epoch_EDF,
)
from .edf_utils import (
    adjust_n_samples,
    find_timestamps,
    extract_signal
)
from .val_utils import (
    get_chn_labels
)
import logging
import pandas as pd
from typing import Literal, List, Tuple, Union, Dict
import pyedflib
import shutil
import os
import re
import snakebids
import bids
import numpy as np
from multiprocessing.pool import Pool
from multiprocessing import get_context
from functools import partial
import json

class cleanSEEG:
    """
    Main class to clean SEEG (Stereo-ElectroEncephaloGraphy) data.

    Attributes:
        edf_path (str): Path to the .edf file containing SEEG data.
        RmTrendMethod (str): Method for removing trends. Options are "HighPass", "LinearDetrend", "Demean".
        methodPLI (str): Method for power line interference removal. Options are "Cleanline", "NotchFilter", "PLIremoval", "Zapline".
        lineFreq (int): Frequency of the power line interference.
        bandwidth (int): Bandwidth for filtering.
        n_harmonics (int): Number of harmonics to consider for removal.
        highpass (List[float]): Highpass filter transition band.
        tfm (List[Tuple[str, bool]]): List of transformations to apply to the contacts coordinates, each defined by a tuple with a string
        identifier and a boolean flag to define whether or not to invert (if set to True, the transform is inverted).
        processes (Optional[int]): Number of processes to use for parallel computation.
    """

    def __init__(
        self,
        edf_path: Union[str, os.PathLike],
        RmTrendMethod: Literal["HighPass", "LinearDetrend", "Demean"] = "HighPass",
        methodPLI: Literal["Cleanline", "NotchFilter", "PLIremoval", "Zapline"] = "Zapline",
        lineFreq: int = 60,
        bandwidth: int = 4,
        n_harmonics:  int = 3,
        highpass: List[float] =[0.5, 1],
        tfm: List[Tuple[str, bool]] = [],
        processes: int = None,
    ) -> None:
        """
        Initializes the cleanSEEG instance with the specified parameters.

        Args:
            edf_path (Union[str, os.PathLike]): Path to the .edf file containing SEEG data.
            RmTrendMethod (Literal["HighPass", "LinearDetrend", "Demean"]): Method for removing trends. Defaults to "HighPass".
            methodPLI (Literal["Cleanline", "NotchFilter", "PLIremoval", "Zapline"]): Method for power line interference removal. Defaults to "Zapline".
            lineFreq (int): Frequency of the power line interference. Defaults to 60.
            bandwidth (int): Bandwidth for filtering. Defaults to 4.
            n_harmonics (int): Number of harmonics to consider for removal. Defaults to 3.
            highpass (List[float]): Highpass filter transition band frequencies. Defaults to [0.5, 1]. If set to None, no filtering is applied.
            tfm (List[Tuple[str, bool]]): List of transformations to apply to each contact coordinate to move it to the parcellation space,
            each defined by a tuple with a string identifier and a boolean flag  to define whether or not to invert (if set to True, the transform is inverted).
            Defaults to an empty list.
            processes (Optional[int]): Number of processes to use for parallel computation. Defaults to None.
        """

        self.edf_path = edf_path
        self.RmTrendMethod = RmTrendMethod
        self.methodPLI = methodPLI
        self.lineFreq = lineFreq
        self.bandwidth = bandwidth
        self.n_harmonics = n_harmonics
        self.highpass = highpass  # Set to None to shut down
        self.tfm = tfm  
        self.processes = processes
        # Extra params to allow computation of steps in series
        self.reref_chn_list = []
        self.reref_chn_df = []
        self.inter_edf = None
        self.clean_edf = None
        self.subject = (None,)
        self.subjects_dir = None

        # Find sample rate: Requires pyedflib <= 0.1.30
        with pyedflib.EdfReader(edf_path) as edf_in:
            self.srate = edf_in.getSampleFrequencies()[0] / edf_in.datarecord_duration
    
    
    # Epoch extraction
    def extract_epochs(
        self,
        n_samples: int = None,
        event_dur: float =None,
        event_start: Union[str, int] = None,
        out_root: Union[str, os.PathLike] = None,
        out_files: Union[List[str], List[os.PathLike]] = None,
        tmpdir: Union[str, os.PathLike] = None,
        snakemake: bool = False,
        return_report: bool = True
    ) -> Union[pd.DataFrame, None] :
        """
        Extract epochs from EDF files based on specified start positions and durations.

        This function extracts clips from EDF files given an initial position (specified by annotations or an index) 
        and a duration (either in seconds or samples). It produces EDF files with the corresponding epochs and can 
        return a small report summarizing the events found.

        Parameters
        ----------
        n_samples : int, optional
            Number of samples to extract. If not specified, `event_dur` must be provided.
        event_dur : float, optional
            Duration of the event in seconds. Used to calculate `n_samples` if not provided.
        event_start : int or str, optional
            Starting point of the event. Can be an index or a label to look for an annotation in the EDF file.
            If not given, the EDF is splitted into multiple epochs of the specified samples or duration.
        out_root : str, optional
            Root directory to save the output EDF files. It is used to create the output files names
            if out_files is set to None.
        out_files : list of str or os.PathLike, optional
            List of output file paths. If not specified, files are named automatically using out_root.
        tmpdir : str or os.PathLike, optional
            Temporary directory to use for intermediate file storage.
        snakemake : bool, optional
            If True, modify the suffix of the output files to include '_tmp'.
        return_report : bool, optional
            If True, return a report summarizing the extracted epochs.

        Returns
        -------
        pandas.DataFrame or None
            DataFrame containing the report with columns for 'Epoch start event', 'Relative start time', 
            'Duration', and 'Clip number'. If `return_report` is False, returns None.
        
        Raises
        ------
        ValueError
            If neither `n_samples` nor `event_dur` is provided.
        
        Notes
        -----
        This function assumes the EDF file is accessible and readable. Make sure to handle file paths and permissions appropriately.
        
        Examples
        --------

        """
        if not n_samples and not event_dur:
            raise ValueError("Either n_samples or event_dur must be provided")
         
        # Initialize report
        report = {
            'Epoch start event': [],
            'Relative start time': [],
            'Duration': [],
            'Clip number': []
        }
        
        with pyedflib.EdfReader(self.edf_path) as f:
            # Number of samples in edf. Use the samples from the first channel as reference.
            N = f.getNSamples()[0]
            # Sampling rate. Use the samples from the first channel as reference.
            srate = f.getSampleFrequencies()[0] / f.datarecord_duration


            if not n_samples:
                # Get number of samples based on event duration
                n_samples = int(event_dur * srate)
            # Check than n_samples is not bigger than N
            if n_samples > N:
                logging.warning('Number of samples bigger than number of points in file, changing to the max.')
                n_samples = N
            
            
            # Adjust n_samples depending on the datarecord duration and the number of samples in file
            # Use the samples from the first channel as reference.
            n_samples = adjust_n_samples(n_samples, f)
            time_stamps_init = find_timestamps(event_start, f, srate, n_samples)

        # Copy file to local scratch if possible
        new_edf = None
        if tmpdir != None and os.path.exists(tmpdir):
            # Copy file to local scratch
            file_name = os.path.basename(self.edf_path)
            new_edf = os.path.join(tmpdir, file_name)
            shutil.copy(self.edf_path, new_edf)

        # Extract entities from input path
        entities = bids.layout.parse_file_entities(self.edf_path)
        # Combine 'extension' with 'suffix' and delete the first one
        if snakemake:
            entities["suffix"] = entities["suffix"] +'_tmp'+ entities["extension"]
        else:
            entities["suffix"] = entities["suffix"] + entities["extension"]
        del entities["extension"]
        # Add 'clip' entity
        entities["rec"] = "clip"

        # Create new EDF files for each event
        report = self._create_epoch_edf_files(time_stamps_init,
                                             n_samples,
                                             out_files,
                                             out_root,
                                             new_edf,
                                             entities,
                                             report,
                                             event_start)

        
        if new_edf is not None:
            os.remove(new_edf)
        
        if return_report:
            return pd.DataFrame(report)

        return None
   
   
    def downsample(
        self,
        channels_tsv: Union[str, os.PathLike],
        target_srate: float,
        out_edf_path: Union[str, os.PathLike, None] = None,
        return_report: bool = True
    ) -> Union[Tuple[np.ndarray, pd.DataFrame], np.ndarray]:
        """
        Downsamples the signal data in an EDF file.

        This function downsamples the signal data in an EDF file to a specified target sampling rate. It can 
        optionally write the downsampled data to a new EDF file and return a report summarizing the downsampling 
        process.

        Parameters
        ----------
        channels_tsv : str or os.PathLike
            Path to a TSV file containing channel information.
        target_srate : float
            Target sampling rate for downsampling.
        out_edf_path : str or os.PathLike, optional
            If not None, this path is used to save the downsampled EDF file.
        return_report : bool, optional
            If True, returns a report summarizing the downsampling process. Default is True.

        Returns
        -------
        Tuple[np.ndarray, pd.DataFrame] or np.ndarray
            If `return_report` is True, returns a tuple with the downsampled data and a DataFrame summarizing 
            the downsampling process. Otherwise, returns only the downsampled data.

        Raises
        ------
        ValueError
            If `channels_tsv` is not provided.

        Examples
        --------
        >>> TODO
        """
        if channels_tsv is None:
            raise ValueError("Please indicate the path to a TSV file with channels information.")
        
        # Open EDF file
        with pyedflib.EdfReader(self.edf_path) as reader:
            # Extract labels to only process signals in tsv file
            labels = reader.getSignalLabels()
            chn_labels, discarded_labels = get_chn_labels(channels_tsv, labels)
            # Number of samples based on the first channel
            N = reader.getNSamples()[0]

        # create a process context. Refer to:
        # https://github.com/dask/dask/issues/3759
        ctx = get_context("spawn")
        with Pool(processes=self.processes, context=ctx) as pool:
            data_list, newSrate = zip(
                *pool.map(
                    partial(
                        downsampling_mne,
                        edf_file=self.edf_path,
                        target_srate=target_srate,
                    ),
                    chn_labels,
                )
            )
        
        # Reformatting to array
        data_dnsampled =  np.array(data_list)
        # n_samples = len(data_list[0])
        # data_dnsampled = np.zeros((len(chn_labels), n_samples))
        # for ch in np.arange(len(data_list)):
        #     data_dnsampled[ch, :] = data_list[ch]

        # Write edf if requested
        if out_edf_path is not None:
            if os.path.exists(out_edf_path):
                logging.warning(f"EDF file {out_edf_path} will be overwritten.")
            create_EDF(
                self.edf_path,
                out_edf_path,
                self.processes,
                signal=data_list,
                new_srate=newSrate[0],
                chn_labels=chn_labels
            )
            del data_list
            self.downsample_edf = out_edf_path
        
        report = {
            'Original sampling rate': self.srate,
            'Target sampling rate': target_srate,
            'New sampling rate': newSrate[0],
            'Discarded channels':discarded_labels
        }
        
        if return_report:
            return data_dnsampled, pd.DataFrame([report])
        
        return data_dnsampled

    def drift_correction(
        self,
        channels_tsv: Union[str, os.PathLike],
        out_edf_path_clean: Union[str, os.PathLike] = None,
    )-> Tuple[np.ndarray, pd.DataFrame, Dict[str, Union[str, List[str], List[float]]]]:
        """
        Applies drift correction to a given EDF file.

        This function applies drift correction to the signal data in an EDF file and produces a new 
        file with the clean signal. A report with relevant information for quality control (QC) of
        the results is also returned.

        Parameters
        ----------
        channels_tsv : str or os.PathLike
            Path to a TSV file containing channel information 
            (*_channels.tsv file according to the BIDS standard).
        out_edf_path_clean : str or os.PathLike, optional
            Path to save the cleaned EDF file. If not provided, the cleaned file is not saved.

        Returns
        -------
        clean : np.ndarray
            The cleaned signal data.
        df_report : pd.DataFrame
            DataFrame containing the original and new mean values of the signals.
        report : dict
            Dictionary containing the method used for detrending and any discarded channels.

        Raises
        ------
        ValueError
            If `channels_tsv` is not provided.
        """
        # Validate channels_tsv
        if channels_tsv is None:
            raise ValueError("Please indicate the path to a TSV file with channels information.")
        
        
        chn_labels, discarded_labels, signal = extract_signal(channels_tsv, self.edf_path)

        clean, filt_params = remove_trend(signal, self.RmTrendMethod, self.srate, self.highpass)

        # Write edfs if requested
        if out_edf_path_clean is not None:
            # Clean signal
            # convert to list
            clean_list = [clean[i, :] for i in range(clean.shape[0])]
            if os.path.exists(out_edf_path_clean):
                logging.warning(
                    f"EDF file {out_edf_path_clean} will be overwritten."
                )
            create_EDF(
                self.edf_path,
                out_edf_path_clean,
                self.processes,
                chn_labels=chn_labels,
                signal=clean_list,
            )
            del clean_list
            self.clean_edf = out_edf_path_clean

        # Create dataframe with relevant info to QC
        df_report = pd.DataFrame({
            'Channel': chn_labels,
            'Original mean': np.mean(signal, axis=1),
            'New mean': np.mean(clean, axis=1)
        })
        
        # Now create report with the characteristics of the method applied
        report = {
            'Method used': self.RmTrendMethod,
            'Discarded channels': discarded_labels
        }
        if self.highpass is not None:
            report['Transition band (Hz)'] = self.highpass
        
        # Finally plot the filter if exists
        if len(filt_params)>0:
            plot_filter(filt_params[0], filt_params[1], os.path.dirname(out_edf_path_clean), self.highpass)
        
        return clean, df_report, report

    # Reference function
    def rereference(
        self,
        electrodes_tsv: Union[str, os.PathLike],
        out_edf_path: Union[str, os.PathLike],
        out_tsv_path: Union[str, os.PathLike] = None,
        return_report: bool = True,
    ) -> Union[Tuple[pd.DataFrame, Dict[str, List[str]]], None]:

        """
        Reference a given monopolar EDF to a bipolar scheme.

        This function converts a monopolar EDF file to a bipolar scheme based on electrode names
        provided in a TSV file. It can also produce a TSV file with the bipolar channel information 
        and return a report summarizing the discarded channels.

        Parameters
        ----------
        electrodes_tsv : str or os.PathLike
            Path to a TSV file containing channel information
            (*electrodes.tsv file according to the BIDS standard).
        out_edf_path : str or os.PathLike
            Path to save the output EDF file with the bipolar scheme.
        out_tsv_path : str or os.PathLike, optional
            Path to save the TSV file with the bipolar channel information. Required if it desired to output the table.
        return_report : bool, optional
            If True, returns a report summarizing the bipolar combinations and discarded channels. Default is True.

        Returns
        -------
        Tuple[pd.DataFrame, Dict[str, List[str]]] or None
            If `return_report` is True, returns a tuple with a DataFrame of bipolar combinations and a dictionary
            with discarded channels. Otherwise, returns None.

        Raises
        ------
        ValueError
            If `electrodes_tsv` is not provided.

        Examples
        --------
        >>> TODO
        """
        # Manage some errors:
        if electrodes_tsv is None:
            raise ValueError("Please indicate the path to a TSV file with channels information.")

        # Extract info about electrodes positions
        elec_pos = pd.read_csv(electrodes_tsv, sep="\t")
        with pyedflib.EdfReader(self.edf_path) as reader:
            # Extract channels present in edf file (might not be all)
            elec_edf = reader.getSignalLabels()
        
        # First check which electrodes from the TSV are in the EDF
        _, discarded_labels = get_chn_labels(electrodes_tsv, elec_edf)
        if len(discarded_labels) >= len(elec_pos):
            raise ValueError("The EDF file is in the wrong format or it's already bipolar.")
        
        # Create bipolar combinations
        bipolar_channels, bipolar_info_df = create_bipolars(
            elec_pos, elec_edf, self.processes
        )

        # Save values in the class
        self.reref_chn_list = bipolar_channels
        self.reref_chn_df = bipolar_info_df

        # Write TSV
        if out_tsv_path is not None: 
            if os.path.exists(out_tsv_path):
               logging.warning(f"tsv file {out_tsv_path} will be overwritten.")
            bipolar_info_df.to_csv(out_tsv_path, index=False, sep="\t")

        # print(bipolar_channels)

        # Create new EDF file
        if out_edf_path is not None:
            if os.path.exists(out_edf_path):
                logging.warning(f"edf file {out_edf_path} will be overwritten.")
            create_EDF(
                self.edf_path, out_edf_path, self.processes, chn_labels=bipolar_channels
            )
            self.rereference_edf = out_edf_path
        if return_report:
            bipolar_combs = {
                'Bipolar channel': [values[0] for values in bipolar_channels],
                'Unipolar channels': [values[1:] for values in bipolar_channels]
            }
            report = {
                'Discarded channels': discarded_labels
            }
            return pd.DataFrame(bipolar_combs), report
        
        return None

    # Function to reject PLI
    def reject_PLI(
        self,
        channels_tsv: Union[str, os.PathLike],
        out_edf_path: Union[str, os.PathLike, None] = None,
        return_report: bool = True
    ) -> Union[Tuple[np.ndarray, Dict[str, Union[List[str], str, int, float]]], np.ndarray]:
        """
        Applies power line interference (PLI) attenuation using one of four methods.

        This function attenuates power line interference (PLI) in the signal data of an EDF file using one of four methods:
        Cleanline, Zapline, rejectPLI, or Notch Filtering. It can optionally write the cleaned data to a new EDF file and
        return a report summarizing the PLI removal process.

        Parameters
        ----------
        channels_tsv : str or os.PathLike
            Path to a TSV file containing channel information.
        write_edf : bool, optional
            If True, writes the cleaned data to a new EDF file. Default is False.
        out_edf_path : str or os.PathLike, optional
            Path to save the cleaned EDF file. Required if `write_edf` is True.
        return_report : bool, optional
            If True, returns a report summarizing the PLI removal process. Default is True.

        Returns
        -------
        Tuple[np.ndarray, Dict[str, Union[List[str], str, int, float]]] or np.ndarray
            If `return_report` is True, returns a tuple with the cleaned data and a dictionary summarizing the PLI removal
            process. Otherwise, returns only the cleaned data.

        Raises
        ------
        ValueError
            If `channels_tsv` is not provided.

        Examples
        --------
        >>>  TODO
        """

        # Manage a few exceptions:
        if channels_tsv is None:
            raise ValueError("Please indicate the path to a TSV file with channels information.")

        # Get signal and labels
        chn_labels, discarded_labels, signal = extract_signal(channels_tsv, self.edf_path)
        
        # Remove line noise
        logging.info("Removing line noise")
        clean = self._wrapper_PLI(signal)
        logging.info("PLI removal completed.")

        # Write edfs if requested
        if out_edf_path is not None:
            # Clean signal
            # convert to list
            clean_list = [clean[i, :] for i in range(clean.shape[0])]
            if os.path.exists(out_edf_path):
                logging.warning(f"EDF file {out_edf_path} will be overwritten.")
            create_EDF(
                self.edf_path,
                out_edf_path,
                self.processes,
                chn_labels=chn_labels,
                signal=clean_list,
            )
            del clean_list
            self.clean_edf = out_edf_path
        
        # Create report if requested
        if return_report:
            report = {
                'Discarded channels': discarded_labels,
                'Method': self.methodPLI,
                'Line frequency (Hz)': self.lineFreq
            }
            if self.methodPLI == "Cleanline":
                report['Bandwidth (Hz)'] = self.bandwidth
            elif self.methodPLI == "NotchFilter":  # TODO: add bandwidth param
                report['Number of Harmonics'] = self.n_harmonics
            elif self.methodPLI == "PLIremoval":
                report.update({
                    'Number of Harmonics': 3,
                    'B0, Initial notch bandwidth of the frequency estimator': 100,
                    'Binf, Asymptotic notch bandwidth of the frequency estimator': 0.01,
                    'Bst, Rate of convergence to 95# of the asymptotic bandwidth Binf': 4,
                    'P0, Initial settling time of the frequency estimator': 0.1,
                    'Asymptotic settling time of the frequency estimator': 2,
                    'Pst, Rate of convergence to 95# of the asymptotic settling time': 5,
                    'W, Settling time of the amplitude and phase estimator': 2
                })
                # Hardcoded for now
            return clean, report
        
        return clean
    
    def identify_regions(
        self,
        electrodes_tsv: Union[str, os.PathLike],
        segmentation_path: Union[str, os.PathLike],
        colortable_file: Union[str, os.PathLike],
        reference: Literal['bipolar', 'unipolar'],
        use_reref: bool = True,
        out_tsv_path: Union[str, os.PathLike, None] = None,
        discard_wm_un: bool = False,
        out_edf_path: Union[str, os.PathLike, None] = None,
        vol_version: bool = False,
        json_out: Union[str, os.PathLike, None] = None,
        mask_out: Union[str, os.PathLike, None] = None,
        return_report: bool = True
    ) -> Union[Tuple[pd.DataFrame, pd.DataFrame, Dict[str, List[str]]], pd.DataFrame]:
        """
        Estimate brain regions for each channel based on a given parcellation file.

        This function estimates brain regions for each channel using a given parcellation file. It supports both
        monopolar and bipolar channels and two methods: volumetric approach and point-based approach.

        Parameters
        ----------
        electrodes_tsv : str or os.PathLike
            Path to a TSV file containing channel information
            (*electrodes.tsv file according to the BIDS standard).
        segmentation_path : str or os.PathLike
            Path to the file containing parcellation information that can be loaded with Nibabel 
            and contains a header and affine transform.
        colortable_file : str or os.PathLike
            Path to a TSV file containing colortable information to assign labels and colors to the
            values in given segmentation file.
        reference : Literal['bipolar', 'unipolar']
            Reference of the channels in the EDF file. Must be either 'bipolar' or 'unipolar'.
        use_reref : bool, optional
            If True, previous information calculated in the re-referencing step is used. 
            Default is True. If True, re-reference must have been ran before.
        out_tsv_path : str, os.PathLike or None, optional
            If not None, indicates the path to save the TSV file with channel location information.
        discard_wm_un : bool, optional
            If True, discards channels from the EDF file estimated to be in white matter and unknown
            regions based on the given parcellation. Default is False.
        out_edf_path : str, os.PathLike or None, optional
            If not None, indicates the path to save the cleaned EDF file with the non-discarded channels.
        vol_version : bool, optional
            If True, uses volumetric approach for region estimation, else, the point-based method
            is used. Default is False.
        json_out : str, os.PathLike or None, optional
           If not None, indicates the path to save the JSON file with region information for each
           channel in the given EDF file. Only applicable if `vol_version` is True.
        mask_out : sstr, os.PathLike or None, optional
            If not None, indicates the path to save the mask file created during the volumetric estimation
            of the regions for each channel. Default is None. Only applicable if `vol_version` is True.
        return_report : bool, optional
            If True, returns a report summarizing the region estimation process. Default is True.

        Returns
        -------
        Tuple[pd.DataFrame, pd.DataFrame, dict] or pd.DataFrame
            If `return_report` is True, returns a tuple with DataFrames for channel coordinates and regions per channel, and a dictionary
            summarizing the discarded channels. Otherwise, returns only the DataFrame for channel coordinates.

        Raises
        ------
        ValueError
            If `reference` is not 'bipolar' or 'unipolar'.
            If `electrodes_tsv` is not provided.

        Examples
        --------
        >>> TODO
        """
        # Manage exception
        if reference not in ['bipolar', 'unipolar']:
            raise ValueError("Reference must be either 'bipolar' or 'unipolar'.")

        if electrodes_tsv is None:
            raise ValueError("Please indicate the path to a TSV file with electrodes information.")

        if use_reref:
            if len(self.reref_chn_list) == 0 or len(self.reref_chn_df) == 0:
                raise ValueError("Please run rereference first or set use_reref to False.")
            chn_info_df = self.reref_chn_df
            chn_list = self.reref_chn_list
        # Extract electrodes information if not using rereference results
        else:
            logging.info(f"Extracting channel positions from {electrodes_tsv}")
            # Extract info about electrodes positions
            elec_pos = pd.read_csv(electrodes_tsv, sep="\t")
            # Filter any channels with nan vals as coordinates
            nan_rows = elec_pos[['x', 'y', 'z']].isna().any(axis=1)
            elec_pos = elec_pos.loc[~nan_rows]
            dropped_chns = elec_pos.loc[nan_rows, 'name'].values.tolist()
            
            # First extract list of channels present in edf file
            with pyedflib.EdfReader(self.edf_path) as edf:
                elec_edf = edf.getSignalLabels()

            # Check if needs to preprocess electrodes.tsv file
            if reference=='bipolar':
                # Needs to convert electrodes df to bipolar
                # Create bipolar combinations
                _, elec_pos = create_bipolars(
                    elec_pos, elec_edf, self.processes, compare_edf=False
                )

            chn_info_df, chn_list, discarded_chns = get_chn_info(
                elec_pos, elec_edf, reference
            )  

        # Create tsv file with information about location of the channels
        df_location, regions_per_chn = extract_location(
            segmentation_path,
            chn_info_df,
            self.tfm,
            colortable_file,
            reference,
            vol_version,
            mask_out
        )


        if out_tsv_path is not None:
            if os.path.exists(out_tsv_path):
                logging.warning(f"TSV file {out_tsv_path} will be overwritten.")
            df_location.to_csv(out_tsv_path, index=False, sep="\t")
        
        if json_out and vol_version:
            with open(json_out, "w") as outfile:
                json.dump(regions_per_chn, outfile)
        
        # Discard data from white matter and unknown if requested
        if discard_wm_un:
            chn_list = apply_bipolar_criteria(df_location, chn_list, self.processes)
            # Overwrite reref info if previously run
            if use_reref:
                self.reref_chn_list = chn_list


            # Write edf file if requested
            if out_edf_path is not None:
                if os.path.exists(out_edf_path):
                    logging.warning(f"EDF file {out_edf_path} will be overwritten.")
                create_EDF(
                    self.edf_path, out_edf_path, self.processes, chn_labels=chn_list
                )
                self.rereference_edf = out_edf_path
        elif (not discard_wm_un) and (out_edf_path is not None):
            logging.warning(
                f"EDF file {out_edf_path} will only be copied from the original EDF as no updates have been made."
            )
            shutil.copyfile(self.edf_path, out_edf_path)
        
        if return_report:
            unique_regions = df_location['region name'].unique()
            region_maps = {'Region': unique_regions, 'Number of channels': [], 'Channels':[]}
            for region in unique_regions:
                chns = df_location.loc[df_location['region name']==region, "name"].values
                region_maps['Number of channels'].append(len(chns))
                region_maps['Channels'].append(chns)
            df_regions = pd.DataFrame(region_maps)
            report = {
                'Discarded channels': discarded_chns+dropped_chns
            }
            return df_location, df_regions, report
        
        return df_location
        

    # Wrapper to manage different PLI methods
    def _wrapper_PLI(self, signal: np.ndarray) -> np.ndarray:
        """
        Applies the specified PLI attenuation method to the signal.

        This function applies the specified power line interference (PLI) attenuation method to the input signal. 
        The method is determined by the `methodPLI` attribute of the class.

        Parameters
        ----------
        signal : np.ndarray
            Input signal data to be processed.

        Returns
        -------
        np.ndarray
            Signal data after PLI attenuation.

        Raises
        ------
        ValueError
            If the specified PLI method is not valid.

        Examples
        --------
        >>> clean_signal = _wrapper_PLI(signal)
        """
        if self.methodPLI == "Cleanline":
            signal = cleanline(
                signal.T, self.srate, processes=self.processes, bandwidth=self.bandwidth
            )
        elif self.methodPLI == "Zapline":
            signal = zapline(signal.T, self.lineFreq / self.srate, self.srate)
        elif self.methodPLI == "NotchFilter":  # TODO: add bandwidth param
            signal = notch_filt(
                signal.T, self.lineFreq, self.srate, n_harmonics=self.n_harmonics
            )
        elif self.methodPLI == "PLIremoval":
            signal = removePLI_chns(
                signal.T,
                self.srate,
                3,
                [100, 0.01, 4],
                [0.1, 2, 5],
                2,
                processes=self.processes,
                f_ac=self.lineFreq,
            )  # Hardcoded for now, TODO
        else:
            raise ValueError("PLI method not valid.")
        signal = signal.T
        return signal

    def _create_epoch_edf_files(
        self,
        time_stamps_init: List[float],
        n_samples: int,
        out_files: Union[List[str], os.PathLike],
        out_root: str,
        new_edf: Union[str, None],
        entities: Dict[str, str],
        report: Dict[str, List],
        event_start: Union[str, int]
    ) -> Dict[str, List]:
        """
        Create new EDF files with the found timestamps given the events of interest.

        This function creates new EDF files with the found timestamps given the events of interest in the original EDF file.
        The function also updates the report with details about each epoch.

        Parameters
        ----------
        time_stamps_init : list of float
            List of initial timestamps for the events of interest.
        n_samples : int
            Number of samples to extract for each epoch.
        out_files : list of str or os.PathLike
            List of output file paths. If not specified, files are named automatically.
        out_root : str
            Root directory to save the output EDF files.
        new_edf : str or None
            Path to a temporary EDF file, if any.
        entities : dict of str
            Dictionary of entities for naming the output files.
        report : dict of list
            Dictionary to store the report details about each epoch.
        event_start : str or int
            Starting point of the event. Can be an index or a label.

        Returns
        -------
        dict of list
            Updated report with details about each epoch.
        """
        for index, event_timestamp in enumerate(time_stamps_init):
            if out_files is None:
                # New file name
                entities["clip"] = f"{index+1:02}"
                out_edf_path = os.path.basename(
                    snakebids.bids(root=out_root, **entities)
                )
                out_edf_path = os.path.join(out_root, out_edf_path)
            else:
                # Filter based on regex
                reg = re.compile(f"clip-{index+1:02}")
                out_edf_path = list(filter(reg.search, out_files))[0]
            
            if new_edf is None:
                time_init, duration = create_epoch_EDF(
                                        self.edf_path,
                                        event_timestamp,
                                        n_samples,
                                        out_edf_path,
                                        self.processes,
                )
            else:
                time_init, duration = create_epoch_EDF(
                    new_edf, event_timestamp, n_samples, out_edf_path, self.processes
                )
            
            # Save report
            report["Epoch start event"].append(event_start)
            report['Relative start time'].append(time_init)
            report["Duration"].append(duration)
            report['Clip number'].append(f"clip-{index+1:02}")
        return report