import logging
import shutil
import os
import datetime
from ansi2html import Ansi2HTMLConverter
import json
import numpy as np
import pandas as pd
import plotly.express as px
from clean_seeg.autoreject.autoreject import read_reject_log
import plotly.express as px
import plotly.graph_objects as go
import pyedflib


def get_filt_figures(figures):
    # figures: list of tuples with file, title
    html_code = f"""
        <h2 class="section__heading">Highpass filter information</h2>
        <p class="section__description">Suspendisse potenti. Quisque blandit urna vitae maximus tempor.
    """ 
    for fig_path, title in figures:
        template = f"""
                <article class="article">
                    <h2 class="article__heading">{title}</h2>
                    <img src="{fig_path}" alt="Article Image" class="article__image" style="display: block; margin-left: auto; margin-right: auto; max-width: 15cm;">
                </article>
            """
        html_code = html_code + template
    return html_code
   

def append_to_template(template_file, div_id, content_to_append, className):
    # Read existing template
    with open(template_file, 'r') as f:
        template_content = f.read()

    # Find the specific div in the template using regex or string matching
    div_start = f'<div id="{div_id}" class="{className}">'
    div_end = f'</div>'
    div_start_index = template_content.find(div_start)

    if div_start_index != -1:
        # Find the closing tag of the div
        div_end_index = template_content.find(div_end, div_start_index)

        if div_end_index != -1:
            # Insert the content inside the div
            updated_template = (
                template_content[: div_start_index + len(div_start)]
                + content_to_append
                + template_content[div_end_index:]
            )

            # Save the updated template back to the file
            with open(template_file, 'w') as f:
                f.write(updated_template)
        else:
            print(f"Error: Closing tag {div_end} not found in the template.")
    else:
        print(f"Error: Opening tag {div_start} not found in the template.")

def process_raw(edf_path):
    edf = pyedflib.EdfReader(edf_path)
    # Duration
    delta = datetime.timedelta(seconds=edf.file_duration)
    base_time = datetime.datetime(1900, 1, 1)
    duration = (base_time + delta).strftime("%H:%M:%S.%f")
    html_code = f"""
                <table class="info">
                    <tbody>
                        <tr>
                            <th>Filename</th>
                            <td>{edf.file_name}</td>
                        </tr>
                        <tr>
                            <th>Measurement Date</th>
                            <td>{str(edf.getStartdatetime())}</td>
                        </tr>
                        <tr>
                            <th>Sampling frequency (Hz)</th>
                            <td>{edf.getSampleFrequencies()[0] / edf.datarecord_duration}</td>
                        </tr>
                        <tr>
                            <th>Signals in file</th>
                            <td>{edf.signals_in_file}</td>
                        </tr>
                        <tr>
                            <th>Duration</th>
                            <td>{duration} (HH:MM:SS)</td>
                        </tr>
                    </tbody>
                </table>
    """
    edf.close()
    return html_code


def process_val_txt(val_txt):
    # Read file
    with open(val_txt, "r") as f:
        val_str = f.read()
    # Convert to html
    conv = Ansi2HTMLConverter()
    my_html = conv.convert(val_str, full=True)
    # Adjust a few details
    my_html = my_html.replace(
        "body_foreground { color: #AAAAAA; }", "body_foreground { color: #000000; }"
    )
    my_html = my_html.replace(
        "body_background { background-color: #000000; }",
        "body_background { background-color: #FFFFFF; }",
    )
    return my_html


def process_extraval(extraval_txt):
    # Read file
    with open(extraval_txt, "r") as f:
        f.seek(0)
        bids_val_results = f.readlines()
    # Delete \n
    for line in bids_val_results:
        line.replace("\n", "")
    # convert to html
    results_val = """<ul class="list">
    """
    passed = []
    errors = []
    for error in bids_val_results:
        error = error.replace(
            "Warning:", '<span>Warning:</span>'
        )
        error = error.replace("Info:", '<span>Info:</span>')
        if "Info" in error:
            passed.append(error)
        else:
            errors.append(error)
    results_val += f""" <div class="list-section success">"""
    for error in passed:
        results_val += f"""
        <li class="li-success">{error}</li>
        """
    results_val += """</div> <div class="list-section error">"""
    for error in errors:
        results_val += f"""
        <li class="li-error">{error}</li>
        """
    results_val += """</div>"""

    results_val += """</ul>
    """
    return results_val


def process_channel_metrics(pyprep_json, pyprep_html):
    # Read file
    with open(pyprep_json) as f:
        results = json.load(f)
    # For the report
    metrics = list(results[0].keys())
    summary = dict()
    n_total = len(results)
    for metric in metrics:
        bad_total = []
        for dictionary in results:
            bad_total += dictionary[metric]
        bad_total = np.array(bad_total)
        n_per_chn = []
        for chn in np.unique(bad_total):
            n_bad = len(bad_total[bad_total == chn])
            if n_bad > int(np.floor(n_total*0.25)):
                n_per_chn.append((chn, n_bad))
        # Sort it
        n_per_chn = sorted(n_per_chn, key=lambda x: x[1], reverse=True)
        if len(n_per_chn) != 0:
            summary[metric] = dict()
            for chn, n in n_per_chn:
                summary[metric][chn] = np.round(n * 100 / n_total, 3)

    if len(summary) > 0:
        # Convert to figures
        fig_list = []
        for metric in summary:
            tmp_dict = {
                "Channel": list(summary[metric].keys()),
                f"Probability of bad by {metric.split('_')[-1]} (%)": list(
                    summary[metric].values()
                ),
            }
            tmp_df = pd.DataFrame(tmp_dict)
            fig = px.bar(
                tmp_df,
                x="Channel",
                y=f"Probability of bad by {metric.split('_')[-1]} (%)",
            )
            fig_list.append(fig)
        # Save html
        with open(pyprep_html, "w") as f:
            for fig in fig_list:
                f.write(fig.to_html(full_html=False, include_plotlyjs="cdn"))
        # Read the html text
        with open(pyprep_html, "r") as f:
            plotly_html = f.read()
    # Case where no channel is tagged
    else:
        plotly_html = """<ul>
        <li>No channels were tagged as noisy by Pyprep.</li>
        </ul>
        """
    return plotly_html


def plot_heatmap(img, time_vec, ch_names, color_title):
    # Create figure
    fig = go.Figure()

    fig_imshow = px.imshow(
        img,
        labels=dict(x="Time (s)", y="Channels", color=color_title),
        x=time_vec,
        y=ch_names,
    )
    fig.add_trace(go.Heatmap(fig_imshow.data[0]))

    fig.update_xaxes(
        showgrid=True,
        ticklabelposition="outside left",
        tickson="boundaries",
        tickformat="%H:%M:%S.%f",
    )
    fig.update_layout(xaxis_rangeslider_visible=True, yaxis=dict(
            fixedrange=False  # Allow zooming along the y-axis
            ))
    return fig

def plot_PSD_channels(df1, df2, names):
    import plotly.express as px
    from plotly.subplots import make_subplots
    import warnings

    # Ignore warnings
    warnings.filterwarnings('ignore')

    # Create subplots with shared y-axis
    fig_subplots = make_subplots(rows=1, cols=2, shared_yaxes=True,
                                 subplot_titles=[f"{name} signals" for name in names])

    showlegend=True
    for i, df in enumerate([df1, df2]):
        # Generate Plotly Express line plot
        fig = px.line(df, y=df.columns[df.columns!='Frequency (Hz)'], x='Frequency (Hz)')

        # Add traces from original plot to the subplots
        for trace in fig.data:
            # Assign unique legend group names to each trace to avoid duplicate labels
            fig_subplots.add_trace(trace.update(showlegend=showlegend), row=1, col=i+1)
        showlegend=False

    min_f = min(np.min(df1.loc[:,'Frequency (Hz)'].values),np.min(df2.loc[:,'Frequency (Hz)'].values))
    max_f = min(np.max(df1.loc[:,'Frequency (Hz)'].values),np.max(df2.loc[:,'Frequency (Hz)'].values))
    # fig_subplots.update_xaxes(range=[min_f, max_f])
    # Update layout
    fig_subplots.update_layout(
        title="Normalized PSD",
        height=400,
        xaxis_range=[min_f, max_f]
    )

    script_html = fig_subplots.to_html(full_html=False, include_plotlyjs="cdn")
    # Restore warnings to default behavior
    warnings.filterwarnings('default')

    return script_html

def get_PSD(edf_path:str, length_segment: float=3.0):
    import pyedflib
    from scipy.signal import welch
    edf = pyedflib.EdfReader(edf_path)
    labels = edf.getSignalLabels()
    srate = edf.getSampleFrequencies()[0]/edf.datarecord_duration

    # Read signals
    signals = []
    for i in range(len(labels)):
        signals.append(edf.readSignal(i))
    signals = np.array(signals)
    edf.close()

    f, welchpow = welch(signals, fs=srate, nperseg=int(length_segment*srate), axis=1)
    welchpow = np.divide(welchpow, np.sqrt(np.sum(welchpow**2, axis=1)).reshape(welchpow.shape[0],1))
    df = pd.DataFrame(welchpow.T, columns=labels)
    df.insert(0, 'Frequency (Hz)', f)
    return df

def process_json(json_file):
    # Read file
    with open(json_file, "r") as f:
        data = json.load(f)
    # convert to html
    html_code = """<ul>
    """
    for metric in data:
        message = f'<span style="font-weight: bold;">{metric}:</span> {data[metric]}'

        html_code += f"""
        <li>{message}</li>
        """
    html_code += """</ul>
    """
    return html_code

def main():
    # Get all inputs
    final_edf = snakemake.input.edf
    # For val
    val_txt = snakemake.input.val_txt
    extraval_txt = snakemake.input.extraval_txt
    
    # Extra
    # pyprep_json = snakemake.input.pyprep_json
    # pyprep_html = snakemake.params.pyprep_html
    out_html = snakemake.output.out_html
    out_dir = os.path.dirname(out_html)
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        # Create copy of files
        with open(os.path.join(current_dir, "../../resources/html_template/template.html"), "r") as f:
            template_content = f.read()
        with open(out_html, "w") as f:
            f.write(template_content)
        
        for filename in ['app.css', 'normalize.css']:
            shutil.copyfile(os.path.join(current_dir, f"../../resources/html_template/{filename}"), os.path.join(out_dir, filename))
        
        # 1. Add raw info
        raw_html = process_raw(final_edf)
        append_to_template(out_html, "raw-info", raw_html, "raw-info")
        del raw_html

        # 2. Add BIDS Val report
        val_html = process_val_txt(val_txt)
        append_to_template(out_html, "validator-results", val_html, "validator-results")
        del val_html

        # 3. Add extra validation report
        extra_val_html = process_extraval(extraval_txt)
        append_to_template(out_html, "extra-validator-results", extra_val_html, "extra-validator-results")
        del extra_val_html
        
        # No run message
        norun_html = """<ul>
        <li>Rule was not called; therefore, no results will be shown.</li>
        </ul>
        """

        # 4. Epoching
        if len(snakemake.input.epoch_inputs)>0:
            df_epoch = snakemake.input.epoch_inputs[0]
            df = pd.read_csv(df_epoch, sep='\t')
            html_table = df.to_html(index=False, classes='styled-table', justify='center')
            append_to_template(out_html, "events-epoch", html_table, "table-container")
            del df, html_table
        else:
            append_to_template(out_html, "norun-epoch", norun_html, "norun-epoch")
        
        # 5. DOWNSAMPLING REPORT
        if len(snakemake.input.dn_inputs)>0:
            dn_df, edf_dn, edf_prev = snakemake.input.dn_inputs
            # 5.1. Sampling rates
            df = pd.read_csv(dn_df, sep='\t')
            html_table = df.to_html(index=False, classes='styled-table', justify='center')
            append_to_template(out_html, "srate-downsample", html_table, "table-container")
            # 5.2 Size comparison
            epoch_size = os.path.getsize(edf_prev)/(10**6) # MB
            dn_size = os.path.getsize(edf_dn)/(10**6) # MB
            size_comp = {
                'Original size (MB)': [np.round(epoch_size,2)],
                'Downsampled size (MB)': [np.round(dn_size,2)],
                'Compression rate': [np.round(epoch_size/dn_size,2)]
            }
            df = pd.DataFrame(size_comp)
            html_table = df.to_html(index=False, classes='styled-table', justify='center')
            append_to_template(out_html, "size-downsample", html_table, "table-container")
            # 5.3 PSD comparison
            df_orig = get_PSD(edf_prev)
            df_new = get_PSD(edf_dn)
            html_code = plot_PSD_channels(df_orig, df_new, ["Original","Downsampled"])
            append_to_template(out_html, "psd-downsample", html_code, "psd-downsample")
            del df, html_table, df_orig, df_new, html_code
        else:
            append_to_template(out_html, "norun-downsample", norun_html, "norun-downsample")

       
        # 6. DETRENDING
        if len(snakemake.input.detrend_inputs)>0:
            df_detrend, json_detrend, edf_clean_detrend, edf_prev = snakemake.input.detrend_inputs
            # 6.1. General info
            html_code = process_json(json_detrend)
            append_to_template(out_html, "info-detrend", html_code, "info-detrend")
            df = pd.read_csv(df_detrend, sep='\t')
            html_table = df.to_html(index=False, classes='styled-table', justify='center')
            append_to_template(out_html, "mean-detrend", html_table, "table-container")
            # 6.2. Filter figures
            if os.path.exists(os.path.join(os.path.dirname(edf_clean_detrend), 'filt_reponse.png')):
                # Move images to report folder
                os.makedirs(os.path.join(out_dir, snakemake.params.clip), exist_ok=True)
                fig_paths = []
                for fig_path in ['filt_reponse.png', 'filt_freq.png']:
                    new_path = os.path.join(out_dir, snakemake.params.clip, fig_path)
                    shutil.copy(os.path.join(os.path.dirname(edf_clean_detrend), fig_path), new_path)
                    fig_paths.append(os.path.join(snakemake.params.clip, fig_path))
                # Generate html
                figures = zip(fig_paths,
                              ["Filter time response", "Filter frequency response"])
                html_code = get_filt_figures(figures)
                append_to_template(out_html, "filters-detrend", html_code, "filters-detrend")
            # 6.3. PSD comparison
            df_orig = get_PSD(edf_prev)
            df_new = get_PSD(edf_clean_detrend)
            html_code = plot_PSD_channels(df_orig, df_new, ["Original", "Detrended"])
            append_to_template(out_html, "psd-detrend", html_code, "psd-detrend")
            del html_code, df, html_table, df_orig, df_new
        else:
            append_to_template(out_html, "norun-detrend", norun_html, "norun-detrend")


        # 7. REREFERENCING
        if len(snakemake.input.reref_inputs)>0:
            df_reref, json_reref = snakemake.input.reref_inputs
            # 7.1 Discarded channels
            html_code = process_json(json_reref)
            append_to_template(out_html, "info-reref", html_code, "info-reref")
            # 7.2. Table with bipolar comb
            df = pd.read_csv(df_reref, sep='\t')
            html_table = df.to_html(index=False, classes='styled-table', justify='center')
            append_to_template(out_html, "channels-reref", html_table, "table-container")
            del df, html_code, html_table
        else:
            append_to_template(out_html, "norun-reref", norun_html, "norun-reref")


        # 8. PLI Reject
        if len(snakemake.input.PLI_inputs)>0:
            json_pli, edf_clean_PLI, edf_prev = snakemake.input.PLI_inputs
            # 8.1 General info
            html_code = process_json(json_pli)
            append_to_template(out_html, "info-pli", html_code, "info-pli")
            # 8.2. PSD comparison
            df_orig = get_PSD(edf_prev)
            df_new = get_PSD(edf_clean_PLI)
            html_code = plot_PSD_channels(df_orig, df_new, ["Original", "Filtered"])
            append_to_template(out_html, "psd-pli", html_code, "psd-pli")
            del html_code, df_orig, df_new
        else:
            append_to_template(out_html, "norun-pli", norun_html, "norun-pli")

        # 9. Region ID 
        if len(snakemake.input.regionID_inputs)>0:
            df_regionID, json_regionID = snakemake.input.regionID_inputs
            # Discarded channels
            html_code = process_json(json_regionID)
            append_to_template(out_html, "info-regions", html_code, "info-regions")
            # Regions report
            df = pd.read_csv(df_regionID, sep='\t')
            html_table = df.to_html(index=False, classes='styled-table', justify='center')
            append_to_template(out_html, "details-regions", html_table, "table-container")
            del html_code, df, html_table
        else:
            append_to_template(out_html, "norun-regions", norun_html, "norun-regions")
        
        # Add pyprep results
        # pyprep_json = snakemake.input.pyprep_json
        # pyprep_html = snakemake.params.pyprep_html
        # pyprep_results = process_channel_metrics(pyprep_json, pyprep_html)
        # report.add_html(title="Noisy channels based on Pyprep", html=pyprep_results)
        # Write html report
    except:
        logging.exception("Got exception on main handler")
        raise


if __name__ == "__main__":
    main()
