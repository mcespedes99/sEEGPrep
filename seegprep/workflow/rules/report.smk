# Define inputs
def epoch_report():
    # If run_all or downsample are called
    if config['run_all'] or config['epoch'] or run_all:
        #print('downsample before filter')
        return [rules.epochs.output.report_file]
    return []

def dn_report():
    out_list = []
    # If run_all or downsample are called, append the report
    if config['run_all'] or config['downsample'] or run_all:
        #print('downsample before filter')
        out_list += [rules.downsample.output.report_file, rules.downsample.output.out_edf]
        # Define the edf previous to dn 
        # If run_all or epochs are called
        if config['run_all'] or config['epoch']:
            #print('epochs before downsample')
            out_list.append(rules.get_epoch_files.output.out_edf)
        # Else if downsample is the first step to execute
        elif config['downsample']:
            #print('downsample is first')
            out_list.append(inputs.path['edf'])
        # run_all by default
        else:
            #print('default filter')
            out_list.append(rules.epochs.output.out_edf)
    return out_list

def detrend_report():
    out_list = []
    # If run_all or downsample are called
    if config['run_all'] or config['filter'] or run_all:
        #print('downsample before filter')
        out_list += [rules.filter_data.output.report_df, rules.filter_data.output.report_json, rules.filter_data.output.out_edf]
        # Define the edf before this one
        # If run_all or downsample are called
        if config['run_all'] or config['downsample']:
            #print('downsample before filter')
            out_list.append(rules.downsample.output.out_edf)
        # Else if filter is executed after epoch extraction
        elif config['epoch']:
            out_list.append(rules.get_epoch_files.output.out_edf)
        # Else if filter is the first step to execute
        elif config['filter']:
            #print('filter is first')
            out_list.append(inputs.path['edf'])
        # run_all by default
        else:
            #print('default filter')
            out_list.append(rules.downsample.output.out_edf)
    return out_list

def reref_report():
    # If run_all or downsample are called
    if config['run_all'] or config['rereference'] or run_all:
        #print('downsample before filter')
        return [rules.rereference.output.report_df, rules.rereference.output.report_json]
    return []

def PLI_report():
    out_list = []
    # If run_all or downsample are called, append the report
    if config['run_all'] or config['PLI_rej'] or run_all:
        #print('downsample before filter')
        out_list+= [rules.PLI_reject.output.report_json, rules.PLI_reject.output.out_edf]
        # Define the edf previous to PLI
        # If run_all or filter are called
        if config['run_all'] or config['rereference'] or run_all:
            #print('reref before regionsID')
            out_list.append(rules.rereference.output.out_edf)
        # Else if filter is called
        elif config['filter']:
            out_list.append(rules.filter_data.output.out_edf)
        # Else if downsample is called
        elif config['downsample']:
            out_list.append(rules.downsample.output.out_edf)
        # Else if filter is executed after epoch extraction
        elif config['epoch']:
            # print('epoch before PLI')
            out_list.append(rules.get_epoch_files.output.out_edf)
        # Else if PLI is called but not any of the previous rules
        else:
            #print('PLI_rej is first')
            out_list.append(inputs.path['edf'])
    return out_list

def regionID_report():
    # If run_all or downsample are called
    if config['run_all'] or config['regions_id'] or run_all:
        #print('downsample before filter')
        return [rules.identify_regions.output.report_df, rules.identify_regions.output.report_json]
    return []

def get_edf():
    # If regions_id is the last step
    if config['run_all'] or config['regions_id']:
        # Only use last output of the rules.merge_tsv.output.out_files, which is the one under bids
        return rules.identify_regions.output.out_edf
    # Power line rejection is the last step and rereference is called
    elif config['PLI_rej']:
        # print('PLI')
        return rules.PLI_reject.output.out_edf
    # Else if rereference is the last step
    elif config['rereference']:
        # print(rules.rereference.output.out_edf, rules.merge_tsv.output.out_files[-1])
        # print(type(rules.rereference.output.out_edf), type(rules.merge_tsv.output.out_files[-1]))
        return rules.rereference.output.out_edf
    # Else if filtering is the last step
    elif config['filter']:
        return rules.filter_data.output.out_edf
    # Else if downsampling is the last step
    elif config['downsample']:
        return rules.downsample.output.out_edf
    # Else if epoch extraction is the last step
    elif config['epoch']:
        return rules.get_epoch_files.output.out_edf
    # Else: default case run_all
    else:
        # print('hola')
        return rules.identify_regions.output.out_edf

rule get_channels_metrics:
    input: 
        edf = get_edf(),
        channels_tsv = rules.rereference.output.out_tsv if config['run_all'] or config['rereference'] or run_all else inputs.path['channels'],
    group:
        "subj"
    output:
        out_json = bids(
                        root='work',
                        datatype='ieeg',
                        suffix=f'ieeg.json',
                        rec='pyprep',
                        **out_edf_wc
                ),
    resources:
        mem_mb = 16000,
    log:
        bids(
            root='logs',
            suffix='pyprep.log',
            **out_edf_wc
        ),
    script: join(workflow.basedir,'scripts/run_pyprep.py')


rule create_report:
    input: 
        edf = get_edf(),
        channels_tsv = inputs.path['channels'],
        val_txt = rules.bids_validator.output.out_txt,
        extraval_txt = rules.extra_validation.output.out_txt,
        epoch_inputs = epoch_report(),
        dn_inputs = dn_report(),
        detrend_inputs = detrend_report(),
        reref_inputs = reref_report(),
        PLI_inputs = PLI_report(),
        regionID_inputs = regionID_report(),
        # pyprep_json = rules.get_channels_metrics.output.out_json,
    params:
        clip ='{clip}',
        # pyprep_html = bids(
        #                 root='work',
        #                 datatype='ieeg',
        #                 suffix=f'ieeg.html',
        #                 rec='pyprep',
        #                 **out_edf_wc
        #         ),
    group:
        "subj"
    output:
        out_html = bids(
                        root='bids',
                        datatype='report',
                        suffix=f'ieeg.html',
                        **out_edf_wc
                ),
    threads: 4
    resources:
        mem_mb = 16000,
        time = 180,
    log:
        bids(
            root='logs',
            suffix='report.log',
            **out_edf_wc
        ),
    script: join(workflow.basedir,'scripts/create_report.py')