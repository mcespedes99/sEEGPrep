# Define output folder
def out_dir_reref():
    # If this is not the last step, place in work
    if config['run_all'] or config['regions_id'] or config['PLI_rej'] or run_all:
        return 'work'
    # If this is the last step, place in bids
    else:
        return 'bids'

# Define inputs 
def reref_inputs():
    # If run_all or filter are called
    if config['run_all'] or config['filter']:
        #print('filter before reref')
        return rules.filter_data.output.out_edf
    # Else if downsample is called
    elif config['downsample']:
        #print('downsample before reref')
        return rules.downsample.output.out_edf
    # Else if reref is executed after epoch extraction
    elif config['epoch']:
        return rules.get_epoch_files.output.out_edf
    # Else if rereference is called but not any of the previous rules
    elif config['rereference']:
        #print('Reref is first')
        return inputs.path['edf']
    else: # Default: run_all
        #print('filter before reref (run all)')
        return rules.filter_data.output.out_edf

# Rule
rule rereference:
    input:
        edf = reref_inputs(),
        # Other parameters
        electrodes_tsv = inputs.path['electrodes'],
    group:
        "subj"
    output:
        out_edf = bids(
                        root=out_dir_reref(),
                        datatype='ieeg',
                        suffix='ieeg.edf',
                        rec='reref',
                        **out_edf_wc
                ),
        out_tsv = bids(
                        root=out_dir_reref(),
                        datatype='ieeg',
                        suffix='reref_native_space.tsv',
                        rec='reref',
                        **out_edf_wc
                ),
        report_df = bids(
                        root='work',
                        datatype='ieeg',
                        suffix='report.tsv',
                        rec='reref',
                        **out_edf_wc
                ),
        report_json = bids(
                        root='work',
                        datatype='ieeg',
                        suffix='report.json',
                        rec='reref',
                        **out_edf_wc
                ),
    resources:
        mem_mb = 16000,
    threads: 16
    benchmark:
       bids(
           root='benchmark',
           suffix='benchmarkReref.txt',
           **out_edf_wc
       ),
    log:
        bids(
            root='logs',
            suffix='rereferencing.log',
            **out_edf_wc
        )
    script: join(workflow.basedir,'scripts/rereference_signal.py')