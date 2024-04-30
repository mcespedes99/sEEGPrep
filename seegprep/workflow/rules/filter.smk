# Define output folder for edf file
def out_dir_filt():
    # If this is not the last step, place in work
    if config['run_all'] or config['rereference'] or config['PLI_rej'] or config['regions_id'] or run_all:
        return 'work'
    # If this is the last step, place in bids
    else:
        return 'bids'

# Define inputs
def filter_inputs():
    # If run_all or downsample are called
    if config['run_all'] or config['downsample']:
        return rules.downsample.output.out_edf
    # Else if filter is executed after epoch extraction
    elif config['epoch']:
        return rules.get_epoch_files.output.out_edf
    # Else if filter is the first step to execute
    elif config['filter']:
        return inputs.path['edf']
    # run_all by default
    else:
        return rules.downsample.output.out_edf

# Rule
rule filter_data:
    input:
        edf = filter_inputs(),
        chn_tsv = inputs.path['channels'],
    group:
        "subj"
    output:
        out_edf = bids(
                        root=out_dir_filt(),
                        datatype='ieeg',
                        suffix='ieeg.edf',
                        rec='denoise',
                        **out_edf_wc
                ),
        report_df = bids(
                        root='work',
                        datatype='ieeg',
                        suffix='report.tsv',
                        rec='denoise',
                        **out_edf_wc
                ),
        report_json = bids(
                        root='work',
                        datatype='ieeg',
                        suffix='report.json',
                        rec='denoise',
                        **out_edf_wc
                ),
    resources:
        mem_mb = 16000,
    threads: 16
    benchmark:
       bids(
           root='benchmark',
           suffix='benchmarkFilter.txt',
           **out_edf_wc
       ),
    log:
        bids(
            root='logs',
            suffix='filtering.log',
            **out_edf_wc
        )
    script: join(workflow.basedir,'scripts/filter_signal.py')
