# Define output folder for edf file
def out_dir_filt():
    # If this is not the last step, place in work
    if config['run_all'] or config['rereference'] or config['regions_id'] or run_all:
        return 'work'
    # If this is the last step, place in bids
    else:
        return 'bids'

# Define inputs
def filter_inputs():
    # If run_all or downsample are called
    if config['run_all'] or config['downsample']:
        print('downsample before filter')
        return rules.downsample.output.out_edf
    # Else if filter is the first step to execute
    elif config['filter']:
        print('filter is first')
        return inputs.path['ieeg']
    # run_all by default
    else:
        print('default filter')
        return rules.downsample.output.out_edf

print(config['freesurf_dir'])
# Rule
rule filter_data:
    input:
        edf = filter_inputs(),
        # Other parameters
        tsv = inputs.path['seega_tsv'],
        tf = inputs.path['tf'],
    params:
        freesurf_dir = config['freesurf_dir'],
        freesurf_patient = config['freesurf_patient_label']
    group:
        "subj"
    output:
        out_edf = bids(
                        root=out_dir_filt(),
                        datatype='ieeg',
                        suffix='filtered.edf',
                        **inputs.wildcards['ieeg']
                ),
        out_tsv = bids(
                        root='bids',
                        datatype='ieeg',
                        suffix='noisy_data.tsv',
                        **inputs.wildcards['ieeg']
                ),
    log:
        bids(
            root='logs',
            suffix='filtering.log',
            **inputs.wildcards['ieeg']
        )
    script: join(workflow.basedir,'scripts/filter_signal.py')
