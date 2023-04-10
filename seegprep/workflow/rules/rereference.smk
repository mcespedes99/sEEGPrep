# Define output folder
def out_dir_reref():
    # If this is not the last step, place in work
    if config['run_all'] or config['regions_id'] or run_all:
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
    # Else if rereference is called but not any of the previous rules
    elif config['rereference']:
        #print('Reref is first')
        return inputs.path['ieeg']
    else: # Default: run_all
        #print('filter before reref (run all)')
        return rules.filter_data.output.out_edf

# Rule
rule rereference:
    input:
        edf = reref_inputs(),
        # Other parameters
        tsv = inputs.path['seega_tsv'],
    group:
        "subj"
    output:
        out_edf = bids(
                        root=out_dir_reref(),
                        datatype='ieeg',
                        suffix='rereferenced.edf',
                        **inputs.wildcards['ieeg']
                ),
        out_tsv = bids(
                        root=out_dir_reref(),
                        datatype='ieeg',
                        suffix='reref_native_space.tsv',
                        **inputs.wildcards['ieeg']
                ),
    log:
        bids(
            root='logs',
            suffix='rereferencing.log',
            **inputs.wildcards['ieeg']
        )
    script: join(workflow.basedir,'scripts/rereference_signal.py')