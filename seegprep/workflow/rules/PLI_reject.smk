# Define output folder
def out_dir_PLI():
    # If this is not the last step, place in work
    if config['run_all'] or config['regions_id'] or run_all:
        return 'work'
    # If this is the last step, place in bids
    else:
        return 'bids'

# Define inputs 
def PLI_inputs():
    # If run_all or filter are called
    if config['run_all'] or config['rereference']:
        #print('reref before regionsID')
        return rules.rereference.output.out_edf, rules.rereference.output.out_tsv
    # Else if filter is called
    elif config['filter']:
        #print('filter before regionsID')
        return rules.filter_data.output.out_edf, inputs.path['seega_tsv']
    # Else if downsample is called
    elif config['downsample']:
        #print('filter before regionsID')
        return rules.downsample.output.out_edf, inputs.path['seega_tsv']
    # Else if filter is executed after epoch extraction
    elif config['epoch']:
        # print('epoch before PLI')
        return rules.get_epoch_files.output.out_edf, inputs.path['seega_tsv']
    # Else if regionsID is called but not any of the previous rules
    elif config['PLI_rej']:
        #print('PLI_rej is first')
        return inputs.path['ieeg'], inputs.path['seega_tsv']
    else: # Default: run_all
        #print('reref before regionsID (run all)')
        return rules.rereference.output.out_edf, rules.rereference.output.out_tsv


#print(config['freesurf_dir'])
# Rule
rule PLI_reject:
    input:
        edf_tsv = PLI_inputs(),
    group:
        "subj"
    output:
        out_edf = bids(
                        root=out_dir_PLI(),
                        datatype='ieeg',
                        suffix='ieeg.edf',
                        rec='PLIreject',
                        **out_edf_wc
                ),
    resources:
        mem_mb = 16000,
    threads: 16
    benchmark:
       bids(
           root='benchmark',
           suffix='PLI.txt',
           **out_edf_wc
       ),
    log:
        bids(
            root='logs',
            suffix='PLI.log',
            **out_edf_wc
        )
    script: join(workflow.basedir,'scripts/reject_PLI.py')