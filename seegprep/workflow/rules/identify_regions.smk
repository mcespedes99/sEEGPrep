# Define inputs 
def region_id_inputs():
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
    # Else if regionsID is called but not any of the previous rules
    elif config['regions_id']:
        #print('RegionsID is first')
        return inputs.path['ieeg'], inputs.path['seega_tsv']
    else: # Default: run_all
        #print('reref before regionsID (run all)')
        return rules.rereference.output.out_edf, rules.rereference.output.out_tsv

# Rule
rule identify_regions:
    input:
        edf_tsv = region_id_inputs(),
        # Other parameters
        parc = inputs.path['parc'],
        tf = inputs.path['tf'],
    params:
        reref_run = config['run_all'] or config['rereference'] or run_all
    group:
        "subj"
    output:
        out_edf = bids(
                        root='bids',
                        datatype='ieeg',
                        suffix='clean.edf',
                        **inputs.wildcards['ieeg']
                ),
        out_tsv = bids(
                        root='bids',
                        datatype='ieeg',
                        suffix='regions_native_space.tsv',
                        **inputs.wildcards['ieeg']
                ),
    log:
        bids(
            root='logs',
            suffix='identify_regions.log',
            **inputs.wildcards['ieeg']
        )
    script: join(workflow.basedir,'scripts/identify_regions.py')
