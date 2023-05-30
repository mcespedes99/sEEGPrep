import os
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
    # Else if filter is executed after epoch extraction
    elif config['epoch']:
        return rules.get_epoch_files.output.out_edf
    # Else if regionsID is called but not any of the previous rules
    elif config['regions_id']:
        #print('RegionsID is first')
        return inputs.path['ieeg'], inputs.path['seega_tsv']
    else: # Default: run_all
        #print('reref before regionsID (run all)')
        return rules.rereference.output.out_edf, rules.rereference.output.out_tsv

def define_parc(wildcards, parc_path):
    parc_options = expand(parc_path, extension=['.mgz','.orig.mgz'], **wildcards)
    # print(parc_options)
    parc_cleaned = []
    for parc_file in parc_options:
        if os.path.exists(parc_file):
            parc_cleaned.append(parc_file)
    # print('Clean')
    # print(parc_cleaned)
    if not parc_cleaned:
        raise ValueError('No parcellation file aparc+aseg was found.')
    return parc_cleaned[0]

# Rule
rule identify_regions:
    input:
        edf_tsv = region_id_inputs(),
        # Other parameters
        parc = lambda wc: define_parc(wc, inputs.path['parc']),
        # expand(inputs.path['parc'], zip, extension=['.mgz','.orig.mgz'], allow_missing=True),
        tf = inputs.path['tf'],
    params:
        reref_run = config['run_all'] or config['rereference'] or run_all
    group:
        "subj"
    output:
        out_edf = bids(
                        root='bids',
                        datatype='ieeg',
                        suffix='ieeg.edf',
                        rec='regionID',
                        **out_edf_wc
                ),
        out_tsv = bids(
                        root='bids',
                        datatype='ieeg',
                        suffix='regions_native_space.tsv',
                        **out_edf_wc
                ),
    resources:
        mem_mb = 16000,
    threads: 16
    benchmark:
       bids(
           root='benchmark',
           suffix='benchmarkRegionID.txt',
           **out_edf_wc
       ),
    log:
        bids(
            root='logs',
            suffix='identify_regions.log',
            **out_edf_wc
        )
    script: join(workflow.basedir,'scripts/identify_regions.py')
