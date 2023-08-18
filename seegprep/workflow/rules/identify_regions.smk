import os
# Define inputs 
def region_id_inputs():
    # If rereference and PLI_rej are called 
    if config['run_all'] or (config['rereference'] and config['PLI_rej']):
        return rules.PLI_reject.output.out_edf, rules.rereference.output.out_tsv
    # If rereference is not called but PLI_rej is called
    elif config['PLI_rej']:
        return rules.PLI_reject.output.out_edf, inputs.path['seega_tsv']
    # If only reref is called (without PLI_rej)
    elif config['rereference']:
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
        return rules.get_epoch_files.output.out_edf, inputs.path['seega_tsv']
    # Else if regionsID is called but not any of the previous rules
    elif config['regions_id']:
        #print('RegionsID is first')
        return inputs.path['ieeg'], inputs.path['seega_tsv']
    else: # Default: run_all
        #print('reref before regionsID (run all)')
        return rules.PLI_reject.output.out_edf, rules.rereference.output.out_tsv

def get_segmentation(wildcards):
    return expand(config['seg_path'], subject=mapping_clinical_to_7T[f"P{wildcards.subject}"])

# Rule
rule identify_regions:
    input:
        edf_tsv = region_id_inputs(),
        parc = get_segmentation,
        tfm = rules.transform_7T_to_clinical.output.tfm, # need to create new rule to go from 7T to 1.5T
    params:
        colortable = os.path.join(workflow.basedir, "..", config['colortable']),
        reref_run = config['run_all'] or config['rereference'] or run_all,
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
                        rec='regionID',
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

