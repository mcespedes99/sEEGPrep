import os
# Define inputs 
def region_id_inputs():
    # If rereference and PLI_rej are called 
    if config['run_all'] or config['PLI_rej']:
        return rules.PLI_reject.output.out_edf
    # If only reref is called (without PLI_rej)
    elif config['rereference']:
        return rules.rereference.output.out_edf
    # Else if filter is called
    elif config['filter']:
        #print('filter before regionsID')
        return rules.filter_data.output.out_edf
    # Else if downsample is called
    elif config['downsample']:
        #print('filter before regionsID')
        return rules.downsample.output.out_edf
    # Else if filter is executed after epoch extraction
    elif config['epoch']:
        return rules.get_epoch_files.output.out_edf
    # Else if regionsID is called but not any of the previous rules
    elif config['regions_id']:
        #print('RegionsID is first')
        return inputs.path['edf']
    else: # Default: run_all
        #print('reref before regionsID (run all)')
        return rules.PLI_reject.output.out_edf

# Rule
rule identify_regions:
    input:
        edf = region_id_inputs(),
        electrodes_tsv = inputs.path['electrodes'],
        parc = inputs.path['parc'],
        # tfm = rules.transform_7T_to_clinical.output.tfm, # need to create new rule to go from 7T to 1.5T
        colortable = config['colortable']
    params:
        reference_edf = 'bipolar' if config['regions_id'] else config['reference_edf'],
        out_json = bids(
                        root='work',
                        datatype='ieeg',
                        suffix='regions.json',
                        rec='regionID',
                        **out_edf_wc
                ),
        out_mask = bids(
                        root='work',
                        datatype='ieeg',
                        suffix='mask.nii.gz',
                        rec='regionID',
                        **out_edf_wc
                ),
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
        report_df = bids(
                        root='work',
                        datatype='ieeg',
                        suffix='report.tsv',
                        rec='regionID',
                        **out_edf_wc
                ),
        report_json = bids(
                        root='work',
                        datatype='ieeg',
                        suffix='report.json',
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

