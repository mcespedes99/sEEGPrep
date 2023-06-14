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
        # parc = lambda wc: define_parc(wc, inputs.path['parc']),
        parc = bids(
                root='work',
                datatype="anat",
                **inputs.wildcards['T1w'],
                desc="synthsegcortparc",
                suffix="dseg.nii.gz"
        ),
        tmp_file = rules.greedy_t1_to_template.output.warped_flo,
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


# # Rule to merge tsv files (there should be only for each initial edf file, not 1 for clip)
# rule merge_regionsID:
#     input:
#         tsv_files = lambda wildcards: expand(rules.identify_regions.output.out_tsv, zip,
#                                         clip=[f'{number:02}' for number in range(1, number_clips[f'subj_{wildcards.subject}'][f'ses_{wildcards.session}']+1)], 
#                                         allow_missing=True),
#     group:
#         "merge_regions"
#     output:
#         out_tsv = bids(
#                         root='bids',
#                         datatype='ieeg',
#                         suffix='regions_native_space.tsv',
#                         rec='regionID',
#                         **inputs.wildcards['ieeg']
#                 ),
#     benchmark:
#        bids(
#            root='benchmark',
#            suffix='benchmarkMergeRegions.txt',
#            **inputs.wildcards['ieeg']
#        ),
#     log:
#         bids(
#             root='logs',
#             suffix='mergeRegions.log',
#             **inputs.wildcards['ieeg']
#         )   
#     script: join(workflow.basedir,'scripts/merge_tsv.py')
