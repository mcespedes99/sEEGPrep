# Define inputs
def input_tsv():
    # If both reref and regions_id are called
    if config['run_all'] or (config['rereference'] and config['regions_id']):
        return rules.rereference.output.out_tsv, rules.identify_regions.output.out_tsv
    # Else if only rereference is called:
    elif config['rereference']:
        return rules.rereference.output.out_tsv
    # If only regions_id is called
    elif config['regions_id']:
        return rules.identify_regions.output.out_tsv
    else: # Default to run_all
        return rules.rereference.output.out_tsv, rules.identify_regions.output.out_tsv

# Define outputs
# If both reref and regions_id are called
out_tsv = []
if config['run_all'] or (config['rereference'] and config['regions_id']):
    # Append reref files
    out_tsv.append(bids(
                        root=out_dir_reref(),
                        datatype='ieeg',
                        suffix='reref_native_space.tsv',
                        rec='reref',
                        **inputs.wildcards['ieeg']
                    )
                )
    # Append regions_id files
    out_tsv.append(bids(
                        root='bids',
                        datatype='ieeg',
                        suffix='regions_native_space.tsv',
                        rec='regionID',
                        **inputs.wildcards['ieeg']
                    )
                )
# Else if only rereference is called:
elif config['rereference']:
    # Append reref files
    out_tsv.append(bids(
                        root=out_dir_reref(),
                        datatype='ieeg',
                        suffix='reref_native_space.tsv',
                        rec='reref',
                        **inputs.wildcards['ieeg']
                    )
                )
# If only regions_id is called
elif config['regions_id']:
    # Append regions_id files
    out_tsv.append(bids(
                        root='bids',
                        datatype='ieeg',
                        suffix='regions_native_space.tsv',
                        rec='regionID',
                        **inputs.wildcards['ieeg']
                    )
                )
else: # Default to run_all
    # Append reref files
    out_tsv.append(bids(
                        root=out_dir_reref(),
                        datatype='ieeg',
                        suffix='reref_native_space.tsv',
                        rec='reref',
                        **inputs.wildcards['ieeg']
                    )
                )
    # Append regions_id files
    out_tsv.append(bids(
                        root='bids',
                        datatype='ieeg',
                        suffix='regions_native_space.tsv',
                        rec='regionID',
                        **inputs.wildcards['ieeg']
                    )
                )

# Rule to merge tsv files (there should be only for each initial edf file, not 1 for clip)
rule merge_tsv:
    input:
        tsv_files = lambda wc: expand(input_tsv(), zip,
                                clip=[f'{number:02}' for number in range(1, number_clips[f'subj_{wc.subject}'][f'ses_{wc.session}']+1)], 
                                allow_missing=True),
    group:
        "merge"
    output:
        out_files = out_tsv,
    benchmark:
       bids(
           root='benchmark',
           suffix='benchmarkMerge.txt',
           **inputs.wildcards['ieeg']
       ),
    log:
        bids(
            root='logs',
            suffix='merge.log',
            **inputs.wildcards['ieeg']
        )   
    script: join(workflow.basedir,'scripts/merge_tsv.py')