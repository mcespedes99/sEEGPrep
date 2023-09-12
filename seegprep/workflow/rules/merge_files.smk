# Rule to merge tsv files (there should be only for each initial edf file, not 1 for clip)
# print(number_clips)
rule merge_files:
    input:
        json_files = lambda wc: expand(rules.identify_regions.output.out_json, zip,
                                clip=[f'{number:02}' for number in range(1, number_clips[f'subj_{wc.subject}'][f'ses_{wc.session}']+1)], 
                                allow_missing=True),
        masks_in = lambda wc: expand(rules.identify_regions.output.out_masks, zip,
                                clip=[f'{number:02}' for number in range(1, number_clips[f'subj_{wc.subject}'][f'ses_{wc.session}']+1)], 
                                allow_missing=True),
        colormask_in = lambda wc: expand(rules.identify_regions.output.out_colormask, zip,
                                clip=[f'{number:02}' for number in range(1, number_clips[f'subj_{wc.subject}'][f'ses_{wc.session}']+1)], 
                                allow_missing=True),
    group:
        "merge"
    output:
        out_json = bids(
                        root='bids',
                        datatype='ieeg',
                        suffix='regions_native_space.json',
                        rec='regionID',
                        **inputs.wildcards['ieeg']
                    ),
        out_masks = get_masks(inputs.wildcards['ieeg']),
        out_colormask = bids(
                        root='bids',
                        datatype='ieeg',
                        suffix='colormask.tsv',
                        rec='regionID',
                        **inputs.wildcards['ieeg']
                ),
    benchmark:
       bids(
           root='benchmark',
           suffix='benchmarkMergeJson.txt',
           **inputs.wildcards['ieeg']
       ),
    log:
        bids(
            root='logs',
            suffix='merge_json.log',
            **inputs.wildcards['ieeg']
        )   
    script: join(workflow.basedir,'scripts/merge_files.py')