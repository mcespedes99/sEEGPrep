rule identify_regions:
    input: 
        edf = bids(
                    root='derivatives',
                    datatype='ieeg',
                    suffix='rereferenced.edf',
                    **inputs.wildcards['ieeg']
                  ),
        tsv = bids(
                    root='derivatives',
                    datatype='ieeg',
                    suffix='reref_native_space.tsv',
                    **inputs.wildcards['ieeg']
                  ),
        parc = inputs.path['parc'],
        tf = inputs.path['tf'],
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