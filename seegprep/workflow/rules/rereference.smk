rule rereference:
    input: 
        edf = bids(
                    root='derivatives',
                    datatype='ieeg',
                    suffix='filtered.edf',
                    **inputs.wildcards['ieeg']
                  ),
        tsv = inputs.path['seega_tsv'],
    group:
        "subj"
    output:
        out_edf = bids(
                        root='derivatives',
                        datatype='ieeg',
                        suffix='rereferenced.edf',
                        **inputs.wildcards['ieeg']
                  ),
        out_tsv = bids(
                        root='derivatives',
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