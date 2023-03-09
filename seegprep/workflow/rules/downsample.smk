rule downsample:
    input: 
        edf = inputs.path['ieeg'],
    group:
        "subj"
    output:
        out_edf = bids(
                        root='derivatives',
                        datatype='ieeg',
                        suffix='downsampled.edf',
                        **inputs.wildcards['ieeg']
                  )
    log:
        bids(
            root='logs',
            suffix='downsampling.log',
            **inputs.wildcards['ieeg']
        )
    script: join(workflow.basedir,'scripts/downsample_signal.py')