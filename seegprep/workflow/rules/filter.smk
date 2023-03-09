rule filter_data:
    input: 
        edf = bids(
                    root='derivatives',
                    datatype='ieeg',
                    suffix='downsampled.edf',
                    **inputs.wildcards['ieeg']
                  ),
        tsv = inputs.path['seega_tsv'],
        tf = inputs.path['tf'],
    params:
        freesurf_dir = config['freesurf_dir'],
        freesurf_patient = config['freesurf_patient_label']
    group:
        "subj"
    output:
        out_edf = bids(
                        root='derivatives',
                        datatype='ieeg',
                        suffix='filtered.edf',
                        **inputs.wildcards['ieeg']
                  ),
        out_tsv = bids(
                        root='derivatives',
                        datatype='ieeg',
                        suffix='noisy_data.tsv',
                        **inputs.wildcards['ieeg']
                  ),
    log:
        bids(
            root='logs',
            suffix='filtering.log',
            **inputs.wildcards['ieeg']
        )
    script: join(workflow.basedir,'scripts/filter_signal.py')