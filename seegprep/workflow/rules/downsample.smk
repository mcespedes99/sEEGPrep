# Define output folder for edf file
def out_dir_dn():
    # If this is not the last step, place in work
    if config['run_all'] or config['filter'] or config['rereference'] or config['regions_id'] or run_all:
        return 'work'
    # If this is the last step, place in bids
    else:
        return 'bids'

print('downsample aqui')
rule downsample:
    input: 
        edf = inputs.path['ieeg'],
    group:
        "subj"
    output:
        out_edf = bids(
                        root=out_dir_dn(),
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