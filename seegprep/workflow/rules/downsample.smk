# Define output folder for edf file
def out_dir_dn():
    # If this is not the last step, place in work
    if config['run_all'] or config['filter'] or config['rereference'] or config['regions_id'] or run_all:
        return 'work'
    # If this is the last step, place in bids
    else:
        return 'bids'


# Define inputs
def dn_inputs():
    # If run_all or epochs are called
    if config['run_all'] or config['epoch']:
        #print('epochs before downsample')
        return rules.get_epoch_files.output.out_edf
    # Else if downsample is the first step to execute
    elif config['downsample']:
        #print('downsample is first')
        return inputs.path['ieeg']
    # run_all by default
    else:
        #print('default filter')
        return rules.epochs.output.out_edf


#print('downsample aqui')
rule downsample:
    input: 
        edf = dn_inputs(),
    group:
        "subj"
    output:
        out_edf = bids(
                        root=out_dir_dn(),
                        datatype='ieeg',
                        suffix='ieeg.edf',
                        rec='dn',
                        **out_edf_wc
                  )
    resources:
        mem_mb = 16000,
    threads: 16
    benchmark:
       bids(
           root='benchmark',
           suffix='benchmarkDN.txt',
           **out_edf_wc
       ),
    log:
        bids(
            root='logs',
            suffix='downsampling.log',
            **out_edf_wc
        )
    script: join(workflow.basedir,'scripts/downsample_signal.py')