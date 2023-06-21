# Define output folder for edf file
def out_dir_filt():
    # If this is not the last step, place in work
    if config['run_all'] or config['rereference'] or config['PLI_rej'] or config['regions_id'] or run_all:
        return 'work'
    # If this is the last step, place in bids
    else:
        return 'bids'

# Define inputs
def filter_inputs():
    # If run_all or downsample are called
    if config['run_all'] or config['downsample']:
        #print('downsample before filter')
        return rules.downsample.output.out_edf
    # Else if filter is executed after epoch extraction
    elif config['epoch']:
        return rules.get_epoch_files.output.out_edf
    # Else if filter is the first step to execute
    elif config['filter']:
        #print('filter is first')
        return inputs.path['ieeg']
    # run_all by default
    else:
        #print('default filter')
        return rules.downsample.output.out_edf

# Get transform file
def get_tf(wildcards):
    if config['t1_to_coords_tfm'] != False:
        tf_path = expand(inputs.path['tf'], **wildcards)[0]
        if os.path.exists(tf_path):
            return tf_path
    return None

#print(config['freesurf_dir'])
# Rule
rule filter_data:
    input:
        edf = filter_inputs(),
        # Other parameters
        tsv = inputs.path['seega_tsv'],
        t1 = 'derivatives/freesurfer/sub-{subject}/mri/T1.mgz',
    params:
        freesurf_dir = config['freesurf_dir'],
        freesurf_patient = config['freesurf_patient_label'],
        tf = get_tf,
    group:
        "subj"
    output:
        out_edf = bids(
                        root=out_dir_filt(),
                        datatype='ieeg',
                        suffix='ieeg.edf',
                        rec='denoise',
                        **out_edf_wc
                ),
        out_tsv = bids(
                        root='bids',
                        datatype='ieeg',
                        suffix='noisy_data.tsv',
                        rec='denoise',
                        **out_edf_wc
                ),
    resources:
        mem_mb = 16000,
    threads: 16
    benchmark:
       bids(
           root='benchmark',
           suffix='benchmarkFilter.txt',
           **out_edf_wc
       ),
    log:
        bids(
            root='logs',
            suffix='filtering.log',
            **out_edf_wc
        )
    script: join(workflow.basedir,'scripts/filter_signal.py')
