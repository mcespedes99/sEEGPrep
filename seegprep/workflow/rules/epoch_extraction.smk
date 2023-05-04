import numpy as np

# Function to get size of edf file in mb
def get_edf_mb(wildcards, input):
    return int(np.ceil(os.path.getsize(input.edf)/10**6))

# Function to get tmpdir
def get_tmpdir():
    if os.environ.get(config['local_scratch_env']) != None:
        return os.environ.get(config['local_scratch_env'])
    elif os.path.exists(config['local_scratch']):
        return config['local_scratch']
    else:
        return None

# Define output folder for edf file
def out_dir_epochs():
    # If this is not the last step, place in work
    if config['run_all'] or config['downsample'] or config['filter'] or config['rereference'] or config['regions_id'] or run_all:
        return 'work'
    # If this is the last step, place in bids
    else:
        return 'bids'

rule epochs:
    input: 
        edf = inputs.path['ieeg'],
        annot = inputs.path['annotations'],
    params:
        script = join(workflow.basedir,'../code/epochs.py'),
        tmpdir = lambda wc: get_tmpdir()
    output:
        out_edf = bids(
                        root=out_dir_epochs(),
                        datatype='ieeg',
                        suffix='ieeg_epoch.edf',
                        **inputs.wildcards['ieeg']
                  )
    # This can be changed using --set-threads epochs=N. It acts as a balance between subject and channel parallelization.
    threads: 16
    group:
        "subj"
    resources:
        mem_mb = 64000,
        storage_mb = get_edf_mb, # Indicate avail space using --resources storage_mb=750000
        time = 180,
    benchmark:
       bids(
           root='benchmark',
           suffix='benchmarkEpoch.txt',
           **inputs.wildcards['ieeg']
       ),
    log:
        bids(
            root='logs',
            suffix='epochs.log',
            **inputs.wildcards['ieeg']
        ),
    script: join(workflow.basedir,'scripts/epochs.py')