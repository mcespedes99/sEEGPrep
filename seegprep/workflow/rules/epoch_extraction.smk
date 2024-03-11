import numpy as np
import copy

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

def get_time(wildcards, number_clips):
    identifier = expand(
                inputs.path['events'],
                zip,
                **wildcards
            )
    # print(subj_id, ses_id, number_clips[subj_id][ses_id])
    return number_clips[identifier[0]]*180

# Define output folder for edf file
def out_dir_epochs():
    # If this is not the last step, place in work
    if config['run_all'] or config['downsample'] or config['filter'] or config['rereference'] or config['PLI_rej'] or config['regions_id'] or run_all:
        # print('test')
        return 'work'
    # If this is the last step, place in bids
    else:
        return 'bids'

rule epochs:
    input: 
        edf = inputs.path['edf'],
        annot = inputs.path['events'],
    output:
        tmp_file = temp(
                        bids(
                            root=out_dir_epochs(),
                            datatype='ieeg',
                            suffix='run.txt',
                            **inputs.wildcards['edf']
                            )
                        ),
        report_file = bids(
                        root=out_dir_epochs(),
                        datatype='ieeg',
                        suffix='report.tsv',
                        rec='clip',
                        **inputs.wildcards['edf']
                        )
    # This can be changed using --set-threads epochs=N. It acts as a balance between subject and channel parallelization.
    threads: 16
    group:
        "epochs"
    resources:
        mem_mb = 64000,
        storage_mb = get_edf_mb, # Indicate avail space using --resources storage_mb=750000
        time = lambda wc: get_time(wc, number_clips),
        tmpdir = lambda wc: get_tmpdir()
    benchmark:
       bids(
           root='benchmark',
           suffix='benchmarkEpoch.txt',
           **inputs.wildcards['edf']
       ),
    log:
        bids(
            root='logs',
            suffix='epochs.log',
            **inputs.wildcards['edf']
        ),
    script: join(workflow.basedir,'scripts/epochs.py')



rule get_epoch_files:
    input: 
        rules.epochs.output.tmp_file
    group:
        "subj"
    output:
        out_edf = bids(
                    root=out_dir_epochs(),
                    datatype='ieeg',
                    suffix='ieeg.edf',
                    rec='clip',
                    **out_edf_wc
                    )
    script: join(workflow.basedir,'scripts/rename_epochs.py')