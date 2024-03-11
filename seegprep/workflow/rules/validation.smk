rule extra_validation:
    input: 
        edf = inputs.path['edf'],
        channels_tsv = inputs.path['channels'],
        tsv_elec = inputs.path['electrodes'],
        json_file = inputs.path['json'],
    group:
        "subj"
    params:
        modality = 'ieeg'
    output:
        out_txt = bids(
                        root='work',
                        datatype='ieeg',
                        suffix='ieeg.txt',
                        rec='validation',
                        **inputs.wildcards['edf']
                ),
    log:
        bids(
            root='logs',
            suffix='validation.log',
            **inputs.wildcards['edf']
        ),
    script: join(workflow.basedir,'scripts/ext_val.py')

rule bids_validator:
    input:
        bids_dir = config['bids_dir'],
    output:
        out_txt = 'bids_validator_results.txt'
    params:
        container = config['singularity']['graham']['bids_validator'],
        flags = config['bids_validator_flags'],
        mode = 'iEEG'
    shell:
        """
        set +e
        singularity run -e -B {input.bids_dir}:/data {params.container} /data {params.flags} > {output.out_txt}
        exitcode=$?
        if [ $exitcode -gt 1 ]
        then
            exit 1
        else
            exit 0
        fi
        """