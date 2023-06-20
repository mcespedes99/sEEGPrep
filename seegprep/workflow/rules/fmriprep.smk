rule fmriprep:
    input:
        t1 = rules.mask_subject_t1w.output.t1,
        fs_license=os.environ["FS_LICENSE"] if config["fs_license"] == False else config["fs_license"],
    params:
        synthstrip_dir=join(config['output_dir'], bids(root="work", suffix="synthstrip")),
        fmriprep_outdir=join(config['output_dir'], bids(root="derivatives")),
        dataset_description=join(
            workflow.basedir, "../resources/dataset_description.json"
        ),
        work_directory=join(config['output_dir'], bids( root="work",suffix="fmriprep")),
        fmriprep_opts='--anat-only --skip_bids_validation --notrack --skull-strip-t1w skip',
        container = config["singularity"]["graham"]["fmriprep"],
    output:
        aparc_aseg = 'derivatives/freesurfer/sub-{subject}/mri/aparc+aseg.mgz',
    # container:
    #     config["singularity"]["graham"]["fmriprep"]
    group:
        "anat"
    threads: 8
    resources:
        mem_mb=16000,
        time=1440,
    log:
        bids(
            root='logs',
            suffix='fmriprep.log',
            **inputs.wildcards['T1w']
        ),
    shell:
        """
        singularity run --cleanenv \
        -B {params.synthstrip_dir}:/data \
        -B {params.fmriprep_outdir}:/out \
        -B {params.work_directory}:/work \
        -B {input.fs_license}:/opt/freesurfer/license.txt \
        {params.container} /data /out participant \
        --participant_label {wildcards.subject} -w /work \
        {params.fmriprep_opts} --nprocs {threads} --mem_mb {resources.mem_mb} &> {log}
        """


# """
# fmriprep {params.synthstrip_dir} {params.fmriprep_outdir} participant \
# --participant_label {wildcards.subject} -w {params.work_directory} {params.fmriprep_opts} \
# --fs-license-file {input.fs_license} --nprocs {threads} --mem_mb {resources.mem_mb} &> {log}
# """



        
