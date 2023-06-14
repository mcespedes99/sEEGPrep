rule synthstrip_t1:
    input:
        t1=inputs.path["T1w"],
    output:
        mask=temp(
            bids(
                root='work',
                datatype="anat",
                **inputs.wildcards['T1w'],
                desc="nofixhdrbrain",
                suffix="mask.nii.gz"
            )
        ),
    group:
        "anat"
    container:
        config["singularity"]["diffparc"]
    threads: 8
    shell:
        "python3 /opt/freesurfer/mri_synthstrip -i {input.t1} -m {output.mask} --no-csf"


rule fixheader_synthstrip:
    input:
        t1=inputs.path["T1w"],
        mask=rules.synthstrip_t1.output.mask,
    output:
        mask=bids(
            root='work',
            datatype="anat",
            **inputs.wildcards['T1w'],
            desc="brain",
            suffix="mask.nii.gz"
        ),
    group:
        "anat"
    container:
        config["singularity"]["diffparc"]
    shell:
        "c3d {input.t1} {input.mask} -copy-transform -o {output.mask}"


rule n4_t1_withmask:
    input:
        t1=inputs.path["T1w"],
        mask=rules.fixheader_synthstrip.output.mask,
    output:
        t1=bids(
            root='work',
            datatype="anat",
            **inputs.wildcards['T1w'],
            desc="preproc",
            suffix="T1w.nii.gz"
        ),
    threads: 8
    container:
        config["singularity"]["diffparc"]
    group:
        "anat"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "N4BiasFieldCorrection -d 3 -i {input.t1} -x {input.mask} -o {output}"


rule mask_subject_t1w:
    input:
        t1=rules.n4_t1_withmask.output.t1,
        mask=rules.fixheader_synthstrip.output.mask,
    output:
        t1=bids(
            root='work',
            datatype="anat",
            **inputs.wildcards['T1w'],
            suffix="T1w.nii.gz",
            desc="masked"
        ),
    container:
        config["singularity"]["diffparc"]
    group:
        "anat"
    shell:
        "c3d {input} -multiply -o {output}"