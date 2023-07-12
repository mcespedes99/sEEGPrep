rule run_synthseg:
    input:
        t1=rules.n4_t1_withmask.output.t1,
    output:
        dseg=temp(
            bids(
                root='work',
                datatype="anat",
                **inputs.wildcards['T1w'],
                desc="synthsegnoresample",
                suffix="dseg.nii.gz"
            )
        ),
    container:
        config["singularity"]["graham"]["diffparc"]
    threads: 8
    resources:
        mem_mb=16000,
    group:
        "anat"
    shell:
        "python /opt/SynthSeg/scripts/commands/SynthSeg_predict.py --i {input} --o {output} --cpu --threads {threads}"



rule run_synthseg_withcortparc:
    input:
        t1=rules.n4_t1_withmask.output.t1,
    output:
        dseg=temp(
            bids(
                root='work',
                datatype="anat",
                **inputs.wildcards['T1w'],
                desc="synthsegcortparcnoresample",
                suffix="dseg.nii.gz"
            )
        ),
    container:
        config["singularity"]["graham"]["diffparc"]
    threads: 8
    resources:
        mem_mb=16000,  # right now these are on the high-end -- could implement benchmark rules to do this at some point..
    group:
        "anat"
    shell:
        "python /opt/SynthSeg/scripts/commands/SynthSeg_predict.py --i {input} --o {output} --cpu --threads {threads} --parc"


rule reslice_synthseg_to_t1:
    input:
        t1=bids(
            root='work',
            datatype="anat",
            **inputs.wildcards['T1w'],
            desc="preproc",
            suffix="T1w.nii.gz"
        ),
        dseg=bids(
            root='work',
            datatype="anat",
            **inputs.wildcards['T1w'],
            desc="{dseg}noresample",
            suffix="dseg.nii.gz"
        ),
    output:
        dseg=bids(
            root='work',
            datatype="anat",
            **inputs.wildcards['T1w'],
            desc="{dseg,synthseg|synthsegcortparc}",
            suffix="dseg.nii.gz"
        ),
    container:
        config["singularity"]["graham"]["diffparc"]
    group:
        "anat"
    shell:
        "c3d -interpolation NearestNeighbor {input.t1} {input.dseg} -reslice-identity -o {output.dseg}"

rule run_synthseg_template:
    input:
        t1=os.path.join(workflow.basedir, "..", config["template_t1w"]),
    output:
        dseg=get_template_prefix(
            root='work', subj_wildcards=inputs.wildcards['T1w'], template=config["template"]
        )
        + "_desc-synthseg_dseg.nii.gz",
    container:
        config["singularity"]["graham"]["diffparc"]
    threads: 8
    group:
        "anat"
    shell:
        "python /opt/SynthSeg/scripts/commands/SynthSeg_predict.py --i {input} --o {output} --cpu --threads {threads}"
