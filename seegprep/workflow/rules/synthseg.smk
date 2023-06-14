def get_template_prefix(root, subj_wildcards, template):
    """creates prefix for template files, including subject/session wildcards
    so that DAGs for each subject/session are kept independent.
        e.g.: sub-001/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym"""

    path_entities = bids(root=root, **subj_wildcards).split("/")[
        :-1
    ]  # leave out the file prefix

    path_entities.append(f"tpl-{template}")  # sub-folder name
    path_entities.append(f"tpl-{template}")  # file prefix

    return "/".join(path_entities)

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
        config["singularity"]["diffparc"]
    threads: 8
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
        config["singularity"]["diffparc"]
    threads: 8
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
        config["singularity"]["diffparc"]
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
        config["singularity"]["diffparc"]
    threads: 8
    group:
        "anat"
    shell:
        "python /opt/SynthSeg/scripts/commands/SynthSeg_predict.py --i {input} --o {output} --cpu --threads {threads}"