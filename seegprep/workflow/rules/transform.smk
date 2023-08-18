def get_7T(wildcards): 
    wc = inputs.wildcards['T1w']
    if 'subject' in wc:
        # print(wc)
        del wc['subject']
    t1_snsx = bids(
            root=config['root_snsx'],
            datatype="anat",
            **wc,
            suffix="T1w.nii.gz",
            acq=config['acq_snsx'],
            run='01',
            subject=mapping_clinical_to_7T[f"P{wildcards.subject}"]
            )
    return t1_snsx
    

rule transform_7T_to_clinical:
    input:
        t1_clinical =inputs.path["T1w"],
        t1_7T = get_7T,
    output:
        tfm=bids(
            root='work',
            datatype="anat",
            **inputs.wildcards['T1w'],
            suffix="tf.txt"
        ),
    group:
        "anat"
    threads: 4
    container:
        config["singularity"]["graham"]["diffparc"]
    shell:
        "greedy -d 3 -a -m NMI -i {input.t1_clinical} {input.t1_7T} -o {output.tfm} -n 100x40x20 -ia-image-centers -dof 6 -threads {threads} "