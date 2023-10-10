rule mask_template_t1w:
    input:
        t1=os.path.join(workflow.basedir, "..", config["template_t1w"]),
        mask=os.path.join(workflow.basedir, "..", config["template_mask"]),
    output:
        t1=get_template_prefix(
            root='work', subj_wildcards=inputs.wildcards['T1w'], template=config["template"]
        )
        + "_desc-masked_T1w.nii.gz",
    container:
        config["singularity"]["graham"]["autotop"]
    group:
        "anat"
    shell:
        "c3d {input} -multiply -o {output}"


rule greedy_t1_to_template:
    input:
        flo=[
            bids(
                root='work/synthstrip',
                datatype="anat",
                **inputs.wildcards['T1w'],
                suffix="T1w.nii.gz",
                desc="masked",
            ),
            bids(
                root='work',
                datatype="anat",
                **inputs.wildcards['T1w'],
                suffix="dseg.nii.gz",
                desc="synthseg",
            ),
        ],
        ref=[
            get_template_prefix(
                root='work', subj_wildcards=inputs.wildcards['T1w'], template=config["template"]
            )
            + "_desc-masked_T1w.nii.gz",
            get_template_prefix(
                root='work', subj_wildcards=inputs.wildcards['T1w'], template=config["template"]
            )
            + "_desc-synthseg_dseg.nii.gz",
        ],
    params:
        input_fixed_moving=lambda wildcards, input: [
            f"-i {fixed} {moving}" for fixed, moving in zip(input.ref, input.flo)
        ],
        input_moving_warped=lambda wildcards, input, output: [
            f"-rm {moving} {warped}"
            for moving, warped in zip(input.flo, output.warped_flo)
        ],
        affine_iterations="100x50x10",
        fluid_iterations="100x50x10",  #default 100x50x10
        gradient_sigma="1.732vox",  #default 1.732vox
        warp_sigma="0.707vox",  #default 0.707vox
        timestep="1.0",  #default 1.0
    output:
        warp=bids(
            root='work',
            datatype="warps",
            suffix="warp.nii.gz",
            from_="subject",
            to=config["template"],
            **inputs.wildcards['T1w']
        ),
        invwarp=bids(
            root='work',
            datatype="warps",
            suffix="invwarp.nii.gz",
            from_="subject",
            to=config["template"],
            **inputs.wildcards['T1w']
        ),
        warped_flo=[
            bids(
                root='work',
                datatype="anat",
                suffix="T1w.nii.gz",
                space=config["template"],
                desc="greedy",
                **inputs.wildcards['T1w']
            )
        ],
        affine=bids(
            root='work',
            datatype="warps",
            suffix="affine.txt",
            from_="subject",
            to=config["template"],
            desc="itk",
            **inputs.wildcards['T1w']
        ),
        affine_xfm_ras=bids(
            root='work',
            datatype="warps",
            suffix="affine.txt",
            from_="subject",
            to=config["template"],
            desc="ras",
            **inputs.wildcards['T1w']
        ),
    threads: 8
    resources:
        mem_mb=16000,  # right now these are on the high-end -- could implement benchmark rules to do this at some point..
        time=60,  # 1 hrs
    container:
        config["singularity"]["graham"]["autotop"]
    group:
        "anat"
    log:
        bids(
            root="logs",
            suffix="greedy.log",
            template=config["template"],
            **inputs.wildcards['T1w']
        ),
    shell:
        #affine first
        "greedy -d 3 -threads {threads} -a -m NCC 2x2x2 {params.input_fixed_moving} -o {output.affine_xfm_ras} -ia-image-centers -n {params.affine_iterations} &> {log} && "

        "greedy -d 3 -threads {threads} -m NCC 2x2x2 {params.input_fixed_moving} -it {output.affine_xfm_ras} -o {output.warp} -oinv {output.invwarp} -n {params.fluid_iterations} -s {params.gradient_sigma} {params.warp_sigma} -e {params.timestep} &>> {log} && "

        "c3d_affine_tool {output.affine_xfm_ras} -oitk {output.affine} &>> {log} && "

        "greedy -d 3 -threads {threads} -rf {input.ref[0]} {params.input_moving_warped} -r {output.warp} {output.affine_xfm_ras} &>> {log}"
        #then deformable:
        #then convert affine to itk format that ants uses
        #and finally warp the moving image