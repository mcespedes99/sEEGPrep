def get_segmentation(wildcards):
    segmentations = []
    for seg in config['seg_path']:
        segmentations += expand(seg, subject=mapping_clinical_to_7T[f"P{wildcards.subject}"])
    # print(segmentations)
    return segmentations

# Rule
rule merge_labels:
    input:
        parc = get_segmentation,
    group:
        "merge"
    output:
        out_seg = bids(
            root='work',
            datatype="anat",
            **inputs.wildcards['T1w'],
            suffix="subfields_atlas-bigbrain_dseg.nii.gz"
        ),
    threads: 1
    benchmark:
       bids(
           root='benchmark',
           suffix='benchmarkMerge.txt',
           **inputs.wildcards['T1w']
       ),
    log:
        bids(
            root='logs',
            suffix='merge.log',
            **inputs.wildcards['T1w']
        )
    script: join(workflow.basedir,'scripts/merge_labels.py')