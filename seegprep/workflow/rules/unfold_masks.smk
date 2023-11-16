def get_unfolded_masks(wc_values):
    from bids.layout import parse_file_entities
    masks = []
    for surf in config['midthickness_path']:
        entities = parse_file_entities(surf)
        hemi = entities['hemi']
        label = entities['label']
        suffix = entities['suffix']
        mask = bids(
                    root='bids',
                    datatype='anat',
                    label=f'{label}mask',
                    suffix=f'{suffix}.shape.gii',
                    rec='unfold',
                    hemi=hemi,
                    **wc_values
                )
        masks.append(mask)
    # print(segmentations)
    return masks

# Rule to merge tsv files (there should be only for each initial edf file, not 1 for clip)
rule unfold_masks:
    input:
        masks_in = rules.merge_files.output.out_masks,
        inner_surfs = partial(get_surfaces, 'inner_path'),
        midthickness_surfs = partial(get_surfaces, 'midthickness_path'),
        outer_surfs = partial(get_surfaces, 'outer_path'),
    group:
        "merge"
    output:
        out_masks = get_unfolded_masks(inputs.wildcards['ieeg']),
    threads: 8
    benchmark:
       bids(
           root='benchmark',
           suffix='benchmarkUnfold.txt',
           **inputs.wildcards['ieeg']
       ),
    log:
        bids(
            root='logs',
            suffix='unfold.log',
            **inputs.wildcards['ieeg']
        )  
    params:
        container = config["singularity"]["graham"]["autotop"] 
    script: join(workflow.basedir,'scripts/unfold_masks.py')