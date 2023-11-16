def get_unfolded_coords(type_coord, wc_values):
    from bids.layout import parse_file_entities
    coords = []
    for nifti in config[type_coord]:
        entities = parse_file_entities(nifti)
        hemi = entities['hemi']
        label = entities['label']
        suffix = entities['suffix']
        desc = entities['desc']
        space = entities['space']
        direction = entities['direction']
        new_nifti = bids(
                    root='work',
                    datatype='coords',
                    label=label,
                    suffix=f'{suffix}unfold.shape.gii',
                    desc = desc,
                    space = space,
                    dir = direction,
                    hemi=hemi,
                    **wc_values
                )
        coords.append(new_nifti)
    # print(coords)
    return coords

def get_coords(type_coord, wildcards):
    coord_files = []
    # print(type_coord)
    for nifti in config[type_coord]:
        coord_files += expand(nifti, subject=mapping_clinical_to_7T[f"P{wildcards.subject}"])
    return coord_files

def get_surfaces(type_surf, wildcards):
    surfaces = []
    for surf in config[type_surf]:
        surfaces += expand(surf, subject=mapping_clinical_to_7T[f"P{wildcards.subject}"])
    # print(segmentations)
    return surfaces

# Rule to merge tsv files (there should be only for each initial edf file, not 1 for clip)
rule coords_to_surf:
    input:
        AP_niftis = partial(get_coords, 'AP_nifti'),
        PD_niftis = partial(get_coords, 'PD_nifti'),
        midthickness_surfs = partial(get_surfaces, 'midthickness_path'),
    group:
        "anat"
    output:
        out_AP_giftis = get_unfolded_coords('AP_nifti', inputs.wildcards['T1w']),
        out_PD_giftis = get_unfolded_coords('PD_nifti', inputs.wildcards['T1w']),
    threads: 8
    benchmark:
       bids(
           root='benchmark',
           suffix='benchmarkCoords.txt',
           **inputs.wildcards['T1w']
       ),
    log:
        bids(
            root='logs',
            suffix='unfoldcoords.log',
            **inputs.wildcards['T1w']
        )  
    params:
        container = config["singularity"]["graham"]["autotop"] 
    script: join(workflow.basedir,'scripts/coords_to_surf.py')