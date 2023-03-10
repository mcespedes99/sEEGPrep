# Define inputs
def define_all_inputs():
    # If regions_id is the last step
    if config['run_all'] or config['regions_id']:
        print('run_all')
        return rules.identify_regions.output
    # Else if rereference is the last step
    elif config['rereference']:
        return rules.rereference.output
    # Else if filtering is the last step
    elif config['filter']:
        return rules.filter_data.output
    # Else if downsampling is the last step
    elif config['downsample']:
        return rules.downsample.output
    # Else: default case run_all
    else:
        print('hola')
        return rules.identify_regions.output

### RULE ALL
# If regions_id is the last step
rule all:
    default_target: True  
    input:
        expand(
            expand(
                define_all_inputs(),
                allow_missing = True,
            ),
            zip,
            **inputs.zip_lists['ieeg']
        ),
