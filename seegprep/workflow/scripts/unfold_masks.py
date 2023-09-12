from pathlib import Path
import sys
import logging
import os
import glob
import re


def main():
    in_masks = snakemake.input.masks_in
    inner_surfs = snakemake.input.inner_surfs
    midthickness_surfs = snakemake.input.midthickness_surfs
    outer_surfs = snakemake.input.outer_surfs
    out_masks = snakemake.output.out_masks
    container = snakemake.params.container
    # print(in_files)
    # print(out_files)
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        for in_mask, inner, midthickness, outer, out_mask in zip(
            in_masks, inner_surfs, midthickness_surfs, outer_surfs, out_masks
        ):
            # For label map
            command = f"singularity exec {container} wb_command -volume-to-surface-mapping {in_mask} \
                {midthickness} {out_mask} -enclosing"  # -ribbon-constrained {inner} {outer} -gaussian 1
            os.system(command)
            # For channel masks
            # To find in masks
            filename = Path(in_mask)
            suffixes = "".join(filename.suffixes)
            filename_no_suffix = in_mask.replace(suffixes, "")
            chn_in_img_name = filename_no_suffix + "_chn-*" + suffixes
            # Out masks name
            filename = Path(out_mask)
            suffixes = "".join(filename.suffixes)
            filename_no_suffix = out_mask.replace(suffixes, "")
            pattern = r"(_chn-.*).nii.gz"
            for chn_in_mask in glob.glob(chn_in_img_name):
                regex_search = re.search(pattern, chn_in_mask)
                chn_name = regex_search.group(1)
                chn_out_mask = filename_no_suffix + chn_name + suffixes
                command = f"singularity exec {container} wb_command -volume-to-surface-mapping {chn_in_mask} \
                {midthickness} {chn_out_mask} -ribbon-constrained {inner} {outer} -gaussian 1"  #
                os.system(command)
    except:
        logging.exception("Got exception on main handler")
        raise


if __name__ == "__main__":
    main()
