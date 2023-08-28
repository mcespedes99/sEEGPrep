import numpy as np
import logging
import nibabel as nb


def main():
    in_files = snakemake.input.parc
    out_file = snakemake.output.out_seg
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        if not isinstance(in_files, str):
            # Read first file
            seg_1 = nb.load(in_files[0])
            # Get data
            merged_data = seg_1.get_fdata()
            affine = seg_1.affine
            header = seg_1.header
            if len(in_files) > 1:
                # Get data from other files
                for file in in_files[1:]:
                    seg = nb.load(file)
                    merged_data += seg.get_fdata()
            # Write out file
            new_img = nb.Nifti1Image(merged_data, affine, header)
            nb.save(new_img, out_file)
        else:
            # Read file
            seg_1 = nb.load(in_files)
            # Get data
            merged_data = seg_1.get_fdata()
            affine = seg_1.affine
            header = seg_1.header
            # Write out file
            new_img = nb.Nifti1Image(merged_data, affine, header)
            nb.save(new_img, out_file)
    except:
        logging.exception("Got exception on main handler")
        raise


if __name__ == "__main__":
    main()
