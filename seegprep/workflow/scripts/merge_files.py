from pathlib import Path
import sys
import logging
import os
import re
import glob


def rename_files(rule, in_files, out_files):
    reg = re.compile(rule)
    in_files_tmp = list(filter(reg.search, in_files))
    out_files_tmp = list(filter(reg.search, out_files))
    print(rule)
    print(in_files_tmp)
    print(out_files_tmp)
    hemi_case = rule.find("hemi") >= 0
    if len(out_files_tmp) > 0:
        # Rename the first file
        os.rename(in_files_tmp[0], out_files_tmp[0])
        if hemi_case:
            # Case of masks (rename first files)
            # To find in masks
            filename = Path(in_files_tmp[0])
            suffixes = "".join(filename.suffixes)
            filename_no_suffix = in_files_tmp[0].replace(suffixes, "")
            chn_in_img_name = filename_no_suffix + "_chn-*" + suffixes
            # Out masks name
            filename = Path(out_files_tmp[0])
            suffixes = "".join(filename.suffixes)
            filename_no_suffix = out_files_tmp[0].replace(suffixes, "")
            pattern = r"(_chn-.*).nii.gz"
            for chn_in_mask in glob.glob(chn_in_img_name):
                regex_search = re.search(pattern, chn_in_mask)
                chn_name = regex_search.group(1)
                chn_out_mask = filename_no_suffix + chn_name + suffixes
                os.rename(chn_in_mask, chn_out_mask)
        # Delete the rest
        for file in in_files_tmp[1:]:
            os.remove(file)
            # Case of masks (delete extra chn masks)
            if hemi_case:
                # To find in masks
                filename = Path(file)
                suffixes = "".join(filename.suffixes)
                filename_no_suffix = file.replace(suffixes, "")
                chn_in_img_name = filename_no_suffix + "_chn-*" + suffixes
                for chn_in_mask in glob.glob(chn_in_img_name):
                    os.remove(chn_in_mask)


def main():
    in_files = snakemake.input
    out_files = snakemake.output
    # print(in_files)
    # print(out_files)
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        # Convert to list if str is provided (not sure if needed)
        if type(out_files) == str:
            out_files = [out_files]
        else:  # It's a snakemake.io.Namedlist
            out_files = list(out_files)
        # Manage possible scenarios by separating tsv files into 2
        # First case of region ID
        print(out_files)
        print(type(out_files))
        print(in_files)
        rename_files("regionID_(.*)(json|space.tsv)", in_files, out_files)
        # Case of reref
        rename_files("reref", in_files, out_files)
        # Case of mask
        rename_files("hemi-L(.*)mask.nii.gz", in_files, out_files)
        rename_files("hemi-R(.*)mask.nii.gz", in_files, out_files)
        # Case of colormask
        rename_files("colormask", in_files, out_files)
    except:
        logging.exception("Got exception on main handler")
        raise


if __name__ == "__main__":
    main()
