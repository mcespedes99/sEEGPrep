from pathlib import Path
import sys
import logging
import os
import glob
import re


def main():
    AP_niftis = snakemake.input.AP_niftis
    PD_niftis = snakemake.input.PD_niftis
    midthickness_surfs = snakemake.input.midthickness_surfs
    out_AP_giftis = snakemake.output.out_AP_giftis
    out_PD_giftis = snakemake.output.out_PD_giftis
    container = snakemake.params.container
    LOG_FILENAME = snakemake.log[0]
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        # For AP direction
        for AP_nifti, midthickness, out_AP_gifti in zip(
            AP_niftis, midthickness_surfs, out_AP_giftis
        ):
            command = f"singularity exec {container} wb_command \
                -volume-to-surface-mapping {AP_nifti} {midthickness} {out_AP_gifti} -enclosing"
            os.system(command)
        # For PD direction
        for PD_nifti, midthickness, out_PD_gifti in zip(
            PD_niftis, midthickness_surfs, out_PD_giftis
        ):
            command = f"singularity exec {container} wb_command \
                -volume-to-surface-mapping {PD_nifti} {midthickness} {out_PD_gifti} -enclosing"
            os.system(command)
    except:
        logging.exception("Got exception on main handler")
        raise


if __name__ == "__main__":
    main()
