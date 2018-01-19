"""This submodule is the main entry point for vcf-parsing functionality, via
import and the command-line.

"""

import sys
import argparse

from ._extensions import pairwise_diffs, count_hets, get_Tv

from .zip_vcf_fasta import zip_vcf_fasta as zip_fasta

# python -m dnaUtilsPy.vcf (function) [args]
if __name__ == "__main__":
    fun = sys.argv[1]
    parser = argparse.ArgumentParser()

    def parse_args():
        return parser.parse_args(sys.argv[2:])

    if fun == "count_hets":
        parser.add_argument("vcf_name")
        args = parse_args()
        print(count_hets(args.vcf_name).to_csv(sep="\t", index=False), end="")

    elif fun == "get_Tv":
        parser.add_argument("in_vcf")
        parser.add_argument("out_vcf")
        args = parse_args()
        get_Tv(**vars(args))

    else:
        raise ValueError(f"Unrecognized command-line function {fun}")
