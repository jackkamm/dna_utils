import argparse
from ._vcf_get_Tv import get_Tv

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("in_vcf")
    parser.add_argument("out_vcf")
    args = parser.parse_args()
    get_Tv(**vars(args))
