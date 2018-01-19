import argparse
from ._vcf_count_hets import count_hets

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf_name")
    args = parser.parse_args()
    print(count_hets(args.vcf_name).to_csv(sep="\t", index=False),
          end="")
