import os
import ast
import pprint
from dnaUtilsPy import vcf_pairwise_diffs

def test_vcf_pairwise_diffs():
    # test.vcf was copied from htslib/test/test-vcf-api.out
    dired = os.path.dirname(os.path.realpath(__file__))
    pairwise_diffs = vcf_pairwise_diffs.pairwise_diffs(
        os.path.join(dired, "test.vcf"))
    pairwise_diffs2 = vcf_pairwise_diffs._test_pairwise_diffs(
        os.path.join(dired, "test.vcf"))
    for k in ('n_comparisons', 'n_differences'):
        pairwise_diffs[k] = pairwise_diffs[k].tolist()
        pairwise_diffs2[k] = pairwise_diffs2[k].tolist()

    pprint.pprint(pairwise_diffs)

    assert pairwise_diffs == pairwise_diffs2
    with open(os.path.join(dired, "test_vcf_pairwise_diffs.out")) as f:
        assert ast.literal_eval(f.read()) == pairwise_diffs
