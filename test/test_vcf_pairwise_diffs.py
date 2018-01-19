import os
import ast
import pprint
import numpy as np
import pysam
from dnaUtilsPy import vcf_pairwise_diffs

# python implementation for prototype/test
def _py_pairwise_diffs(vcf_name):
    bcf_in = pysam.VariantFile(vcf_name)
    samples = bcf_in.header.samples
    n_diffs = np.zeros([len(samples)]*2)
    n_comps = np.zeros([len(samples)]*2)
    for rec in bcf_in.fetch():
        rec_samples = rec.samples
        for i in range(len(samples)):
            for j in range(len(samples)):
                s_i, s_j = rec.samples[i], rec.samples[j]
                idx_i = s_i.allele_indices
                idx_j = s_j.allele_indices
                if idx_i is None or idx_j is None:
                    continue

                for ii in range(len(idx_i)):
                    for jj in range(len(idx_j)):
                        # don't compare to self
                        if i == j and jj == ii:
                            continue
                        a_i = idx_i[ii]
                        a_j = idx_j[jj]
                        if a_i is None or a_j is None:
                            continue
                        n_comps[i,j] += 1
                        if a_i != a_j:
                            n_diffs[i,j] += 1

    return {"samples": list(samples),
            "n_differences": n_diffs,
            "n_comparisons": n_comps}

def test_vcf_pairwise_diffs():
    # test.vcf was copied from htslib/test/test-vcf-api.out
    dired = os.path.dirname(os.path.realpath(__file__))
    pairwise_diffs = vcf_pairwise_diffs.pairwise_diffs(
        os.path.join(dired, "test.vcf"))
    pairwise_diffs2 = _py_pairwise_diffs(
        os.path.join(dired, "test.vcf"))
    for k in ('n_comparisons', 'n_differences'):
        pairwise_diffs[k] = pairwise_diffs[k].tolist()
        pairwise_diffs2[k] = pairwise_diffs2[k].tolist()

    pprint.pprint(pairwise_diffs)

    assert pairwise_diffs == pairwise_diffs2
    with open(os.path.join(dired, "test_vcf_pairwise_diffs.out")) as f:
        assert ast.literal_eval(f.read()) == pairwise_diffs
