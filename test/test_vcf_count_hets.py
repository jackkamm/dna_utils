import os
import ast
import pprint
import numpy as np
import pysam
import dnaUtilsPy.vcf_count_hets


def _py_count_hets(vcf_name):
    #vcf_name="test.vcf"
    bcf_in = pysam.VariantFile(vcf_name)
    samples = bcf_in.header.samples

    n_hets = np.zeros(len(samples), dtype=int)
    n_missing = np.zeros(len(samples), dtype=int)
    total = 0
    for rec in bcf_in.fetch():
        if "GT" not in rec.format:
            continue
        total += 1
        for i in range(len(samples)):
            s = rec.samples[i]
            ai = s.allele_indices
            ai_set = set(ai)
            if None in ai_set:
                n_missing[i] += 1
                continue
            elif len(ai_set) > 1:
                n_hets[i] += 1
    return {
        "samples": list(samples),
        "n_hets": list(n_hets),
        "n_missing": list(n_missing),
        "total": total,
    }


def test_vcf_count_hets():
    dired = os.path.dirname(os.path.realpath(__file__))
    vcf = os.path.join(dired, "test.vcf")
    het_counts1 = _py_count_hets(vcf)
    het_counts2 = dnaUtilsPy.vcf_count_hets.count_hets(vcf)
    assert het_counts1 == het_counts2
    assert het_counts1 == {
        'samples': ['NA00001', 'NA00002', 'NA00003'],
        'n_hets': [0, 3, 0],
        'n_missing': [0, 0, 3],
        "total": 6
    }
