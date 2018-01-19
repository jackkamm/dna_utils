import os
import numpy as np
import pandas as pd
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
    n_hom = total - n_hets - n_missing
    return pd.DataFrame(
        list(zip(samples, n_hets, n_hom, n_missing)),
        columns=["Samples", "Hets", "Homs", "Missing"])


def test_vcf_count_hets():
    dired = os.path.dirname(os.path.realpath(__file__))
    vcf = os.path.join(dired, "test.vcf")
    het_counts1 = dnaUtilsPy.vcf_count_hets.count_hets(vcf)
    het_counts2 = _py_count_hets(vcf)
    assert het_counts1.equals(het_counts2)
    assert het_counts1.equals(
        pd.DataFrame(
            {
                'Samples': ['NA00001', 'NA00002', 'NA00003'],
                'Hets': [0, 3, 0],
                "Homs": [6, 3, 3],
                'Missing': [0, 0, 3]
            },
            columns=["Samples", "Hets", "Homs", "Missing"]))
