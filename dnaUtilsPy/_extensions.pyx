cimport numpy as np
import numpy as np
import pandas as pd


from pysam.libchtslib cimport (
    htsFile, hts_open, hts_close, bcf_hdr_t, bcf_hdr_read,
    bcf_hdr_nsamples, bcf_hdr_destroy)


cdef extern from "dna_utils.h":
    void vcfGetTv(const char *in_vcf, const char *out_vcf)
    void vcfCountHets(htsFile *fp, bcf_hdr_t *hdr,
                      int nsmpl, int* nhet, int* nmiss, int* ntot)
    void vcfPairDiffs(htsFile *fp, bcf_hdr_t *hdr,
                      int nsmpl, int* ndiff, int* ntot)


def get_Tv(in_vcf, out_vcf):
    vcfGetTv(in_vcf.encode(), out_vcf.encode())


def count_hets(vcf_name):
    vcf_name = vcf_name.encode()
    cdef const char* fname = vcf_name
    cdef htsFile *fp = hts_open(fname, "r")
    cdef bcf_hdr_t *hdr = bcf_hdr_read(fp)

    nsmpl = bcf_hdr_nsamples(hdr)
    samples = [hdr.samples[i].decode() for i in range(nsmpl)]

    cdef int[::1] nhet = np.zeros([nsmpl], dtype=np.intc)
    cdef int[::1] nmiss = np.zeros([nsmpl], dtype=np.intc)
    cdef int ntot = 0
    vcfCountHets(fp, hdr, nsmpl, &nhet[0], &nmiss[0], &ntot)

    bcf_hdr_destroy(hdr)
    hts_close(fp)

    n_hets=np.asarray(nhet)
    n_missing=np.asarray(nmiss)
    n_hom = ntot - n_hets - n_missing
    return pd.DataFrame(list(zip(samples, n_hets, n_hom, n_missing)),
                        columns=["Samples", "Hets", "Homs", "Missing"])


def pairwise_diffs(vcf_name):
    vcf_name = vcf_name.encode()
    cdef const char* fname = vcf_name
    cdef htsFile *fp = hts_open(fname, "r")
    cdef bcf_hdr_t *hdr = bcf_hdr_read(fp)

    nsmpl = bcf_hdr_nsamples(hdr)
    samples = [hdr.samples[i].decode() for i in range(nsmpl)]

    cdef int[:,::1] ndiff = np.zeros([nsmpl, nsmpl], dtype=np.intc)
    cdef int[:,::1] ntot = np.zeros([nsmpl, nsmpl], dtype=np.intc)
    vcfPairDiffs(fp, hdr, nsmpl, &ndiff[0,0], &ntot[0,0])

    bcf_hdr_destroy(hdr)
    hts_close(fp)

    return {"samples": samples,
            "n_differences": np.asarray(ndiff),
            "n_comparisons": np.asarray(ntot)}
