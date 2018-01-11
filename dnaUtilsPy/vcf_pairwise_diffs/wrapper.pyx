cimport numpy as np
import numpy as np

from pysam.libchtslib cimport htsFile, hts_open, hts_close, bcf_hdr_t, bcf_hdr_read, bcf_hdr_nsamples, bcf_hdr_destroy

cdef extern from "vcfpairdiffs.h":
    void vcfpairdiffs(htsFile *fp, bcf_hdr_t *hdr,
                      int nsmpl, int* ndiff, int* ntot)

def pairwise_diffs(vcf_name):
    vcf_name = vcf_name.encode()
    cdef const char* fname = vcf_name
    cdef htsFile *fp = hts_open(fname, "r")
    cdef bcf_hdr_t *hdr = bcf_hdr_read(fp)

    nsmpl = bcf_hdr_nsamples(hdr)
    samples = [hdr.samples[i].decode() for i in range(nsmpl)]

    cdef int[:,::1] ndiff = np.zeros([nsmpl, nsmpl], dtype=np.intc)
    cdef int[:,::1] ntot = np.zeros([nsmpl, nsmpl], dtype=np.intc)
    vcfpairdiffs(fp, hdr, nsmpl, &ndiff[0,0], &ntot[0,0])

    bcf_hdr_destroy(hdr)
    hts_close(fp)

    return {"samples": samples,
            "n_differences": np.asarray(ndiff),
            "n_comparisons": np.asarray(ntot)}
