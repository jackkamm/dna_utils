cimport numpy as np
import numpy as np

from pysam.libchtslib cimport htsFile, hts_open, hts_close, bcf_hdr_t, bcf_hdr_read, bcf_hdr_nsamples, bcf_hdr_destroy

cdef extern from "vcfcounthets.h":
    void vcfcounthets(htsFile *fp, bcf_hdr_t *hdr, int nsmpl, int* nhet, int* nmiss, int* ntot);

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
    vcfcounthets(fp, hdr, nsmpl, &nhet[0], &nmiss[0], &ntot)

    bcf_hdr_destroy(hdr)
    hts_close(fp)

    n_hets=np.asarray(nhet)
    n_missing=np.asarray(nmiss)
    return {"samples": samples,
            "n_hets": n_hets.tolist(),
            "n_missing": n_missing.tolist(),
            "total": ntot}
