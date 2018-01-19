cimport numpy as np
import numpy as np

cdef extern from "vcfGetTv.h":
    void vcfGetTv(const char *in_vcf, const char *out_vcf);

def get_Tv(in_vcf, out_vcf):
    vcfGetTv(in_vcf.encode(), out_vcf.encode())
