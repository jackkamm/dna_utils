#include <htslib/hts.h>
#include <htslib/vcf.h>

void vcfcounthets(htsFile *fp, bcf_hdr_t *hdr, int nsmpl, int* nhet, int* nmiss, int* ntot);
