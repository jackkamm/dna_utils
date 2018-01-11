#include <htslib/hts.h>
#include <htslib/vcf.h>

void vcfpairdiffs(htsFile *fp, bcf_hdr_t *hdr, int nsmpl, int* ndiff, int* ntot);
