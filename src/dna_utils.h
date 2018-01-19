#include <htslib/hts.h>
#include <htslib/vcf.h>

void vcfGetTv(const char *in_vcf, const char *out_vcf);
void vcfCountHets(htsFile *fp, bcf_hdr_t *hdr, int nsmpl, int* nhet, int* nmiss, int* ntot);
void vcfPairDiffs(htsFile *fp, bcf_hdr_t *hdr, int nsmpl, int* ndiff, int* ntot);
