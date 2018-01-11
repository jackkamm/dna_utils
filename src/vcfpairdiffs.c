#include "vcfpairdiffs.h"

void vcfpairdiffs(htsFile *fp, bcf_hdr_t *hdr,
                    int nsmpl, int* ndiff, int* ntot) {
  int i, j, ii, jj, iAllele, jAllele, ngt;
  int32_t *iptr, *jptr, ismpl, jsmpl, *gt_arr = NULL, ngt_arr = 0;
  bcf1_t *line = bcf_init();
  while (bcf_read(fp, hdr, line) == 0)
    {
      // note this uses realloc, so free gt_arr after the loop
      ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
      if ( ngt <= 0) continue; // GT not present
      int max_ploidy = ngt/nsmpl;
      for (i=0; i<nsmpl; i++) {
        iptr = gt_arr + i*max_ploidy;
        for (ii=0; ii<max_ploidy; ii++){
          ismpl = iptr[ii];
          // smaller ploidy
          if (ismpl == bcf_int32_vector_end) continue;
          // missing allele
          if (bcf_gt_is_missing(ismpl)) continue;

          iAllele = bcf_gt_allele(ismpl);

          for (j=0; j<=i; j++) {
            jptr = gt_arr + j*max_ploidy;
            for (jj=0; jj<max_ploidy; jj++){
              // don't compare to self
              if (i == j && ii >= jj) continue;

              jsmpl = jptr[jj];
              // smaller ploidy
              if (jsmpl == bcf_int32_vector_end) continue;
              // missing allele
              if (bcf_gt_is_missing(jsmpl)) continue;

              jAllele = bcf_gt_allele(jsmpl);

              ntot[i*nsmpl + j]++;
              if (iAllele != jAllele){
                ndiff[i*nsmpl + j]++;
              }
            }
          }
        }
      }
    }
  // symmetrize
  for (i=0; i<nsmpl; i++) {
    for (j=0; j<=i; j++) {
      ndiff[j*nsmpl+i] += ndiff[i*nsmpl+j];
      ntot[j*nsmpl+i] += ntot[i*nsmpl+j];
    }
  }
  free(gt_arr);
  bcf_destroy(line);
}
