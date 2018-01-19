#include "vcfcounthets.h"

void vcfcounthets(htsFile *fp, bcf_hdr_t *hdr, int nsmpl, int* nhet, int* nmiss, int* ntot){
  int i,ii,ngt,max_ploidy;
  int32_t *iptr, *gt_arr = NULL, ngt_arr = 0;
  bcf1_t *line = bcf_init();
  while (bcf_read(fp, hdr, line) == 0)
    {
      // note this uses realloc, so free gt_arr after the loop
      ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
      if ( ngt <= 0) continue; // GT not present
      (*ntot)++;
      max_ploidy = ngt/nsmpl;
      for (i=0; i<nsmpl; i++) {
        iptr = gt_arr + i*max_ploidy;
        for (ii=0; ii<max_ploidy; ii++){
          // smaller ploidy
          if (iptr[ii] == bcf_int32_vector_end) break;
          // missing allele
          if (bcf_gt_is_missing(iptr[ii])) {
            nmiss[i]++;
            break;
          }
          // het
          if (bcf_gt_allele(iptr[0])!=bcf_gt_allele(iptr[ii])) {
            nhet[i]++;
            break;
          }
        }
      }
    }
  free(gt_arr);
  bcf_destroy(line);
}
