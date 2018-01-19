#include <ctype.h>
#include <string.h>
#include "dna_utils.h"

void vcfGetTv(const char *in_vcf, const char *out_vcf){
  htsFile *in_fp = hts_open(in_vcf, "r");
  bcf_hdr_t *hdr_in = bcf_hdr_read(in_fp);
  bcf1_t *rec    = bcf_init1();

  htsFile *out_fp = hts_open(out_vcf, "w");
  bcf_hdr_t *hdr_out = bcf_hdr_dup(hdr_in);
  bcf_hdr_write(out_fp, hdr_out);

  int i;
  char a;
  char alleles[] = "ZZ";
  while ( bcf_read1(in_fp, hdr_in, rec)>=0 )
  {
    if (!bcf_is_snp(rec)) continue;
    if (rec->n_allele!=2) continue;
    for (i=0; i<2; i++) {
      a = toupper(rec->d.allele[i][0]);
      if (a!='A' && a!='C' && a !='G' && a!='T') {
        a = 'Z';
        break;
      }
      alleles[i] = a;
    }
    if (a=='Z') continue;
    if (alleles[0] > alleles[1]){
      a = alleles[1];
      alleles[1] = alleles[0];
      alleles[0] = a;
    }
    if (strncmp(alleles,"AG",2)==0 || strncmp(alleles,"CT",2)==0)
      continue;

    bcf_write1(out_fp, hdr_out, rec);
  }
  bcf_destroy1(rec);
  bcf_hdr_destroy(hdr_in);
  bcf_hdr_destroy(hdr_out);
  hts_close(in_fp);
  hts_close(out_fp);
}
