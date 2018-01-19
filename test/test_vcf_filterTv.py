import os
import tempfile
import pysam
from dnaUtilsPy.vcf_get_Tv import get_Tv

def _get_Tv(in_vcf, out_vcf):
    bcf_in = pysam.VariantFile(in_vcf)
    bcf_out = pysam.VariantFile(out_vcf, "w", header=bcf_in.header)
    acgt_tuple = tuple("ACGT")
    for rec in bcf_in.fetch():
        if len(rec.alleles) != 2:
            # not biallelic
            continue
        if any(a not in acgt_tuple for a in rec.alleles):
            # not a SNP
            continue
        alleles = "".join(sorted(rec.alleles)).upper()
        if alleles != "AG" and alleles != "CT":
            # not a transition
            bcf_out.write(rec)
    bcf_in.close()
    bcf_out.close()

expected = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20090805
##unused=<XX=AA,Description="Unused generic">
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
##FORMAT=<ID=TS,Number=1,Type=String,Description="Test String">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
20	14370	rs6054257	G	C	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ:TS	0|0:48:1:51,51:String1	1|0:48:8:51,51:SomeOtherString2	1/1:43:5:.,.:YetAnotherString3
20	14370	.	T	A	29	PASS	NS=3;DP=99;AF=0.5;DB;H2	GT:DP:HQ:TS	0|0:9:51,51:String1	1|0:9:51,51:SomeOtherString2	1/1:9:.,.:YetAnotherString3
"""

dired = os.path.dirname(os.path.realpath(__file__))
def test_filter_Tv():
    _, fname = tempfile.mkstemp()
    #_get_Tv(os.path.join(dired, "test.vcf"), fname)
    get_Tv(os.path.join(dired, "test.vcf"), fname)
    with open(fname) as f:
        actual = "".join(f)
    os.remove(fname)
    assert actual == expected
