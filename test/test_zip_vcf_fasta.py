import os
from dnaUtilsPy.zip_vcf_fasta import zip_vcf_fasta

outdir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                      "test_zip_vcf_fasta.d")

expected1 = """##fileformat=VCFv4.2
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
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003	new1	new2
20	1	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT	0|0	1|0	1/1	1	.
20	2	.	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT	0|0	1|0	1/1	.	.
20	3	.	G	A	29	PASS	NS=3;DP=99;AF=0.5;DB;H2	GT	0|0	1|0	1/1	0	.
20	4	.	A	G,T	67	.	NS=2;DP=10;AF=0.333,.;AA=T;DB	GT	2	1	./.	2	.
20	5	.	A	G,T	67	.	NS=2;DP=10;AF=0.333,.;AA=T;DB	GT	2	1	./.	.	.
20	6	.	G	A,T	67	.	NS=2;DP=99;AF=0.333,.;AA=T;DB	GT	2	1	./.	.	0
20	7	.	G	A,T	67	.	NS=2;DP=99;AF=0.333,.;AA=T;DB	GT	2	1	./.	.	.
20	8	.	G	A	29	PASS	NS=3;DP=99;AF=0.5;DB;H2	GT	0|0	1|0	1/1	.	1
"""

def test_zip_vcf_fasta():
    outfile = os.path.join(outdir, "out.vcf")
    if os.path.exists(outfile):
        os.remove(outfile)
    zip_vcf_fasta(os.path.join(outdir, "test.vcf"),
                  outfile,
                  [os.path.join(outdir, "fa1.fa"),
                   os.path.join(outdir, "fa2.fa")],
                  ["new1", "new2"])
    assert "".join(open(outfile)) == expected1

expected2 = """##fileformat=VCFv4.2
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
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003	new1	new2
20	1	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT	0|0	1|0	1/1	1	.
20	2	.	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT	0|0	1|0	1/1	.	0
20	3	.	G	A	29	PASS	NS=3;DP=99;AF=0.5;DB;H2	GT	0|0	1|0	1/1	0	.
20	4	.	A	G,T	67	.	NS=2;DP=10;AF=0.333,.;AA=T;DB	GT	2	1	./.	2	0
20	5	.	A	G,T	67	.	NS=2;DP=10;AF=0.333,.;AA=T;DB	GT	2	1	./.	0	2
20	6	.	G	A,T	67	.	NS=2;DP=99;AF=0.333,.;AA=T;DB	GT	2	1	./.	.	0
20	7	.	G	A,T	67	.	NS=2;DP=99;AF=0.333,.;AA=T;DB	GT	2	1	./.	0	.
20	8	.	G	A	29	PASS	NS=3;DP=99;AF=0.5;DB;H2	GT	0|0	1|0	1/1	.	1
"""

def test_zip_vcf_fasta2():
    outfile = os.path.join(outdir, "out2.vcf")
    if os.path.exists(outfile):
        os.remove(outfile)
    zip_vcf_fasta(os.path.join(outdir, "test.vcf"),
                  outfile,
                  [os.path.join(outdir, "fa1.fa"),
                   os.path.join(outdir, "fa2.fa")],
                  ["new1", "new2"], skip_lower=False)
    assert "".join(open(outfile)) == expected2
