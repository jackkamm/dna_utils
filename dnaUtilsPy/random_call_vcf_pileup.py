import random

_minimal_vcf_header = [
    "##fileformat=VCFv4.2\n",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
]

vcf_info_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]


def random_call_vcf_pileup(pileup_stream, template_vcf_stream,
                           sample_name, quality=30, pileup_format="GATK"):
    """
    Uses pileup file to make a random haploid call at every position in
    the template vcf.

    Pileup file and template vcf should consist of a single chromosome.

    Returns an iterator over the lines of a new vcf with a single sample
    consisting of the random calls.
    """
    if pileup_format.upper() != "GATK":
        raise NotImplementedError

    for vcf_line in template_vcf_stream:
        if vcf_line.startswith("##"):
            continue
        elif vcf_line.startswith("#CHROM"):
            for header_line in _minimal_vcf_header:
                yield header_line
            yield "\t".join(vcf_info_fields + [sample_name]) + "\n"
            break
        else:
            raise ValueError("Malformed VCF header")

    def next_vcf_line():
        try:
            return next(template_vcf_stream).split()
        except StopIteration:
            return None

    def next_pileup_line():
        try:
            line = next(pileup_stream)
            if line.startswith("[REDUCE RESULT]"):
                return None
            return line.split()
        except StopIteration:
            return None

    def output_line(chosen):
        return "\t".join(vcf_line[:5] + [".", ".", ".", "GT", chosen]) + "\n"

    pileup_line = next_pileup_line()
    vcf_line = next_vcf_line()

    while vcf_line is not None:
        vcf_chrom, vcf_pos = vcf_line[:2]

        if pileup_line is not None:
            pileup_chrom, pileup_pos = pileup_line[:2]

            if pileup_chrom != vcf_chrom:
                raise NotImplementedError("Pileup and VCF should consist of a single chromosome")

            vcf_pos = int(vcf_pos)
            pileup_pos = int(pileup_pos)
            if pileup_pos > vcf_pos:
                yield output_line(".")
                vcf_line = next_vcf_line()
            elif pileup_pos < vcf_pos:
                pileup_line = next_pileup_line()
            else:
                assert pileup_pos == vcf_pos

                ref = vcf_line[3]
                alts = vcf_line[4].split(",")
                vcf_alleles = [ref] + alts

                pileup_alleles = []
                for a, q in zip(pileup_line[3], pileup_line[4]):
                    if a in vcf_alleles and ord(q)-33 >= quality:
                        pileup_alleles.append(a)

                if pileup_alleles:
                    chosen = str(vcf_alleles.index(random.choice(pileup_alleles)))
                else:
                    chosen = "."
                yield output_line(chosen)

                pileup_line = next_pileup_line()
                vcf_line = next_vcf_line()
        else:
            yield output_line(".")
            vcf_line = next_vcf_line()
