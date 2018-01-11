import pysam

# TODO rewrite in C because this uses buggy/unstable parts of pysam
def zip_vcf_fasta(in_vcf, out_vcf, fasta_paths, sample_names,
                  skip_lower=True):
    """
    Adds samples from fasta to VCF at positions in the VCF.
    Currently this discards all sample info aside from GT.
    """
    if len(sample_names) != len(fasta_paths):
        raise ValueError(
            "len(sample_names) != len(fasta_paths)")

    fastas = [pysam.FastaFile(path) for path in fasta_paths]

    bcf_in = pysam.VariantFile(in_vcf)

    header = bcf_in.header.copy()
    for name in sample_names:
        header.add_sample(name)
    bcf_out = pysam.VariantFile(out_vcf, "w", header=header)

    for rec in bcf_in.fetch():
        alleles = [a.upper() for a in rec.alleles]
        new_gts = []
        for f in fastas:
            raw_gt = f.fetch(rec.chrom, rec.start, rec.stop)
            if not skip_lower:
                raw_gt = raw_gt.upper()
            try:
                idx = alleles.index(raw_gt)
            except ValueError:
                idx = None
            new_gts.append(idx)
        # FIXME: new_record() arguments experimental; use with caution and expect changes
        new_rec = bcf_out.new_record(
            contig=rec.contig,
            start=rec.start,
            stop=rec.stop,
            alleles=alleles,
            id=rec.id, qual=rec.qual, filter=rec.filter,
            info=rec.info)
        # TODO: add other fields besides GT
        for i, s in enumerate(rec.samples.values()):
            new_rec.samples[i].allele_indices = s.allele_indices
            new_rec.samples[i].phased = s.phased
        for i, gt in enumerate(new_gts, start=len(rec.samples)):
            new_rec.samples[i].allele_indices = (gt,)
        bcf_out.write(new_rec)

    bcf_in.close()
    bcf_out.close()
    for f in fastas:
        f.close()
