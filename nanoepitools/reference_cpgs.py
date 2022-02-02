import pyfaidx
from pathlib import Path

def format_cpg(chrom, pos):
    return f"{chrom}:{pos}"


class ReferenceCpGs:
    def __init__(self, reference_fasta: Path):
        self.ref = pyfaidx.Fasta(reference_fasta)
    
    def get_CGs(self, chrom, start, end, upper=True, formatted=False):
        seq = str(self.ref[chrom][start: end + 1])
        if upper:
            seq = seq.upper()
        cg = [format_cpg(chrom, start + i) if formatted else start + i for i, first, second in
            zip(range(len(seq) - 1), seq[:-1], seq[1:]) if first == "C" and second == "G"]
        return cg
    
    def add_cpgs_to_dmr_hits(self, hits, upper=True):
        for hit in hits:
            hit["CpGs"] = self.get_CGs(hit["chrom"], hit["start"], hit["end"], upper=upper)
