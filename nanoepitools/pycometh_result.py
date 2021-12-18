import gzip
import numpy as np
import tqdm

from nanoepitools.annotations.annotations import GFFAnnotationsReader, GFFFeature
from nanoepitools.annotations.enhancers import Enhancers
import nanoepiseg.math as net_m

from mb_analysis.ase_asm_analysis.mapping_utils import MapToPromoter
from nanoepitools.ranges import intersect_ranges_by_chromosome, intersect_ranges


def merge_duplicate_diffmet_hits(oldhits):
    # Merge so we don't count them double
    newhits = []
    a = 0
    b = 0
    for i, hiti in enumerate(oldhits):
        a += 1
        duplicate = -1
        for j, hitj in enumerate(newhits):
            if hiti["start"] < hitj["end"] and hitj["start"] < hiti["end"] and hiti["chrom"] == hitj["chrom"]:
                duplicate = j
                break
        if duplicate >= 0:
            newhits[duplicate] = {
                "chrom": hiti["chrom"],
                "start": min(hiti["start"], newhits[duplicate]["start"]),
                "end": max(hiti["end"], newhits[duplicate]["end"]),
                "diff": np.mean([hiti["diff"], newhits[duplicate]["diff"]]),
            }
        else:
            b += 1
            newhits.append(hiti)
    return newhits

class PycomethOutput:
    def __init__(self, met_comp_file):
        self.met_comp_file = met_comp_file
    
    def read_pvals(self):
        with gzip.open(self.met_comp_file, "rt") if self.met_comp_file.endswith(".gz") else open(
            self.met_comp_file, "rt"
        ) as met_comp_fp:
            met_comp_header = met_comp_fp.readline().strip().split("\t")
            for i, line in enumerate(met_comp_fp):
                line = {k: v for k, v in zip(met_comp_header, line.strip().split("\t"))}
                yield float(line["pvalue"])
    
    def read_file(
        self,
        drop_insignificant=True,
        b_minus_a=False,
        retest_fun=None,
        pval_threshold=0.05,
        min_diff=0.25,
        progress=False,
    ):
        N = 0
        if progress:
            with gzip.open(self.met_comp_file, "rt") if self.met_comp_file.endswith(".gz") else open(
                self.met_comp_file, "rt"
            ) as met_comp_fp:
                N = len(met_comp_fp.readlines())
        
        with gzip.open(self.met_comp_file, "rt") if self.met_comp_file.endswith(".gz") else open(
            self.met_comp_file, "rt"
        ) as met_comp_fp:
            met_comp_header = met_comp_fp.readline().strip().split("\t")
            with tqdm.tqdm(disable=(N == 0), total=N) as pbar:
                for i, line in enumerate(met_comp_fp):
                    pbar.update(1)
                    line = {k: v for k, v in zip(met_comp_header, line.strip().split("\t"))}
                    if drop_insignificant and "Significant" not in line["comment"]:
                        continue
                    
                    if retest_fun is not None:
                        recomputed_pval = retest_fun(line["raw_llr_list"])
                        if recomputed_pval > pval_threshold:
                            continue
                    
                    line["adj_pvalue"] = float(line["adj_pvalue"])
                    if line["adj_pvalue"] > pval_threshold:
                        continue
                    
                    line["start"] = int(line["start"])
                    line["end"] = int(line["end"])
                    
                    line["diff"] = eval(line["difference"])
                    if len(line["diff"]) == 0:
                        continue
                    line["diff"] = line["diff"][0]
                    if abs(line["diff"]) < min_diff:
                        continue
                    if not b_minus_a:
                        line["diff"] = -line["diff"]
                    yield line
    
    def load_promoters_hit(self, gff: GFFAnnotationsReader, promoter_before_tss, promoter_after_tss, **kwargs):
        diff_met_table = {}
        map_to_promoter = MapToPromoter(
            gff,
            promoter_before_tss=promoter_before_tss,
            promoter_after_tss=promoter_after_tss,
            input_will_be_sorted=True,
        )
        import tqdm
        lines = list(self.read_file(**kwargs))
        for line in tqdm.tqdm(lines):
            # find promoter
            genes_recorded = set()
            for overlapping_promoter in map_to_promoter(line["chromosome"], line["start"], line["end"]):
                # From transcript to gene
                gene = overlapping_promoter.parent
                # Record every gene just once (if multiple promoters per gene are hit)
                if gene.id in genes_recorded:
                    continue
                genes_recorded.update({gene.id})
                
                diff_met_list = diff_met_table.get(gene.sanitized_id(), [])
                diff_met_entry = {}
                diff_met_entry["chrom"] = line["chromosome"]
                diff_met_entry["start"] = line["start"]
                diff_met_entry["end"] = line["end"]
                diff_met_entry["diffmet"] = line["diff"]
                diff_met_entry["gene_name"] = gene.name
                diff_met_list.append(diff_met_entry)
                diff_met_table[gene.sanitized_id()] = diff_met_list
        return diff_met_table
    
    def load_gene_bodies_hit(self, gff: GFFAnnotationsReader, **kwargs):
        diff_met_table = {}
        for line in self.read_file(**kwargs):
            gff_chrom: GFFFeature = gff.chromosomes[line["chromosome"]]
            # find promoter
            genes = list(gff_chrom.get_in_range(line["start"], line["end"], max_recursion=0))
            for gene in genes:
                diff_met_list = diff_met_table.get(gene.id, [])
                diff_met_entry = {}
                diff_met_entry["chrom"] = line["chromosome"]
                diff_met_entry["start"] = line["start"]
                diff_met_entry["end"] = line["end"]
                diff_met_entry["diffmet"] = line["diff"]
                diff_met_entry["gene_name"] = gene.name
                diff_met_list.append(diff_met_entry)
                diff_met_table[gene.id] = diff_met_list
        return diff_met_table
    
    def load_enhancers_hit(self, enhancers: Enhancers, **kwargs):
        diff_met_table = {}
        for line in self.read_file(**kwargs):
            enhancers_hit = enhancers.enhancers_df.loc[
                (enhancers.enhancers_df["chr"] == line["chromosome"])
                & (enhancers.enhancers_df["start"] < line["end"])
                & (line["start"] < enhancers.enhancers_df["end"])
            ]
            for _, enhancers_row in enhancers_hit.iterrows():
                gene = enhancers_row["nearest_gene"]
                
                diff_met_list = diff_met_table.get(gene.sanitized_id(), [])
                diff_met_entry = {}
                diff_met_entry["chrom"] = line["chromosome"]
                diff_met_entry["start"] = line["start"]
                diff_met_entry["end"] = line["end"]
                diff_met_entry["diffmet"] = line["diff"]
                diff_met_entry["gene_name"] = gene.name
                diff_met_list.append(diff_met_entry)
                diff_met_table[gene.sanitized_id()] = diff_met_list
        return diff_met_table
