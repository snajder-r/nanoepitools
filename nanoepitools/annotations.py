import re
from typing import Union, IO

from pathlib import Path
import pandas as pd


class GFFAnnotations:
    def __init__(self, gff_file: Union[Path, str]):
        self.gff_file = gff_file
        self.gene_regex = re.compile("ID=gene:.*;Name=(.*);biotype.*")
        self.transcript_regex = re.compile(".*Parent=gene:[^;]*;Name=([^;]*);.*")

    def read_dataframe(self, **kwargs):
        return pd.read_csv(
            self.gff_file,
            sep="\t",
            comment="#",
            usecols=[0, 2, 3, 4, 6, 8],
            header=None,
            names=["chr", "type", "start", "end", "direction", "info"],
            dtype={"chr": str},
            **kwargs
        )

    def read_genes(self, protein_coding=True):
        annotation_df = self.read_dataframe(chunksize=1000)
        gene_df = []
        for chunk in annotation_df:
            chunk = chunk.loc[chunk["type"] == "gene"]
            if protein_coding:
                chunk = chunk.loc[chunk["info"].map(lambda x: "biotype=protein_coding" in x)]
            chunk = chunk.drop("type", axis=1).copy()
            chunk["gene"] = chunk["info"].map(lambda x: self.gene_regex.match(x).group(1))
            chunk = chunk.drop("info", axis=1).copy()
            gene_df.append(chunk)
        gene_df = pd.concat(gene_df)
        return gene_df.set_index("gene")

    def read_promoters(self, before=2000, after=500):
        """
        :param before: number of basepairs to include before TSS
        :param after: number of basepairs to include after TSS
        :return:
        """
        annotation_df = self.read_dataframe(chunksize=1000)
        promoters_df = []
        for chunk in annotation_df:
            chunk = chunk.loc[
                chunk.apply(lambda x: ("mRNA" in x["type"]) and "Parent=gene" in x["info"], axis=1)
            ].copy()
            # Get the transcript name, e.g. "DDX11L1-202
            chunk["transcript"] = chunk["info"].map(lambda x: self.transcript_regex.match(x).group(1))
            # Translate coordinates to area around TSS

            def transcripts_to_promoter(row):
                if row["direction"] == "+":
                    row["end"] = row["start"] + after
                    row["start"] = row["end"] - before
                elif row["direction"] == "-":
                    row["start"] = row["end"] - after
                    row["end"] = row["end"] + before
                else:
                    assert False
                return row

            chunk = chunk.apply(transcripts_to_promoter, axis=1)
            chunk = chunk.drop(["type", "info"], axis=1)
            promoters_df.append(chunk)
        promoters_df = pd.concat(promoters_df)
        return promoters_df.set_index("transcript")

    def read_tss(self):
        """
        :param before: number of basepairs to include before TSS
        :param after: number of basepairs to include after TSS
        :return:
        """
        annotation_df = self.read_dataframe(chunksize=1000)
        tss_df = []
        for chunk in annotation_df:
            chunk = chunk.loc[
                chunk.apply(
                    lambda x: ("RNA" in x["type"] or "transcript" in x["type"]) and "Parent=gene" in x["info"], axis=1
                )
            ].copy()
            # Get the transcript name, e.g. "DDX11L1-202
            chunk["transcript"] = chunk["info"].map(lambda x: self.transcript_regex.match(x).group(1))

            chunk["tss"] = chunk.apply(lambda x: x["end"] if x["direction"] == "-" else x["start"], axis=1)
            chunk = chunk.drop(["type", "info", "start", "end"], axis=1)
            tss_df.append(chunk)
        tss_df = pd.concat(tss_df)
        return tss_df.set_index("transcript")


class GFFTranscripts:
    def __init__(self, ensemble_id):
        self.ensemble_id = ensemble_id
        self.exons: List[List[int]] = []


class GFFFeature:
    def __init__(self, id=None, name=None):
        self.name: str = name
        self.id: str = id
        self.children: Dict[str, GFFFeature] = {}


class GFFAnnotationChromosome:
    def __init__(self, name):
        self.name = name
        self.tlos: Dict[str, GFFFeature] = {}


class GFFAnnotationsObject:
    def __init__(self):
        self.gene_regex = re.compile("ID=gene:(ENSG[^;]*);Name=([^;]*);.*")
        self.transcript_regex = re.compile(".*Parent=gene:([^;]*);Name=([^;]*);.*")
        self.exon_parent_regex = re.compile(".*Parent=transcript:([^;]*);.*")

        self.id_regex = re.compile("ID=([^:]*):([^;]*).*")

        self.chromosomes: Dict[str, GFFAnnotationChromosome] = {}

        self.included_features = ["gene", "mRNA", "exon"]

    def read_genes(self, gff_file, protein_coding=True):
        annotation_df = pd.read_csv(
            gff_file,
            sep="\t",
            comment="#",
            usecols=[0, 2, 3, 4, 6, 8],
            header=None,
            names=["chr", "type", "start", "end", "direction", "info"],
            dtype={"chr": str},
            chunksize=1000,
        )
        gene_df = []
        for chunk in annotation_df:
            for index, row in chunk.iterrows():
                if row["chr"] != cur_chrom.name:
                    cur_chrom = GFFAnnotationChromosome(row["chr"])
                    self.chromosomes[row["chr"]] = cur_chrom

                if row["type"] not in self.included_features:
                    continue

                has_parent = "Parent" in row["info"]
                if has_parent:
                    if cur_tlo is None:
                        continue

                id_mitch = self.id_regex.match(row["info"])
                id = None if id_mitch is None else id_mitch.group(2)
                cur_feature = GFFFeature(id=id)

                if nothas_parent:
                    cur_tlo.children[id] = cur_feature
                else:
                    cur_tlo.children[id] = cur_feature
                    cur_tlo = cur_feature


path = "/hps/research1/birney/users/adrien/misc/analyses/Medulloblastoma_DNA_promethion/DNA_pipeline_2/results/input/annotation/annotation.gff3"
gff = GFFAnnotationsObject()
gff.read_genes(path)
