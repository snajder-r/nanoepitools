import re
from typing import Union, IO

from pathlib import Path
import pandas as pd
import itertools


def range_comparison_dist_fn(x, start, end):
    if x.start > end or x.end < start:
        # Not overlapping:
        return min(abs(x.end - start), abs(x.start - end))
    else:
        # Overlapping:
        return 0

class GFFFeature:
    def __init__(self, start, end, type, id, direction, name=None, parent=None):
        self.start = start
        self.end = end
        self.type: str = type
        self.id: str = id
        self.direction = direction
        self.index_width = 100000
        if name is None and id is not None:
            if ":" in id:
                name = id.split(":")[1]
            else:
                name = id
        self.name: str = name
        self.children: Dict[str, GFFFeature] = {}
        self.leafs: List[GFFFeature] = []
        
        self.sorted_children: List[GFFFeature] = []
        self.parent = parent
        
        
        
    def build_sorted(self):
        self.sorted_children = sorted(itertools.chain(self.children.values(), self.leafs), key=lambda x: x.start)
        for child in self.sorted_children:
            child.build_sorted()
        
    def get_in_range(self, start, end, max_recursion=0):
        for feature in self.sorted_children:
            if feature.start > end:
                break
            elif feature.end < start:
                continue
            else:
                yield feature
                if max_recursion > 0:
                    for x in feature.get_in_range(start, end, max_recursion-1):
                        yield x


    def get_nearest_feature(self, start=None, end=None, dist_fn=None, nearest=None, nearest_dist=10e10, max_recursion=0):
        if dist_fn is None:
            dist_fn = lambda x: range_comparison_dist_fn(x, start, end)
            
        for feature in self.sorted_children:
            dist = dist_fn(feature)
            if dist < nearest_dist:
                nearest_dist = dist
                nearest = feature
            else:
                continue
                
        if max_recursion > 0:
            return nearest.get_nearest_feature(dist_fn=dist_fn, max_recursion=max_recursion-1)
        else:
            return nearest, nearest_dist

class GFFAnnotationsReader:

    def __init__(self):
        self.chromosomes: Dict[str, GFFFeature] = {}
        self.included_features = ["gene", "mRNA", "exon"]
        self.index = {}

    def read(self, gff_file, only_protein_coding=True):
        with open(gff_file, "rt") as f:
            
            rowit = (
                {"chr": l[0], "type": l[2], "start": int(l[3]), "end": int(l[4]), "direction": l[6], "info": l[8]}
                for l in (l.split("\t") for l in f if l[0] != "#")
            )

            cur_chrom = None
            for row in rowit:
                if cur_chrom is None or cur_chrom.id != row["chr"]:
                    cur_chrom = GFFFeature(0, 0, type="chrom", direction="+", id=row["chr"])
                    self.chromosomes[row["chr"]] = cur_chrom
                    open_features: Dict[str, GFFFeature] = {}

                if row["type"] not in self.included_features:
                    continue

                info = [kv.split("=") for kv in row["info"].split(";")]
                info = {k: v for k, v in info}
                id = info.get("ID", None)
                name = info.get("Name", None)
                parent_id = info.get("Parent", None)

                if "biotype" in info:
                    if only_protein_coding and not info["biotype"] == "protein_coding":
                        continue

                cur_feature = GFFFeature(row["start"], row["end"], type=row["type"], id=id, direction=row["direction"], name=name)

                if parent_id is not None:
                    if parent_id not in open_features:
                        continue
                    cur_parent = open_features[parent_id]
                else:
                    cur_parent = cur_chrom

                cur_feature.parent = cur_parent

                if id is None:
                    # Only leafs have no id
                    cur_parent.leafs.append(cur_feature)
                else:
                    cur_parent.children[id] = cur_feature
                    open_features[id] = cur_feature

    def build_index(self):
        for chrom_container in self.chromosomes.values():
            chrom_container.build_sorted()

    def _print(self, feature, depth=2):
        for fid, f in feature.children.items():
            print("".join(["-"] * depth), fid)
            self._print(f, depth + 2)
        for f in feature.leafs:
            print("".join(["-"] * depth), f.start, f.end)

    def print(self):
        for chromname, chrom in self.chromosomes.items():
            print(chromname)
            self._print(chrom)
            
