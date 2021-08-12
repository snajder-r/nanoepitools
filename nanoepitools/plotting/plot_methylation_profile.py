from typing import Dict
from typing import List, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Patch

from meth5.sparse_matrix import SparseMethylationMatrixContainer
from nanoepitools.math import p_to_llr

def plot_met_profile(
    matrix: np.ndarray,
    samples: np.ndarray = None,
    sample_order: List[str] = None,
    sample_colors: Dict = None,
    site_genomic_pos=None,
    site_genomic_pos_end=None,
    marker_height=0.75,
    segment: np.array = None,
    highlights: List[Tuple] = None,
    highlight_color: Union[str, List[str]] = None,
    min_marker_width_relative: float = 0.002,
    sample_hatch: Dict = None,
    aggregate_samples=False
):
    def val_to_color(val):
        return 1 - np.exp(-np.abs(val) * 0.5)
    
    if samples is None:
        samples = np.array(["_" for _ in range(matrix.shape[0])])
        sample_order = ["_"]
    
    if sample_order is None:
        sample_order = sorted(list(set(samples)))
    
    if site_genomic_pos_end is not None:
        x_range = np.max(site_genomic_pos_end) - np.min(site_genomic_pos)
        min_marker_width = min_marker_width_relative * x_range
    
    y_off = 0
    start = 0
    x_max = 0
    end = matrix.shape[1]
    sample_marker_range = {}
    
    
    for s in sample_order:
        x = np.arange(start, end)
        part_matrix = matrix[:, x][(samples == s)]
        
        if part_matrix.shape[0] <= 1:
            continue
        
        if aggregate_samples:
            part_matrix = (part_matrix > 2.5).sum(axis=0) / (np.abs(part_matrix) > 2.5).sum(axis=0)
            part_matrix = p_to_llr(part_matrix)
            part_matrix[np.isnan(part_matrix)] = 0
            part_matrix = part_matrix[np.newaxis,:]
        
        active_reads = np.array((part_matrix != 0).sum(axis=1)).flatten() > 0
        
        part_matrix = part_matrix[active_reads]
        hasval = np.array(part_matrix != 0).flatten()
        y = np.arange(part_matrix.shape[0]) + y_off
        
        x, y = np.meshgrid(x, y)
        x = x.flatten()[hasval]
        y = y.flatten()[hasval]
        matrix_data = np.array(part_matrix).flatten()[hasval]
        
        color = [[24/255, 3/255, 248/255, val_to_color(-v)] if v < 0 else [252/255, 127/255, 44/255, val_to_color(v)] for v in matrix_data]
        
        if site_genomic_pos is not None:
            if site_genomic_pos_end is not None:
                x_end = site_genomic_pos_end[x]
            x = site_genomic_pos[x]  # Translate to actual pos on chrom
            if site_genomic_pos_end is not None:
                # Makes it so very short blocks are still visible
                marker_adjust = (min_marker_width - x_end + x) // 2
                marker_adjust[marker_adjust < 0] = 0
                x = x - marker_adjust
                x_end = x_end + marker_adjust
        
        if site_genomic_pos_end is None:
            plt.scatter(x, y, c=color, marker="|", s=15 * (0.25 + marker_height))
        else:
            # if the end location for each marker is given, we need to plot
            # rectangles
            
            patches = [Rectangle((x[i], y[i]), x_end[i] - x[i] + 1, marker_height) for i in range(len(x))]
            
            patch_collection = PatchCollection(patches)
            patch_collection.set_color(color)
            patch_collection.set_edgecolor(None)
            plt.gca().add_collection(patch_collection)
            plt.gca().autoscale_view()
        
        if sample_colors is not None:
            sample_marker_range[s] = (y_off, part_matrix.shape[0] + y_off)
            x_max = max(x.max(), x_max)
        
        y_off += part_matrix.shape[0]
    
    if sample_colors is not None:
        xlim = plt.xlim()
        marker_width = (xlim[1] - xlim[0]) * 0.015
        for s, marker in sample_marker_range.items():
            kwargs = {}
            if sample_hatch is not None:
                kwargs.update({"hatch": sample_hatch[s], "edgecolor": "w"})
            
            patch = Rectangle(
                (xlim[1] - marker_width * 2, marker[0]),
                marker_width,
                marker[1] - marker[0],
                facecolor=sample_colors[s],
                **kwargs
            )
            plt.gca().add_patch(patch)
        if sample_hatch is not None:
            legend_elements = [
                Patch(facecolor=sample_colors[s], edgecolor="w", hatch=sample_hatch[s], label=s) for s in sample_order
            ]
        else:
            legend_elements = [Patch(facecolor=sample_colors[s], edgecolor="w", label=s) for s in sample_order]
        
        plt.legend(handles=legend_elements,  bbox_to_anchor=(1.05, 1), loc='upper left')
    
    if segment is not None:
        plot_segment_lines(segment, site_genomic_pos=site_genomic_pos)
    
    if highlights is not None:
        plot_vertical_highlights(
            highlights,
            highlight_color=highlight_color,
            site_genomic_pos=site_genomic_pos,
        )


def plot_vertical_highlights(highlights: List[Tuple], highlight_color="b", site_genomic_pos=None):
    if highlights is not None:
        for i, highlight_range in enumerate(highlights):
            if site_genomic_pos is not None:
                highlight_range = [site_genomic_pos[p] for p in highlight_range]
            
            if isinstance(highlight_color, list):
                color = highlight_color[i]
            else:
                color = highlight_color
            
            plt.axvspan(
                highlight_range[0],
                highlight_range[1],
                color=color,
                alpha=0.25,
            )


def plot_segment_lines(segment: np.array, site_genomic_pos=None):
    ymax = plt.gca().get_ylim()[1]
    for i in range(1, len(segment)):
        if segment[i] == segment[i - 1]:
            continue
        if site_genomic_pos is not None:
            x = (site_genomic_pos[i - 1] + site_genomic_pos[i]) / 2
        else:
            x = i - 1 + 0.5
        plt.plot((x, x), (0, ymax - 1), c="k")


def plot_met_profile_from_matrix(matrix: SparseMethylationMatrixContainer, **kwargs):
    plot_met_profile(
        matrix.met_matrix.todense(),
        samples=matrix.read_samples,
        site_genomic_pos=matrix.genomic_coord,
        site_genomic_pos_end=matrix.genomic_coord_end,
        **kwargs
    )
