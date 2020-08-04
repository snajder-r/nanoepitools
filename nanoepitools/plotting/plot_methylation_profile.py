import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List


def plot_met_profile(matrix: np.array, samples: np.array, sample_order:
                     List[str], sample_colors: Dict, site_genomic_pos=None,
                     segment: np.array = None):

    def val_to_color(val): return 1 - np.exp(-np.abs(val) * 0.5)

    y_off = 0
    start = 0
    end = matrix.shape[1]

    for s in sample_order:
        x = np.arange(start, end)
        part_matrix = matrix[:, x][(samples == s)]
        if site_genomic_pos is not None:
            x = site_genomic_pos[x]  # Translate to actual pos on chrom
        active_reads = np.array((part_matrix != 0).sum(axis=1)).flatten() > 0

        part_matrix = part_matrix[active_reads]
        hasval = np.array(part_matrix != 0).flatten()
        y = np.arange(part_matrix.shape[0]) + y_off

        x, y = np.meshgrid(x, y)
        x = x.flatten()[hasval]
        y = y.flatten()[hasval]
        matrix_data = np.array(part_matrix).flatten()[hasval]

        color = [[0, 1, 0, val_to_color(-v)] if v < 0 else [1, 0, 0,
                                                            val_to_color(v)]
                 for v in matrix_data]

        plt.scatter(x, y, c=color, marker='|', s=15)

        x = np.ones(part_matrix.shape[0]) * (x.max() + 20)
        y = np.arange(part_matrix.shape[0]) + y_off
        plt.scatter(x, y, c=sample_colors[s])

        y_off += part_matrix.shape[0]

    if segment is not None:
        plot_segment_lines(segment, site_genomic_pos=site_genomic_pos)


def plot_segment_lines(segment: np.array, site_genomic_pos=None):
    ymax = plt.gca().get_ylim()[1]
    for i in range(1, len(segment)):
        if segment[i] == segment[i - 1]:
            continue
        if site_genomic_pos is not None:
            i = site_genomic_pos[i]
        plt.plot((i - 1 + 0.5, i - 1 + 0.5), (0, ymax - 1), c='k')