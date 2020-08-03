from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import nanoepitools.math as nem
import nanoepitools.nanopolish_calls as npc

default_base_colors = {'C': 'r', 'T': 'g', 'G': 'b', 'A': 'y'}


def plot_met_multiple_hist(metcalls: List[pd.DataFrame], bins=200,
                           bound=20, alpha=0.5, colors=None, labels=None,
                           normalize_histograms=False,
                           title='Methylation log likelihood ratios'):
    """
    Plots a histogram of log-likelihood ratios from nanopolish output.
    Plots multiple histograms over each other, to compare different
    samples
    :param metcalls: a list of dataframes as produced by nanopolish
    :param bins: The number of bins for the histogram (default:200)
    :param bound: Clip llrs with an absolute value larger than this (
    default:20)
    :param alpha: Opacity of the histogram (default:0.5)
    :param colors: Color for each group
    :param labels: Label for each group (for legend)
    :param normalize_histograms: whether to rescale histograms so they have
    the same total
    :param title: Title for the plot
    :return:
    """
    bins = np.arange(-bound, bound+1, (2 * (bound + 1)) / bins)
    for i, metcall in enumerate(metcalls):
        llr_type = metcall['log_lik_ratio']
        llr_type = np.clip(llr_type, -bound, bound)
        weights = np.ones(len(llr_type))
        if normalize_histograms:
            weights = weights / len(llr_type)
        plt.hist(llr_type, weights=weights, bins=bins, alpha=alpha,
                 color=(colors[i] if colors is not None else None),
                 label=(labels[i] if labels is not None else None))
    plt.xlabel("Methylation log-likelihood ratio")
    plt.ylabel('Frequency')
    plt.title(title)
    plt.legend()
    plt.tight_layout()


def plot_met_hist_types(metcall: pd.DataFrame, typecol='type', typedict=None,
                        typecolor=None, **kwargs):
    """
    Plots a histogram of log-likelihood ratios from nanopolish output.
    Plots multiple histograms over each other, to compare different
    samples. kwargs will be passed on to plot_met_multiple_hist
    :param metcall: the dataframe as produced by nanopolish
    :param typecol: which column distinguishes the types/samples (
    default='type')
    :param typedict: Provide a readable type/sample name for the legend
    :param typecolor: Provide a color for each type/sample
    :return:
    """
    types = list(set(metcall[typecol]))
    if typedict is None:
        typedict = {t: t for t in types}
    if typecolor is None:
        typecolor = {t: None for t in types}

    metcalls = [metcall.loc[metcall[typecol] == t] for t in types]
    labels = [typedict[t] for t in types]
    colors = [typecolor[t] for t in types]
    plot_met_multiple_hist(metcalls, labels=labels, colors=colors, **kwargs)


def plot_read_bs_dist(metcall: pd.DataFrame, llr_threshold=2.5, min_calls=20):
    """
    Binarizes the methylation calls and then computes a beta-score
    (methylation rate) per read. Then plots a histogram of the distribution
    of read methylation rates.
    :param metcall: the dataframe as produced by nanopolish
    :param llr_threshold: exclude calls that are closer than this threshold
    to zero
    :param min_calls: exclude reads with fewer (included) calls
    """
    metcall = metcall[['read_name', 'log_lik_ratio']].copy()
    metcall['ismet'] = metcall['log_lik_ratio'] > llr_threshold
    metcall['isunmet'] = metcall['log_lik_ratio'] < -llr_threshold
    metcall = metcall.groupby('read_name').sum()
    metcall['bs'] = metcall['ismet'] / (metcall['ismet'] + metcall['isunmet'])
    metcall = metcall.loc[(metcall['ismet'] + metcall['isunmet']) > min_calls]
    plt.hist(metcall['bs'])


def plot_kmer_lenght_vs_uncertainty(metcall: pd.DataFrame,
                                    has_correct_col=False,
                                    uncertainty_method='linear'):
    """
    Plots 3 (or 4) plots as subplots:
     * Distribution of subsequence length
     * Call uncertainty dependent on sequence length
     * Cumulative uncertainty over sequence length
     * (if 'has_correct_col' is True) Prediction errors on equence length
    :param metcall: the dataframe as produced by nanopolish
    :param has_correct_col: whether or not an additional "correct" column
    has been added, that is needed for the last plot and tells whether a
    call is correct
    :param uncertainty_method: the uncertainty method. See docstring of
    llr_to_uncertainty. Default: 'linear'
    :return: figure object
    """
    uncertainty = nem.llr_to_uncertainty(metcall['log_lik_ratio'],
                                         method=uncertainty_method)
    kmerlen = metcall['sequence'].map(lambda x: len(x))

    uncertainty_len_df = pd.DataFrame({'kmerlen': kmerlen,
                                       'uncertainty': uncertainty})

    len_uncertainty_np_cum = np.cumsum(uncertainty_len_df.groupby(
        'kmerlen').sum())
    len_uncertainty_np = uncertainty_len_df.groupby('kmerlen').mean()
    len_incidents_np = uncertainty_len_df.groupby('kmerlen').count()

    fig = plt.figure(figsize=(15, 10))
    fig.patch.set_facecolor('w')
    axes = fig.subplots(2, 2)
    axes[0, 0].scatter(len_uncertainty_np.index, np.log10(
        len_incidents_np.uncertainty))
    axes[0, 0].set_title('Sequence length frequency (log transformed)')
    axes[0, 0].set_xlabel("Subsequence length")
    axes[0, 0].set_ylabel("Frequency log10")

    axes[0, 1].scatter(len_uncertainty_np.index, len_uncertainty_np.uncertainty)
    axes[0, 1].set_title('Mean uncertainty per sequence length')
    axes[0, 1].set_xlabel("Subsequence length")
    axes[0, 1].set_ylabel("Uncertainty")

    axes[1, 0].scatter(len_uncertainty_np_cum.index,
                       len_uncertainty_np_cum.uncertainty)
    axes[1, 0].set_title('Cumulative uncertainty over sequence length')
    axes[1, 0].set_xlabel("Subsequence length")
    axes[1, 0].set_ylabel("Cumulative Uncertainty")

    if has_correct_col:
        len_correct_df = pd.DataFrame({'kmerlen': kmerlen,
                                       'correct': metcall['correct']})
        len_correct = len_correct_df.groupby('kmerlen').sum()
        axes[1, 1].scatter(len_correct.index,
                           1 - (
                                   len_correct.correct / len_incidents_np.uncertainty))
        axes[1, 1].set_title('Error rate per sequence length')
        axes[1, 1].set_xlabel("Subsequence length")
        axes[1, 1].set_ylabel("Error rate")

    return fig


def plot_kmer_function(kmer_x_axis: pd.Series,
                       kmer_y_axis: pd.Series, base_colors,
                       x_axis_label, y_axis_label):
    # Tuples represent x axis, y axis, and charcter in sequence
    index_map = [(0, 0, 0), (0, 1, 1), (1, 0, 4), (1, 1, 5)]

    fig = plt.figure(figsize=(15, 10))
    fig.patch.set_facecolor('w')
    axes = fig.subplots(2, 2)
    for index_tuple in index_map:
        ax = axes[index_tuple[0], index_tuple[1]]

        for base in base_colors.keys():
            index = kmer_x_axis.index.map(lambda x: x[index_tuple[2]] == base)
            ax.scatter(kmer_x_axis[index], kmer_y_axis[index],
                       c=base_colors[base], label=base)
        ax.legend()
        ax.set_title('K-mer per base %d' % (index_tuple[2] + 1))
        ax.set_xlabel(x_axis_label)
        ax.set_ylabel(y_axis_label)

    return fig


def plot_kmer_uncertainty(metcall: pd.DataFrame,
                          base_colors=default_base_colors,
                          uncertainty_method='linear'):
    """
    Plot uncertainty of prediction per kmer
    :param metcall: the dataframe as produced by nanopolish
    :param base_colors: Dictionary that assigns color to each base
    :param uncertainty_method: the uncertainty method. See docstring of
    llr_to_uncertainty. Default: 'linear'
    :return: figure object
    """
    kmer_uncertainty = \
        npc.compute_kmer_uncertainty(metcall,
                                     uncertainty_method=uncertainty_method)
    kmer_incidents = npc.count_kmer_incidents(metcall)

    plot_kmer_function(kmer_incidents, kmer_uncertainty, base_colors,
                       'k-mer frequency', 'Mean Uncertainty')


def plot_kmer_error_rate(metcall: pd.DataFrame,
                         base_colors=default_base_colors,
                         error_method='llr'):
    """
    Plot error rate of prediction per kmer
    :param metcall: the dataframe as produced by nanopolish. Must have
    additional 'correct' column
    :param base_colors: Dictionary that assigns color to each base
    :param error_method: How to quantify error. Default: 'llr'
    :return: figure object
    """

    kmer_error = npc.compute_kmer_error(metcall, error_method=error_method)
    kmer_incidents = npc.count_kmer_incidents(metcall)

    plot_kmer_function(kmer_incidents, kmer_error, base_colors,
                       'k-mer frequency', 'Mean log likelihood ratio error')


def plot_kmer_error_vs_uncertainty(metcall: pd.DataFrame,
                                   base_colors=default_base_colors,
                                   uncertainty_method='linear',
                                   error_method='llr'):
    """
    Plot error rate of prediction per kmer
    :param metcall: the dataframe as produced by nanopolish. Must have
    additional 'correct' column
    :param base_colors: Dictionary that assigns color to each base
    :param uncertainty_method: the uncertainty method. See docstring of
    llr_to_uncertainty. Default: 'linear'
    :param error_method: How to quantify error. Default: 'llr'
    :return: figure object
    """

    kmer_error = npc.compute_kmer_error(metcall, error_method=error_method)
    kmer_uncertainty = \
        npc.compute_kmer_uncertainty(metcall,
                                     uncertainty_method=uncertainty_method)

    plot_kmer_function(kmer_uncertainty, kmer_error, base_colors,
                       'Mean Uncertainty', 'Mean error rate')
