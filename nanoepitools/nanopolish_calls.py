import pandas as pd
import numpy as np
import nanoepitools.math as nem


def get_only_single_cpg_calls(metcall: pd.DataFrame):
    """
    Specifically for methylation calls for CpGs. Since Nanopolish groups calls
    that are close to each other, this function filters out only single-cpg
    calls (for things like k-mer uncertainty analysis)
    :param metcall: the dataframe as produced by nanopolis
    :return:
    """
    return metcall.loc[metcall['sequence'].map(lambda x: len(x) == 11)]


def add_sixmer_column(metcall_onecpg: pd.DataFrame):
    """
    Specifically for methylation calls for CpGs. Adds a column with the
    sixmer around the CpG site. Assumes only single cpg calls
    :param metcall_onecpg: the dataframe as produced by nanopolis
    """
    assert(all(metcall_onecpg['sequence'].map(lambda x: len(x) == 11)))
    metcall_onecpg['sixmer'] = metcall_onecpg['sequence'].map(lambda x: x[3:-2])


def compute_kmer_uncertainty(metcall: pd.DataFrame,
                             uncertainty_method='linear'):
    """
    :param metcall: the dataframe as produced by nanopolis
    :param uncertainty_method: the uncertainty method. See docstring of
    llr_to_uncertainty. Default: 'linear'
    :return: series with kmers as index and uncertainty as values
    """
    metcall_onecpg = get_only_single_cpg_calls(metcall).copy()
    add_sixmer_column(metcall_onecpg)
    metcall_onecpg['uncertainty'] = \
        nem.llr_to_uncertainty(metcall_onecpg['log_lik_ratio'],
                               method=uncertainty_method)

    metcall_onecpg = metcall_onecpg[['sixmer', 'uncertainty']]
    kmer_uncertainty = metcall_onecpg.groupby('sixmer').mean()
    return kmer_uncertainty['uncertainty']


def count_kmer_incidents(metcall: pd.DataFrame):
    """
    :param metcall: the dataframe as produced by nanopolis
    :return: series with kmers as index and counts as values
    """
    metcall_onecpg = get_only_single_cpg_calls(metcall).copy()
    add_sixmer_column(metcall_onecpg)
    metcall_onecpg = metcall_onecpg[['sixmer', 'log_lik_ratio']]
    kmer_incidents = metcall_onecpg.groupby('sixmer').count()
    return kmer_incidents['log_lik_ratio']


def compute_kmer_error(metcall: pd.DataFrame, error_method='llr'):
    """
    :param metcall: the dataframe as produced by nanopolis
    :param error_method: what method to use to compute the error.
    Default: 'llr'
    :return: series with kmers as index and errors as values
    """
    metcall_onecpg = get_only_single_cpg_calls(metcall).copy()
    add_sixmer_column(metcall_onecpg)

    if error_method == 'llr':
        p = np.abs(metcall_onecpg['log_lik_ratio'])
        p[metcall_onecpg['correct']] = -p
        metcall_onecpg['error'] = p
    elif error_method == 'confusion_rate':
        metcall_onecpg['error'] = ~metcall_onecpg['correct']
    else:
        raise ValueError('Invalid error method')
    kmer_error = metcall_onecpg.groupby('sixmer').mean()
    return kmer_error['error']