import re
from pathlib import Path
from typing import List, Union, Dict

import numpy as np
import pandas as pd

import nanoepitools.math as nem


def load_nanopore_metcalls_from_tsv(input_folder: Union[str, Path],
                                    samples: List[str],
                                    mettypes=['cpg', 'gpc', 'dam']):
    # matches filenames like "footprinting_24_met_cpg.tsv" or
    # "invitro12_2_0_met_dam.tsv"
    filename_regex = re.compile('^(.*)_([0-9]*)_met_([^_]*)\\.tsv$')
    input_folder = Path(input_folder)

    all_sample_met = dict()
    for sample in samples:
        all_sample_met[sample] = dict()
        for mettype in mettypes:
            all_sample_met[sample][mettype] = None

    for f in input_folder.iterdir():
        mitch = filename_regex.match(f.name)
        if mitch is not None:
            sample, batch, mettype = mitch.groups()
            # Only load select sample and mettypes
            if sample not in samples or mettype not in mettypes:
                continue

            met_part = pd.read_csv(f, sep='\t')
            if all_sample_met[sample][mettype] is None:
                all_sample_met[sample][mettype] = met_part
            else:
                all_sample_met[sample][mettype] = all_sample_met[sample][
                    mettype].append(met_part)

    return all_sample_met


def load_merged_nanopore_metcalls(input_folder: Union[str, Path],
                                  samples: List[str], chroms: List[str],
                                  mettypes: List[str] = ['cpg', 'gpc',
                                                         'dam']) -> \
        Dict[str, Dict[str, Dict[str, pd.DataFrame]]]:
    """
    Loads pickled nanopolish methylation calls as pandas dataframes and
    organizes them by sample, chromosome, and methylation type.

    The output is a cascading tree of dictionaries:
        samplename -> (chromosome -> (methylation type -> dataframe))

    :param input_folder: folder containing subfolders for each sample
    :param samples: which samples to include
    :param chroms: which chromosomes to include
    :param mettypes: which methylation types to include (default: cpg,
    gpc and dam)
    :return: tree of dictionaries with dataframes as leafes
    """
    base_filename = '{chrom}_met_{mettype}.pkl'
    input_folder = Path(input_folder)

    all_sample_met = dict()
    for sample in samples:
        all_sample_met[sample] = dict()

        sample_dir = input_folder.joinpath(sample)

        for chrom in chroms:
            all_sample_met[sample][chrom] = dict()
            for mettype in mettypes:
                all_sample_met[sample][chrom][mettype] = None

                filename = base_filename.format(chrom=chrom, mettype=mettype)
                filepath = sample_dir.joinpath(filename)

                all_sample_met[sample][chrom][mettype] = pd.read_pickle(
                    filepath, compression='gzip')

    return all_sample_met


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
    assert (all(metcall_onecpg['sequence'].map(lambda x: len(x) == 11)))
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
    metcall_onecpg['uncertainty'] = nem.llr_to_uncertainty(
        metcall_onecpg['log_lik_ratio'], method=uncertainty_method)

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
