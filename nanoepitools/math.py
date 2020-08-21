import numpy as np
from scipy.stats import rankdata


def llr_to_p(llr, prior=0.5):
    """
    Convert log-likelihood ratios log(p(x|a)/p(x|~a)) to posterior
    probabilty p(a|x) given a prior p(a). For unbiased prediction,
    leave prior at 0.5
    """
    return 1 / (1 + np.exp(-llr) * (1 / prior - 1))


def p_to_llr(p, prior=0.5):
    """
    Converts the posterior probability p(a|x) into a log-likelihood ratio
    log(p(x|a)/p(x|~a)) given a prior pa(a)
    """
    return -np.log(prior * (1 - p) / (p * (1 - prior)))


def llr_to_uncertainty(llr, method='linear'):
    if method == 'linear':
        p = llr_to_p(llr)
        return 0.5 - np.abs(0.5 - p)


def fdr_from_pvals(p_vals: np.ndarray) -> np.ndarray:
    """
    Computes FDR from p-values using the Benjamini-Hochberg method.
    :param p_vals: numpy array of p-values
    :return: numpy array of adjusted p-values
    """
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr
