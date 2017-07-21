# -*- coding: utf-8 -*-

# Copied from https://github.com/KarchinLab/probabilistic2020/blob/master/prob2020/python/mymath.py

import warnings
import numpy as np
from scipy.misc import logsumexp

import ksiga.sparse_util as su


def shannon_entropy(p):
    """Calculates shannon entropy in bits.

    Parameters
    ----------
    p : np.array
        array of probabilities

    Returns
    -------
    shannon entropy in bits
    """
    return -np.sum(np.where(p!=0, p * np.log2(p), 0))


def log_shannon_entropy(log_p):
    """Calculates shannon entropy in bits.

    Parameters
    ----------
    p : np.array
        array of probabilities

    Returns
    -------
    shannon entropy in bits
    """
    out = -logsumexp(log_p + np.log(log_p))
    return out


def max_shannon_entropy(n):
    """Returns max possible entropy given "n" mutations.

    The maximum possible entropy is the entropy of the
    uniform distribution. The uniform distribution has
    entropy equal to log(n) (I will use base 2).

    Parameters
    ----------
    n : int
        total mutation counts

    Returns
    -------
    max possible shannon entropy in bits
    """
    if n <= 0:
        return 0.
    return float(np.log2(n))


def normalized_mutation_entropy(counts, total_cts=None):
    """Calculate the normalized mutation entropy based on a list/array
    of mutation counts.

    Note: Any grouping of mutation counts together should be done before hand

    Parameters
    ----------
    counts : np.array_like
        array/list of mutation counts

    Returns
    -------
    norm_ent : float
        normalized entropy of mutation count distribution.
    """
    cts = np.asarray(counts, dtype=float)
    if total_cts is None:
        total_cts = np.sum(cts)
    if total_cts > 1:
        p = cts / total_cts
        ent = shannon_entropy(p)
        max_ent = max_shannon_entropy(total_cts)
        norm_ent = ent / max_ent
    else:
        norm_ent = 1.0
    return norm_ent


def kl_divergence(p, q):
    """Compute the Kullback-Leibler (KL) divergence for discrete distributions.

    Parameters
    ----------
    p : np.array
        "Ideal"/"true" Probability distribution
    q : np.array
        Approximation of probability distribution p

    Returns
    -------
    kl : float
        KL divergence of approximating p with the distribution q
    """
    # make sure numpy arrays are floats
    # p = p.astype(float)
    # q = q.astype(float)

    # compute kl divergence
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        kl = np.sum(np.where(p != 0, p*np.log2(p/q), 0))
        # warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')
    return kl


def js_divergence(p, q):
    """Compute the Jensen-Shannon Divergence between two discrete distributions.

    Parameters
    ----------
    p : np.array
        probability mass array (sums to 1)
    q : np.array
        probability mass array (sums to 1)

    Returns
    -------
    js_div : float
        js divergence between the two distrubtions
    """
    m = .5 * (p+q)
    js_div = .5*kl_divergence(p, m) + .5*kl_divergence(q, m)
    return js_div


def js_distance(p, q):
    """Compute the Jensen-Shannon distance between two discrete distributions.

    NOTE: JS divergence is not a metric but the sqrt of JS divergence is a
    metric and is called the JS distance.

    Parameters
    ----------
    p : np.array
        probability mass array (sums to 1)
    q : np.array
        probability mass array (sums to 1)

    Returns
    -------
    js_dist : float
        Jensen-Shannon distance between two discrete distributions
    """
    js_dist = np.sqrt(js_divergence(p, q))
    return js_dist


def sparse_kl_divergence(p, q):
    """Compute the Kullback-Leibler (KL) divergence on sparse matrix.

    Parameters
    ----------
    p : scipy csr_matrix (with 1 row)
        "Ideal"/"true" Probability distribution
    q : scipy csr_matrix (with 1 row)
        Approximation of probability distribution p

    Returns
    -------
    kl : float
        KL divergence of approximating p with the distribution q
    """
    # Get index that appear in both sparse matrix (assume that both are sorted)
    p_idx = su.searchsorted(p.indices, q.indices)
    q_idx = su.searchsorted(q.indices, p.indices)

    # Get only part where BOTH has value.
    p_val = p.data[p_idx]
    q_val = q.data[q_idx]

    kl = np.sum(p_val * np.log2(p_val/q_val))

    return kl


def sparse_js_divergence(p, q):
    """Compute the Jensen-Shannon Divergence between two discrete distributions.

    Parameters
    ----------
    p : np.array
        probability mass array (sums to 1)
    q : np.array
        probability mass array (sums to 1)

    Returns
    -------
    js_div : float
        js divergence between the two distrubtions
    """
    m = .5 * (p+q)
    js_div = .5*sparse_kl_divergence(p, m) + .5*sparse_kl_divergence(q, m)
    return js_div


def sparse_js_distance(p, q):
    """Compute the Jensen-Shannon distance between two discrete distributions.

    NOTE: JS divergence is not a metric but the sqrt of JS divergence is a
    metric and is called the JS distance.

    Parameters
    ----------
    p : np.array
        probability mass array (sums to 1)
    q : np.array
        probability mass array (sums to 1)

    Returns
    -------
    js_dist : float
        Jensen-Shannon distance between two discrete distributions
    """
    js_dist = np.sqrt(sparse_js_divergence(p, q))
    return js_dist
