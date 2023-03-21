import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from tqdm import tqdm
from loguru import logger

from pyunicorn.timeseries import RecurrencePlot
from pyunicorn.timeseries.surrogates import Surrogates

np.seterr(all="ignore")


def pearson(yi, yj):
    return sp.stats.pearsonr(yi, yj)[0]


def rmd(yi, yj, recurrence_rate=0.1, **kwargs):
    """
    Recurrence based measure of dependence between to time series `yi` and `yj`.

    Parameters
    ----------
    yi: (n,) array
        First time series.
    yj: (n,) array
        Second time series.
    recurrence_rate: float
        Recurrence rate of the recurrence plots that should be used.
    kwargs: dict
        Additional arguments to `pyunicorn.RecurrenceRate`

    Returns
    -------
    rmd: float
        The recurrence based measure of dependence between the two time series.
        Returns zero when there is no dependence.

    """

    # generate the recurrence plots
    # TODO: the function uses the supremum norm, is this okay?
    irp = RecurrencePlot(yi, recurrence_rate=recurrence_rate, silence_level=10, **kwargs)
    jrp = RecurrencePlot(yj, recurrence_rate=recurrence_rate, silence_level=10, **kwargs)

    # get the recurrence matrices
    irp = irp.recurrence_matrix().copy()
    jrp = jrp.recurrence_matrix().copy()

    # calculate the column recurrence rate
    Pij = np.mean(irp*jrp, axis=1)
    Pi  = np.mean(irp, axis=1)
    Pj  = np.mean(jrp, axis=1)

    # calculate the recurrence based measure of dependence
    rmd = np.log2(np.sum(Pij/(Pi*Pj)))

    return rmd if np.isfinite(rmd) else 0


def create_surrogates(y, n_surrogates=1, threshold=None, dimension=1, delay=0):
    """
    Generates twin surrogates of a time series `y`.

    The surrogates are generated using the function `pyunicorn.twin_surrogates`.

    Parameters
    ----------
    y: (n,) array
        Time series for which the surrogates should be generated.
    n_surrogates: int
        Number of surrogates to generate.
    threshold: float
        Distance threshold in phase space to determine twins.
    dimension: int, default 0
        Dimension for embedding
    delay: int, default 0
        Delay for twin creation.

    Returns
    -------
    surrogates: (n_surrogates, n) or (n)
        Sequence of surrogates or only generated surrogate.
    """

    # use 1% of the phase space extension as a threshold for the recurrence plot
    if threshold is None:
        threshold = 0.01*np.max(np.diff(y))

    # set up surrogate generator
    surrogate_generator = Surrogates(y, silence_level=100)
    surrogate_generator.clear_cache()

    # generate S surrogates
    surrogates = np.array([surrogate_generator.twin_surrogates(np.array([y]), dimension, delay, threshold) for s in range(n_surrogates)]).squeeze()

    return surrogates if n_surrogates > 1 else np.array([surrogates])


def correlation_matrix(ys, n_surrogates=10, recurrence_rate=0.1, metric=pearson,
                       qlower=0.025, qupper=0.975, **kwargs):

    # how many individual time series
    N = len(ys)

    # empty matrix for correlation metric
    C = np.zeros((N, N))
    # empty matrix for significance
    S = np.zeros((N, N))

    logger.info(f"Network with {N} nodes can have {N*(N+1)/2 - N} possible links")
    logger.info(f"Using significance interval [{qlower:.3f}, {qupper:.3f}) based on {n_surrogates} surrogates")
    logger.warning(f"Calculating correlation matrix, get a cup of coffee...")

    N = len(ys)
    with tqdm(total=N*(N+1)/2 - N) as pbar:

        for i ,yi in enumerate(ys):
            for j, yj in enumerate(ys):

                if i > j:

                    # calculate the metric for this pair of yi, yj
                    C[i, j] = metric(yi, yj, **kwargs)

                    # retrieve the recurrence threshold
                    rp = RecurrencePlot(yi, recurrence_rate=recurrence_rate, silence_level=10)
                    threshold = rp.threshold_from_recurrence_rate(rp.distance_matrix(1, "supremum"), recurrence_rate)
                    # TODO: the function uses the supremum norm, is this okay?

                    # create S surrogates of the series yj
                    surrogates = create_surrogates(yj, n_surrogates, threshold=threshold)

                    # for S surrogates, calculate the metric
                    m = [metric(yi, s) for s in surrogates]
                    q = np.quantile(m, [qlower, qupper])

                    # map whether effect is statistically significant
                    S[i, j] = (C[i, j] < q[0]) ^ (q[1] < C[i, j])
                    S[i, j] *= abs(np.mean(m) - C[i, j])

                    pbar.update()

    return C + C.T, S + S.T


def adjacency_matrix(C, S, p = 0.05, binary=True):
    S = S > np.quantile(S, 1 - p)
    return C*S != 0 if binary else C*S


def link_distances(adjacency_matrix, distance_matrix):
    """

    Returns a flattened list of all link lengths in a network.

    For a network represented as a adjacency matrix, link lengths are calculated based on an distance matrix.

    Parameters
    ----------
    adjacency_matrix: (n, n) matrix
        The adjacency matrix where `True` entries represent connected nodes.
    distance_matrix: (n, n) matrix
        The distance matrix where entries `Dij` represent the distance of nodes `i` and `j`.

    Returns
    -------
    link lengths: array
        A flattened list of all link lengths in the network.

    """

    # rename the
    A = adjacency_matrix
    D = distance_matrix

    # take the upper triangle without the diagonal
    A = A[np.triu_indices(len(A), k=1)]
    D = D[np.triu_indices(len(D), k=1)]

    # return a flattened version of the link distances
    return D[A == True].flatten()


def create_surrogate_adjacency_matrix(adjacency_matrix, distance_matrix, n_surrogates=1, plot=None, **kwargs):
    """
    Generates surrogate adjacency matrices from a given adjacency matrix.
    Surrogates share both network link density and link distance distribution with the original network.

    Based on a given distance matrix the algorithm derives a distribution of link lengths in the original network.
    Equally, the distribution that would be possible if all links in the network are connected is calculated in the same
    distance interval. This interval is limited to the minimal and maximal length in the original network.

    The surrogate network is constructed by treating the ratio of seen link length count vs. possible link length count
    as a probability of connecting the links. This preserves both the overall link density and the distribution.

    Parameters
    ----------
    adjacency_matrix: (n, n) matrix
        The adjacency matrix of the network for which surrogates should be created.
    distance_matrix: (n, n) matrix
        Distance matrix that holds the distances between nodes in the adjacency matrix.
    n_surrogates: int
        Number of surrogate matrices that should be created. Quicker to set this here
        than running the whole function multiple times.
    plot: bool or str
        Whether to plot the link length distribution. When a boolean is given, the plot
        is shown. Strings are interpreted as file names under which the figure is saved.
    kwargs: dict
        Additional parameters passed to `np.histogram` to generate the
        link length distribution

    Returns
    -------
    surrogates: (n_surrogates, n, n) array or (n, n) array
        A sequence of surrogates or the surrogate when only one is generated.

    """

    # calculate the link distances in the adjacency matrix
    d = link_distances(adjacency_matrix, distance_matrix)

    logger.info(f"Calculated link lengths for {len(d)} links")

    # calculate probability that a link distance is seen
    counts_seen, bins = np.histogram(d, **kwargs)
    counts_base, _    = np.histogram(distance_matrix[np.triu_indices(distance_matrix.shape[0], k=1)], bins=bins)

    # maximum allowed distance range in the surrogate
    dmin = np.min(d)
    dmax = np.max(d)

    # center the bins
    bins = (bins[:-1] + bins[1:])/2

    logger.info(f"Calculated seen and possible distributions in the interval [{bins[0]:.2f}, {bins[-1]:.2f})")

    # calculate the link probability per distance
    probs = np.nan_to_num(counts_seen/counts_base, nan=0, posinf=0, neginf=0)

    # set up the n_surrogate matrices
    S = np.zeros((n_surrogates, *adjacency_matrix.shape), dtype=bool)

    logger.info("Generating surrogates, drink a cup of tea...")

    N = adjacency_matrix.shape[0]
    with tqdm(total=N*(N+1)/2 - N) as pbar:
        for i in range(N):
            for j in range(N):

                if i > j:
                    # get the bin index of the distance matrix at this location
                    bin_index = np.argmin(np.abs(bins - distance_matrix[i, j]))

                    # get the probability that for a link distance
                    prob = probs[bin_index]
                    prob *= (dmin < distance_matrix[i, j] < dmax)

                    # randomly connect the nodes based on the distance
                    S[:, i, j] = np.where(np.random.random(size=n_surrogates) < prob, True, False)

                    # update the pbar
                    pbar.update()

    # transpose all matrices
    S = S + np.transpose(S, axes=(0, 2, 1))

    # do the plots
    try:
        if plot is not None:
            logger.info("Preparing plots")
            plt.bar(bins, counts_base, width=bins[1] - bins[0], color="lightgray", label="possible")
            plt.bar(bins, counts_seen, width=bins[1] - bins[0], color="darkgray", label="seen")
            plt.title("Distance Distribution in Network")
            plt.legend()
            if type(plot) is str:
                plt.savefig(plot)
            else:
                plt.show()
    except Exception as e:
        logger.error(f"Error during plotting: {e}")

    return S if n_surrogates > 1 else S[0]