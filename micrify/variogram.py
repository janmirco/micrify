import logging

import numba
import numpy as np
from numpy.typing import NDArray


def semivariogram(
    data: NDArray[np.float64],
    num_bins: int,
    max_lag: float,
    bin_edges: NDArray[np.float64],
    estimation_approach: str = "numba",
    hx: float = 1.0,
    hy: float = 1.0,
) -> tuple[list[np.float64], list[np.float64]]:
    if estimation_approach == "numba":
        sqdiff_sums, bin_counts = estimate_numba(data, num_bins=num_bins, max_lag=max_lag, bin_edges=bin_edges, hx=hx, hy=hy)
    elif estimation_approach == "vectorized":
        sqdiff_sums, bin_counts = estimate_vectorized(data, num_bins=num_bins, max_lag=max_lag, bin_edges=bin_edges, hx=hx, hy=hy)
    elif estimation_approach == "manual":
        sqdiff_sums, bin_counts = estimate_manual(data, num_bins=num_bins, max_lag=max_lag, bin_edges=bin_edges, hx=hx, hy=hy)
    else:
        raise NotImplementedError(f"Chosen estimation approach is not implemented. {estimation_approach = }. Choose one of ['numba', 'vectorized', 'manual'].")

    lags = []
    semivariogram_vals = []
    for bin_idx in range(num_bins):
        if bin_counts[bin_idx] > 0:
            lag_center = 0.5 * (bin_edges[bin_idx] + bin_edges[bin_idx + 1])
            lags.append(lag_center)
            semivariogram_vals.append(0.5 * sqdiff_sums[bin_idx] / bin_counts[bin_idx])
    return lags, semivariogram_vals


def estimate_manual(
    data: NDArray[np.float64],
    num_bins: int,
    max_lag: float,
    bin_edges: NDArray[np.float64],
    hx: float = 1.0,
    hy: float = 1.0,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    sqdiff_sums = np.zeros(num_bins)
    bin_counts = np.zeros(num_bins)
    for x1 in range(data.shape[0]):
        logging.info(f"[Manual] Start round {x1 + 1} out of {data.shape[0]}...")
        for y1 in range(data.shape[1]):
            for x2 in range(data.shape[0]):
                for y2 in range(data.shape[1]):
                    if (x1 < x2) or ((x1 == x2) and (y1 < y2)):  # unique pairs, no repeats
                        lag = np.sqrt(hx * (x2 - x1) ** 2.0 + hy * (y2 - y1) ** 2.0)
                        if (lag < 1.0e-08) or (lag > max_lag):
                            continue
                        bin_idx = np.searchsorted(bin_edges, lag, side="right") - 1
                        if (bin_idx < 0) or (bin_idx >= num_bins):
                            continue
                        sqdiff_sums[bin_idx] += (data[x1, y1] - data[x2, y2]) ** 2.0
                        bin_counts[bin_idx] += 1
    return sqdiff_sums, bin_counts


@numba.njit
def estimate_numba(
    data: NDArray[np.float64],
    num_bins: int,
    max_lag: float,
    bin_edges: NDArray[np.float64],
    hx: float = 1.0,
    hy: float = 1.0,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    sqdiff_sums = np.zeros(num_bins)
    bin_counts = np.zeros(num_bins)
    for x1 in range(data.shape[0]):
        for y1 in range(data.shape[1]):
            for x2 in range(data.shape[0]):
                for y2 in range(data.shape[1]):
                    if (x1 < x2) or ((x1 == x2) and (y1 < y2)):  # unique pairs, no repeats
                        lag = np.sqrt(hx * (x2 - x1) ** 2.0 + hy * (y2 - y1) ** 2.0)
                        if (lag < 1.0e-08) or (lag > max_lag):
                            continue
                        bin_idx = np.searchsorted(bin_edges, lag, side="right") - 1
                        if (bin_idx < 0) or (bin_idx >= num_bins):
                            continue
                        sqdiff_sums[bin_idx] += (data[x1, y1] - data[x2, y2]) ** 2.0
                        bin_counts[bin_idx] += 1
    return sqdiff_sums, bin_counts


def estimate_vectorized(
    data: NDArray[np.float64],
    num_bins: int,
    max_lag: float,
    bin_edges: NDArray[np.float64],
    hx: float = 1.0,
    hy: float = 1.0,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    num_data_points = data.shape[0] * data.shape[1]
    if num_data_points > 2_500:
        logging.info(f"Warning: Large grid size will fill up your memory! {num_data_points = }. Watch your memory and scale up carefully.")
    if num_data_points > 16_000:
        raise RuntimeError(f"Too large of a grid size for 32 GB memory! {num_data_points = }. Choose less than 16_000!")

    # Flatten data and create coordinate arrays
    coords = np.array([(x, y) for x in range(data.shape[0]) for y in range(data.shape[1])])
    vals = data.ravel()

    # Compute all pairwise distances and squared differences
    diff = vals[:, None] - vals[None, :]
    diff2 = diff**2

    # Pairwise distances
    dcoords = coords[:, None, :] - coords[None, :, :]
    dists = np.sqrt(hx * (dcoords[:, :, 0] ** 2.0) + hy * (dcoords[:, :, 1] ** 2.0))

    # Select unique pairs (upper triangle, excluding diagonal)
    iu = np.triu_indices(len(vals), k=1)

    lags_all = dists[iu]
    diff2_all = diff2[iu]

    # Filter by max_lag
    mask = (lags_all > 0) & (lags_all <= max_lag)
    lags_all = lags_all[mask]
    diff2_all = diff2_all[mask]

    # Bin and calculate semivariogram
    bin_indices = np.searchsorted(bin_edges, lags_all, side="right") - 1
    sqdiff_sums = np.zeros(num_bins)
    bin_counts = np.zeros(num_bins)
    for idx, d2 in zip(bin_indices, diff2_all):
        if 0 <= idx < num_bins:
            sqdiff_sums[idx] += d2
            bin_counts[idx] += 1

    return sqdiff_sums, bin_counts
