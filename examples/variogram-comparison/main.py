import logging
from pathlib import Path

import gstools as gs
import matplotlib.pyplot as plt
import micrify as mc
import numpy as np
import rippl as rp

# Set Parameters
GRID_SIZE = 60
TOTAL_NUM_DATA_POINTS = GRID_SIZE * GRID_SIZE
MAX_LAG = GRID_SIZE / 2.0
MEAN = 1.0
STD = 0.05
LEN_SCALE = 4.0
SPATIALLY_CORRELATED = True
GSTOOLS_USED = True
NUM_BINS = int(GRID_SIZE / 1.0)  # only when `GSTOOLS_USED = False`
ESTIMATION_APPROACH = "numba"  # ["numba", "vectorized", "manual"]


def main() -> None:
    # Set up output dir and logging
    output_dir = rp.path.set_up()
    rp.log.set_up(output_dir)

    # Log paths
    logging.info(f"{Path.cwd() = }")
    logging.info(f"{output_dir = }")

    # Log parameters
    logging.info(f"{GRID_SIZE = }")
    logging.info(f"{TOTAL_NUM_DATA_POINTS = }")
    logging.info(f"{MAX_LAG = }")
    logging.info(f"{MEAN = }")
    logging.info(f"{STD = }")
    logging.info(f"{LEN_SCALE = }")
    logging.info(f"{SPATIALLY_CORRELATED = }")
    logging.info(f"{GSTOOLS_USED = }")
    logging.info(f"{ESTIMATION_APPROACH = }")

    if not SPATIALLY_CORRELATED:
        # Without spatial correlation
        data = np.random.normal(loc=MEAN, scale=STD, size=(GRID_SIZE, GRID_SIZE))
    else:
        # With spatial correlation
        x = np.arange(GRID_SIZE)
        y = np.arange(GRID_SIZE)
        xx, yy = np.meshgrid(x, y, indexing="ij")
        model = gs.Gaussian(dim=2, var=STD**2.0, len_scale=LEN_SCALE)
        srf = gs.SRF(model, mean=MEAN)
        data = srf.structured([x, y])
        logging.info(f"{model = }")
        logging.info(f"{srf = }")

    # logging.info(data)
    logging.info(f"{data.min() = :.4f}")
    logging.info(f"{data.max() = :.4f}")
    logging.info(f"{data.mean() = :.4f}")
    logging.info(f"{data.std() = :.4f}")
    logging.info(f"{data.var() = :.4f}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    if GSTOOLS_USED:
        x = np.arange(GRID_SIZE)
        y = np.arange(GRID_SIZE)
        xx, yy = np.meshgrid(x, y)
        coords = np.vstack([xx.ravel(), yy.ravel()])  # shape (2, N)
        vals = data.ravel()  # shape (N,)
        logging.info("[GSTools] Start variogram estimation...")
        with mc.timer.Timer("GSTools") as t1:
            lags, semivariogram_vals = gs.vario_estimate(coords, vals, max_dist=MAX_LAG)
        t1.print_time(milli=False)
        logging.info("[GSTools] Done.")
        nonzero_mask = semivariogram_vals > 0
        lags = lags[nonzero_mask]
        semivariogram_vals = semivariogram_vals[nonzero_mask]
        axes[0].scatter(lags, semivariogram_vals, color="tab:green", marker="x", zorder=10, label="GSTools semivariogram")

        num_bins = len(lags)
    else:
        num_bins = NUM_BINS

    logging.info(f"{num_bins = }")
    bin_edges = np.linspace(0.0, MAX_LAG, num_bins + 1)
    logging.debug(f"{bin_edges = }")

    logging.info("[Micrify] Start variogram estimation...")
    with mc.timer.Timer("Micrify") as t2:
        lags, semivariogram_vals = mc.variogram.semivariogram(
            data,
            num_bins=num_bins,
            max_lag=MAX_LAG,
            bin_edges=bin_edges,
            estimation_approach=ESTIMATION_APPROACH,
        )
    t2.print_time(milli=False)
    logging.info("[Micrify] Done.")

    logging.info(f"Time comparison: GSTools / Micrify = {t1.time / t2.time:.4f}")

    axes[0].scatter(lags, semivariogram_vals, color="tab:blue", label="Micrify semivariogram")
    axes[0].plot(lags, data.var() * np.ones_like(lags), linestyle="dashed", color="tab:red", label="Sample variance (sill)")

    axes[0].set_xlim(0, MAX_LAG)
    axes[0].set_xlabel("Lag (distance)")
    axes[0].set_ylabel("Semivariance")
    axes[0].legend()
    axes[0].grid(linestyle="dotted")

    axes[1].set_title("Random Field")
    im = axes[1].imshow(data, origin="lower", cmap="viridis", aspect="auto")
    fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04)

    plt.savefig(output_dir / Path("variogram_comparison.png"))
    plt.savefig(output_dir / Path("variogram_comparison.svg"))
    # plt.show()


if __name__ == "__main__":
    main()
