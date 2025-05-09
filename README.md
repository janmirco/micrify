# 🔬 Micrify

Micrify is a toolbox developed for my research in computational mechanics and multiscale analysis of materials.
By applying concepts from [spatial statistics](https://en.wikipedia.org/wiki/Spatial_statistics), I analyze spatial correlations within materials to better understand statistical fluctuations arising from, for example, production variability or microstructural heterogeneity.
Combined with the [finite element method](https://en.wikipedia.org/wiki/Finite_element_method), this toolbox enables the simulation and analysis of how microscopic spatial randomness influences the macroscopic behavior and properties of materials.

## Usage

How to use Micrify can be found in the [examples](examples) directory.
For example, a comparison of a variogram estimation implementation with the implementation of [GSTools](https://geostat-framework.readthedocs.io/projects/gstools) can be found [here](examples/variogram-comparison).

## Installation

### Directly via GitHub

Go to your `uv` project and install `micrify`.

```bash
uv add git+https://github.com/janmirco/micrify.git
```

For future updates, simply run the following code block.

```bash
uv lock --upgrade
uv sync
```

### Locally in editable mode

1. Clone repository.

```bash
git clone https://github.com/janmirco/micrify.git
```

2. Go to your `uv` project and install `micrify` in editable mode.

```bash
uv add --editable <PATH_TO_HELLO_REPO>
```

## Software used

- Language: [Python](https://www.python.org/)
- Code execution: [GNU Make](https://www.gnu.org/software/make/)
- Package and project manager: [uv](https://docs.astral.sh/uv/)
- Python packages:
    - [Utly](https://github.com/janmirco/utly)
    - [GSTools](https://geostat-framework.readthedocs.io/projects/gstools)
    - [Matplotlib](https://matplotlib.org/)
    - [Numba](https://numba.pydata.org/)
    - [NumPy](https://numpy.org/)
