# Variogram estimation comparison

In this project, equidistant grids of data points are used.
For this special purpose, together with [Numba](https://numba.pydata.org/), variogram estimation is even faster than that of [GSTools](https://geostat-framework.readthedocs.io/projects/gstools) when using large enough grids (over ~2500 data points in my tests).

## Examplary study

### Parameters

```python
TOTAL_NUM_DATA_POINTS = GRID_SIZE * GRID_SIZE
MAX_LAG = GRID_SIZE / 2.0
MEAN = 1.0
STD = 0.05
LEN_SCALE = 4.0
SPATIALLY_CORRELATED = True
GSTOOLS_USED = True
ESTIMATION_APPROACH = "numba"  # ["numba", "vectorized", "manual"]
```

### Time comparison

| Grid size (-) | Num data points (-) | GSTools (s) | Micrify (s) | Ratio (-) |
| ------------: | ------------------: | ----------: | ----------: | --------: |
|             5 |                  25 |    0.000196 |    0.441334 |    0.0004 |
|            10 |                 100 |    0.000619 |    0.438565 |    0.0014 |
|            20 |                 400 |    0.007060 |    0.436542 |    0.0162 |
|            40 |                1600 |    0.132431 |    0.458424 |    0.2889 |
|            52 |                2704 |    0.727744 |    0.718794 |    1.0125 |
|            80 |                6400 |    2.428380 |    0.538877 |    4.5064 |
|           160 |               25600 |   69.134364 |    3.163753 |   21.8520 |
|           320 |              102400 | 1305.199375 |   39.589920 |   32.9680 |
