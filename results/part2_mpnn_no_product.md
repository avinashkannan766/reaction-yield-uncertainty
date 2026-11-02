# tdc_buchwald: benchmark (3 seeds, mean ± std)

## split: random

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| mpnn_ensemble3 | 0.816 ± 0.008 | 11.688 ± 0.233 | 8.276 ± 0.116 | 0.894 ± 0.008 | 0.900 ± 0.007 | 33.974 ± 0.762 |

## split: group:scaffold

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| mpnn_ensemble3 | 0.214 ± 0.023 | 26.464 ± 0.397 | 19.889 ± 0.404 | 0.705 ± 0.011 | 0.698 ± 0.042 | 44.929 ± 3.615 |

## split: group:role:aryl_halide

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| mpnn_ensemble3 | 0.613 ± 0.190 | 15.710 ± 4.657 | 11.680 ± 3.456 | 0.769 ± 0.116 | 0.874 ± 0.070 | 41.126 ± 5.063 |

## split: group:role:additive

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| mpnn_ensemble3 | 0.606 ± 0.050 | 16.940 ± 1.342 | 12.366 ± 0.965 | 0.790 ± 0.050 | 0.856 ± 0.022 | 41.369 ± 2.703 |
