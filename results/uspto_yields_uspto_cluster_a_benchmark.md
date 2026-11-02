# uspto_yields_uspto_cluster_a: benchmark (mean ± std over seeds / leave-one-out folds)

## split: group:cluster:100

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| mean_baseline | -0.005 ± 0.006 | 25.699 ± 0.234 | 21.771 ± 0.128 | nan ± nan | 0.911 ± 0.013 | 75.000 ± 0.816 |
| ridge_morgan | 0.071 ± 0.007 | 24.704 ± 0.220 | 20.272 ± 0.146 | 0.318 ± 0.013 | 0.900 ± 0.005 | 71.886 ± 0.465 |
| lgbm_morgan | 0.131 ± 0.015 | 23.902 ± 0.143 | 19.757 ± 0.116 | 0.357 ± 0.020 | 0.895 ± 0.005 | 70.032 ± 0.405 |

## split: group:cluster:100+gval

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| mean_baseline | -0.004 ± 0.005 | 25.689 ± 0.244 | 21.739 ± 0.133 | nan ± nan | 0.913 ± 0.009 | 75.333 ± 0.471 |
| ridge_morgan | 0.066 ± 0.010 | 24.770 ± 0.252 | 20.266 ± 0.188 | 0.318 ± 0.014 | 0.908 ± 0.007 | 72.784 ± 0.448 |
| lgbm_morgan | 0.131 ± 0.014 | 23.896 ± 0.165 | 19.790 ± 0.107 | 0.356 ± 0.017 | 0.906 ± 0.004 | 71.723 ± 0.195 |

## split: group:cluster:100+gkf

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| mean_baseline | -0.005 ± 0.005 | 25.698 ± 0.235 | 21.764 ± 0.124 | nan ± nan | 0.909 ± 0.011 | 75.291 ± 0.275 |
| ridge_morgan | 0.079 ± 0.012 | 24.602 ± 0.202 | 20.161 ± 0.135 | 0.325 ± 0.016 | 0.909 ± 0.005 | 72.712 ± 0.181 |
| lgbm_morgan | 0.136 ± 0.016 | 23.829 ± 0.139 | 19.656 ± 0.151 | 0.364 ± 0.019 | 0.906 ± 0.005 | 71.274 ± 0.044 |
