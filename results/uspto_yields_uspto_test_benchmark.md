# uspto_yields_uspto_test: benchmark (mean ± std over seeds / leave-one-out folds)

## split: group:cluster:100

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| lgbm_morgan | 0.131 ± 0.015 | 23.902 ± 0.143 | 19.757 ± 0.116 | 0.357 ± 0.020 | 0.895 ± 0.005 | 70.032 ± 0.405 |

## split: group:cluster:100+gval

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| lgbm_morgan | 0.131 ± 0.014 | 23.896 ± 0.165 | 19.790 ± 0.107 | 0.356 ± 0.017 | 0.906 ± 0.004 | 71.723 ± 0.195 |

## split: group:cluster:100+gkf

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| lgbm_morgan | 0.136 ± 0.016 | 23.829 ± 0.139 | 19.656 ± 0.151 | 0.364 ± 0.019 | 0.906 ± 0.005 | 71.274 ± 0.044 |
