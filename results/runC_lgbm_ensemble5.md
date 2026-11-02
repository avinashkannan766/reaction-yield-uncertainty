# tdc_buchwald: benchmark (mean ± std over seeds / leave-one-out folds)

## split: random

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| lgbm_ensemble5_morgan | 0.942 ± 0.001 | 6.584 ± 0.097 | 4.653 ± 0.077 | 0.965 ± 0.004 | 0.910 ± 0.011 | 18.768 ± 1.000 |

## split: group:scaffold

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| lgbm_ensemble5_morgan | 0.490 ± 0.010 | 21.324 ± 0.210 | 15.822 ± 0.238 | 0.713 ± 0.006 | 0.711 ± 0.035 | 36.568 ± 2.703 |

## split: group:role:aryl_halide

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| lgbm_ensemble5_morgan | 0.642 ± 0.231 | 14.630 ± 2.779 | 10.613 ± 1.250 | 0.833 ± 0.116 | 0.756 ± 0.019 | 27.879 ± 2.156 |

## split: group:role:additive

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| lgbm_ensemble5_morgan | 0.731 ± 0.055 | 13.952 ± 1.627 | 9.799 ± 0.946 | 0.870 ± 0.044 | 0.739 ± 0.031 | 24.204 ± 0.232 |

## split: loo:role:ligand

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| lgbm_ensemble5_morgan | 0.398 ± 0.800 | 13.910 ± 6.630 | 10.209 ± 4.925 | 0.915 ± 0.065 | 0.751 ± 0.193 | 23.946 ± 3.444 |

## split: loo:role:base

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| lgbm_ensemble5_morgan | 0.557 ± 0.213 | 17.035 ± 3.708 | 11.888 ± 2.142 | 0.849 ± 0.085 | 0.671 ± 0.111 | 24.064 ± 3.548 |

## leave-one-out folds (R² per held-out group)

| split | experiment | held out | n_test | R² | RMSE | coverage |
|---|---|---|---|---|---|---|
| loo:role:ligand | lgbm_ensemble5_morgan | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C(C)(C)C)C( | 986 | 0.756 | 13.89 | 0.70 |
| loo:role:ligand | lgbm_ensemble5_morgan | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2 | 990 | -0.984 | 24.69 | 0.46 |
| loo:role:ligand | lgbm_ensemble5_morgan | COc1ccc(OC)c(P(C(C)(C)C)C(C)(C)C)c1-c1c( | 989 | 0.930 | 7.59 | 0.91 |
| loo:role:ligand | lgbm_ensemble5_morgan | COc1ccc(OC)c(P([C@]23C[C@H]4C[C@H](C[C@H | 990 | 0.889 | 9.47 | 0.94 |
| loo:role:base | lgbm_ensemble5_morgan | CCN=P(N=P(N(C)C)(N(C)C)N(C)C)(N(C)C)N(C) | 1320 | 0.261 | 21.34 | 0.59 |
| loo:role:base | lgbm_ensemble5_morgan | CN(C)C(=NC(C)(C)C)N(C)C | 1317 | 0.752 | 12.29 | 0.83 |
| loo:role:base | lgbm_ensemble5_morgan | CN1CCCN2CCCN=C12 | 1318 | 0.660 | 17.47 | 0.59 |
