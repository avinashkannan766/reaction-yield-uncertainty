# tdc_buchwald: benchmark (mean ± std over seeds / leave-one-out folds)

## split: random

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| rf_morgan | 0.920 ± 0.007 | 7.707 ± 0.283 | 5.231 ± 0.192 | 0.953 ± 0.002 | 0.902 ± 0.018 | 19.939 ± 0.622 |
| lgbm_morgan | 0.941 ± 0.001 | 6.615 ± 0.091 | 4.696 ± 0.078 | 0.965 ± 0.004 | 0.917 ± 0.005 | 20.217 ± 0.924 |
| lgbm_ensemble5_morgan | 0.942 ± 0.001 | 6.584 ± 0.097 | 4.653 ± 0.077 | 0.965 ± 0.004 | 0.910 ± 0.011 | 18.768 ± 1.000 |

## split: group:scaffold

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| rf_morgan | 0.522 ± 0.014 | 20.637 ± 0.293 | 16.178 ± 0.282 | 0.720 ± 0.009 | 0.668 ± 0.058 | 37.104 ± 2.874 |
| lgbm_morgan | 0.474 ± 0.025 | 21.639 ± 0.523 | 15.933 ± 0.327 | 0.712 ± 0.008 | 0.420 ± 0.017 | 19.299 ± 0.555 |
| lgbm_ensemble5_morgan | 0.490 ± 0.010 | 21.324 ± 0.210 | 15.822 ± 0.238 | 0.713 ± 0.006 | 0.711 ± 0.035 | 36.568 ± 2.703 |

## split: group:role:aryl_halide

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| rf_morgan | 0.590 ± 0.301 | 15.269 ± 3.779 | 10.742 ± 2.022 | 0.830 ± 0.124 | 0.773 ± 0.048 | 28.917 ± 4.575 |
| lgbm_morgan | 0.588 ± 0.283 | 15.534 ± 3.286 | 11.120 ± 1.430 | 0.819 ± 0.131 | 0.595 ± 0.039 | 18.471 ± 0.562 |
| lgbm_ensemble5_morgan | 0.642 ± 0.231 | 14.630 ± 2.779 | 10.613 ± 1.250 | 0.833 ± 0.116 | 0.756 ± 0.019 | 27.879 ± 2.156 |

## split: group:role:additive

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| rf_morgan | 0.763 ± 0.069 | 13.031 ± 2.074 | 8.932 ± 1.276 | 0.884 ± 0.044 | 0.759 ± 0.061 | 23.630 ± 0.906 |
| lgbm_morgan | 0.721 ± 0.052 | 14.224 ± 1.537 | 9.955 ± 0.788 | 0.866 ± 0.041 | 0.637 ± 0.012 | 18.611 ± 0.823 |
| lgbm_ensemble5_morgan | 0.731 ± 0.055 | 13.952 ± 1.627 | 9.799 ± 0.946 | 0.870 ± 0.044 | 0.739 ± 0.031 | 24.204 ± 0.232 |

## split: group:role:aryl_halide+gval

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| rf_morgan | 0.614 ± 0.076 | 16.394 ± 1.619 | 11.940 ± 1.466 | 0.867 ± 0.042 | 0.842 ± 0.067 | 39.200 ± 3.235 |
| lgbm_morgan | 0.541 ± 0.169 | 17.856 ± 4.326 | 14.315 ± 4.889 | 0.839 ± 0.034 | 0.896 ± 0.093 | 56.407 ± 4.497 |
| lgbm_ensemble5_morgan | 0.536 ± 0.180 | 17.912 ± 4.478 | 14.201 ± 5.015 | 0.849 ± 0.041 | 0.885 ± 0.033 | 55.119 ± 16.738 |

## split: group:role:additive+gval

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| rf_morgan | 0.741 ± 0.091 | 13.573 ± 2.619 | 9.310 ± 1.749 | 0.877 ± 0.041 | 0.866 ± 0.090 | 37.816 ± 8.930 |
| lgbm_morgan | 0.729 ± 0.088 | 13.914 ± 2.487 | 9.805 ± 1.811 | 0.863 ± 0.049 | 0.864 ± 0.096 | 43.471 ± 15.108 |
| lgbm_ensemble5_morgan | 0.720 ± 0.086 | 14.152 ± 2.430 | 10.093 ± 1.899 | 0.863 ± 0.047 | 0.863 ± 0.091 | 42.196 ± 11.467 |

## split: loo:role:ligand

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| rf_morgan | 0.269 ± 0.833 | 16.499 ± 6.901 | 12.334 ± 5.383 | 0.903 ± 0.064 | 0.645 ± 0.231 | 23.098 ± 3.733 |
| lgbm_morgan | 0.406 ± 0.799 | 13.716 ± 6.635 | 10.018 ± 4.895 | 0.914 ± 0.064 | 0.666 ± 0.171 | 18.402 ± 1.842 |
| lgbm_ensemble5_morgan | 0.398 ± 0.800 | 13.910 ± 6.630 | 10.209 ± 4.925 | 0.915 ± 0.065 | 0.751 ± 0.193 | 23.946 ± 3.444 |

## split: loo:role:base

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| rf_morgan | 0.553 ± 0.146 | 17.439 ± 3.003 | 12.172 ± 1.873 | 0.837 ± 0.071 | 0.572 ± 0.112 | 20.651 ± 3.041 |
| lgbm_morgan | 0.520 ± 0.257 | 17.603 ± 4.338 | 12.291 ± 2.542 | 0.844 ± 0.088 | 0.583 ± 0.075 | 18.772 ± 2.304 |
| lgbm_ensemble5_morgan | 0.557 ± 0.213 | 17.035 ± 3.708 | 11.888 ± 2.142 | 0.849 ± 0.085 | 0.671 ± 0.111 | 24.064 ± 3.548 |

## leave-one-out folds (R² per held-out group)

| split | experiment | held out | n_test | R² | RMSE | coverage |
|---|---|---|---|---|---|---|
| loo:role:ligand | rf_morgan | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C(C)(C)C)C( | 986 | 0.482 | 20.24 | 0.46 |
| loo:role:ligand | rf_morgan | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2 | 990 | -1.146 | 25.67 | 0.37 |
| loo:role:ligand | rf_morgan | COc1ccc(OC)c(P(C(C)(C)C)C(C)(C)C)c1-c1c( | 989 | 0.924 | 7.92 | 0.89 |
| loo:role:ligand | rf_morgan | COc1ccc(OC)c(P([C@]23C[C@H]4C[C@H](C[C@H | 990 | 0.816 | 12.17 | 0.86 |
| loo:role:ligand | lgbm_morgan | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C(C)(C)C)C( | 986 | 0.773 | 13.40 | 0.62 |
| loo:role:ligand | lgbm_morgan | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2 | 990 | -0.975 | 24.63 | 0.41 |
| loo:role:ligand | lgbm_morgan | COc1ccc(OC)c(P(C(C)(C)C)C(C)(C)C)c1-c1c( | 989 | 0.925 | 7.85 | 0.87 |
| loo:role:ligand | lgbm_morgan | COc1ccc(OC)c(P([C@]23C[C@H]4C[C@H](C[C@H | 990 | 0.900 | 8.98 | 0.77 |
| loo:role:ligand | lgbm_ensemble5_morgan | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C(C)(C)C)C( | 986 | 0.756 | 13.89 | 0.70 |
| loo:role:ligand | lgbm_ensemble5_morgan | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2 | 990 | -0.984 | 24.69 | 0.46 |
| loo:role:ligand | lgbm_ensemble5_morgan | COc1ccc(OC)c(P(C(C)(C)C)C(C)(C)C)c1-c1c( | 989 | 0.930 | 7.59 | 0.91 |
| loo:role:ligand | lgbm_ensemble5_morgan | COc1ccc(OC)c(P([C@]23C[C@H]4C[C@H](C[C@H | 990 | 0.889 | 9.47 | 0.94 |
| loo:role:base | rf_morgan | CCN=P(N=P(N(C)C)(N(C)C)N(C)C)(N(C)C)N(C) | 1320 | 0.360 | 19.86 | 0.53 |
| loo:role:base | rf_morgan | CN(C)C(=NC(C)(C)C)N(C)C | 1317 | 0.713 | 13.21 | 0.73 |
| loo:role:base | rf_morgan | CN1CCCN2CCCN=C12 | 1318 | 0.587 | 19.25 | 0.46 |
| loo:role:base | lgbm_morgan | CCN=P(N=P(N(C)C)(N(C)C)N(C)C)(N(C)C)N(C) | 1320 | 0.163 | 22.71 | 0.50 |
| loo:role:base | lgbm_morgan | CN(C)C(=NC(C)(C)C)N(C)C | 1317 | 0.759 | 12.10 | 0.68 |
| loo:role:base | lgbm_morgan | CN1CCCN2CCCN=C12 | 1318 | 0.639 | 18.00 | 0.57 |
| loo:role:base | lgbm_ensemble5_morgan | CCN=P(N=P(N(C)C)(N(C)C)N(C)C)(N(C)C)N(C) | 1320 | 0.261 | 21.34 | 0.59 |
| loo:role:base | lgbm_ensemble5_morgan | CN(C)C(=NC(C)(C)C)N(C)C | 1317 | 0.752 | 12.29 | 0.83 |
| loo:role:base | lgbm_ensemble5_morgan | CN1CCCN2CCCN=C12 | 1318 | 0.660 | 17.47 | 0.59 |
