# tdc_buchwald_gkf: benchmark (mean ± std over seeds / leave-one-out folds)

## split: group:scaffold+gkf

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| mean_baseline | -0.017 ± 0.001 | 30.101 ± 0.012 | 26.176 ± 0.005 | nan ± nan | 0.919 ± 0.002 | 83.403 ± 0.171 |
| ridge_morgan | 0.504 ± 0.001 | 21.029 ± 0.031 | 17.530 ± 0.038 | 0.743 ± 0.001 | 0.992 ± 0.002 | 87.164 ± 0.317 |
| rf_morgan | 0.513 ± 0.010 | 20.818 ± 0.217 | 16.363 ± 0.166 | 0.716 ± 0.007 | 0.996 ± 0.001 | 94.215 ± 0.243 |
| lgbm_morgan | 0.476 ± 0.059 | 21.572 ± 1.203 | 15.741 ± 0.748 | 0.713 ± 0.023 | 0.972 ± 0.009 | 84.687 ± 1.369 |
| lgbm_ensemble5_morgan | 0.484 ± 0.012 | 21.437 ± 0.247 | 15.875 ± 0.205 | 0.711 ± 0.005 | 0.989 ± 0.005 | 91.402 ± 3.114 |
| lgbm_morgan_desc | 0.462 ± 0.059 | 21.861 ± 1.222 | 15.807 ± 0.728 | 0.707 ± 0.028 | 0.975 ± 0.008 | 86.334 ± 1.144 |
| mpnn_ensemble3 | 0.216 ± 0.034 | 26.421 ± 0.577 | 19.558 ± 0.519 | 0.693 ± 0.008 | 0.869 ± 0.055 | 63.802 ± 4.651 |
| mpnn_noproduct | 0.225 ± 0.019 | 26.276 ± 0.323 | 19.618 ± 0.352 | 0.708 ± 0.004 | 0.901 ± 0.037 | 66.552 ± 5.395 |
| id_embed_noproduct | 0.323 ± 0.003 | 24.563 ± 0.057 | 19.890 ± 0.041 | 0.578 ± 0.002 | 0.877 ± 0.092 | 65.531 ± 10.253 |
| mpnn_transformer_noproduct | 0.139 ± 0.039 | 27.685 ± 0.620 | 19.889 ± 0.487 | 0.653 ± 0.014 | 0.822 ± 0.071 | 61.536 ± 10.307 |

## split: group:role:aryl_halide+gkf

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| mean_baseline | -0.030 ± 0.040 | 27.081 ± 2.838 | 23.334 ± 2.656 | nan ± nan | 0.918 ± 0.044 | 76.432 ± 1.462 |
| ridge_morgan | 0.500 ± 0.145 | 18.420 ± 1.781 | 15.138 ± 1.705 | 0.756 ± 0.111 | 0.915 ± 0.053 | 56.242 ± 2.831 |
| rf_morgan | 0.586 ± 0.293 | 15.470 ± 3.578 | 10.887 ± 1.773 | 0.827 ± 0.116 | 0.990 ± 0.006 | 61.744 ± 11.263 |
| lgbm_morgan | 0.592 ± 0.284 | 15.400 ± 3.401 | 10.948 ± 1.538 | 0.822 ± 0.126 | 0.924 ± 0.067 | 50.246 ± 3.899 |
| lgbm_ensemble5_morgan | 0.646 ± 0.217 | 14.641 ± 2.553 | 10.629 ± 1.156 | 0.835 ± 0.110 | 0.977 ± 0.017 | 55.099 ± 10.647 |
| lgbm_morgan_desc | 0.669 ± 0.209 | 14.074 ± 2.813 | 9.973 ± 1.525 | 0.838 ± 0.103 | 0.938 ± 0.054 | 52.349 ± 4.683 |
| mpnn_ensemble3 | 0.625 ± 0.190 | 15.416 ± 4.821 | 11.220 ± 3.682 | 0.776 ± 0.118 | 0.909 ± 0.095 | 48.257 ± 7.606 |
| mpnn_noproduct | 0.603 ± 0.201 | 15.748 ± 4.637 | 11.531 ± 3.521 | 0.764 ± 0.119 | 0.940 ± 0.049 | 52.005 ± 5.329 |
| id_embed_noproduct | 0.269 ± 0.123 | 22.846 ± 3.804 | 17.814 ± 3.393 | 0.529 ± 0.112 | 0.906 ± 0.080 | 68.974 ± 5.068 |
| mpnn_transformer_noproduct | 0.666 ± 0.179 | 14.487 ± 4.876 | 10.259 ± 3.199 | 0.804 ± 0.104 | 0.914 ± 0.088 | 46.311 ± 7.826 |

## split: group:role:additive+gkf

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| mean_baseline | -0.066 ± 0.056 | 27.856 ± 0.541 | 23.698 ± 0.444 | nan ± nan | 0.892 ± 0.010 | 74.998 ± 0.696 |
| ridge_morgan | 0.553 ± 0.048 | 18.048 ± 1.247 | 14.542 ± 0.955 | 0.778 ± 0.047 | 0.912 ± 0.032 | 53.419 ± 2.330 |
| rf_morgan | 0.770 ± 0.070 | 12.825 ± 2.131 | 8.753 ± 1.363 | 0.886 ± 0.044 | 0.950 ± 0.037 | 44.663 ± 4.141 |
| lgbm_morgan | 0.724 ± 0.056 | 14.124 ± 1.642 | 9.839 ± 0.877 | 0.866 ± 0.046 | 0.926 ± 0.045 | 46.992 ± 4.303 |
| lgbm_ensemble5_morgan | 0.731 ± 0.059 | 13.954 ± 1.715 | 9.786 ± 0.983 | 0.869 ± 0.046 | 0.927 ± 0.035 | 44.541 ± 3.656 |
| lgbm_morgan_desc | 0.728 ± 0.067 | 13.996 ± 1.914 | 9.793 ± 0.971 | 0.872 ± 0.048 | 0.938 ± 0.044 | 48.988 ± 4.605 |
| mpnn_ensemble3 | 0.629 ± 0.065 | 16.415 ± 1.683 | 11.852 ± 1.258 | 0.800 ± 0.058 | 0.920 ± 0.068 | 52.418 ± 10.720 |
| mpnn_noproduct | 0.612 ± 0.066 | 16.782 ± 1.676 | 12.188 ± 1.196 | 0.794 ± 0.053 | 0.908 ± 0.074 | 52.063 ± 9.615 |
| id_embed_noproduct | 0.591 ± 0.062 | 17.232 ± 1.565 | 12.455 ± 0.995 | 0.819 ± 0.064 | 0.896 ± 0.042 | 46.156 ± 2.813 |
| mpnn_transformer_noproduct | 0.740 ± 0.055 | 13.701 ± 1.665 | 9.373 ± 0.878 | 0.864 ± 0.039 | 0.945 ± 0.030 | 51.313 ± 6.895 |

## split: loo:role:ligand+gkf

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| mean_baseline | -0.289 ± 0.393 | 28.103 ± 2.141 | 24.196 ± 1.509 | nan ± nan | 0.904 ± 0.059 | 76.688 ± 2.050 |
| ridge_morgan | 0.124 ± 0.792 | 20.287 ± 3.813 | 16.779 ± 3.577 | 0.825 ± 0.062 | 0.859 ± 0.160 | 57.286 ± 8.100 |
| rf_morgan | 0.290 ± 0.850 | 15.878 ± 6.986 | 11.854 ± 5.381 | 0.908 ± 0.066 | 0.881 ± 0.180 | 53.688 ± 13.907 |
| lgbm_morgan | 0.384 ± 0.832 | 13.882 ± 6.847 | 10.140 ± 5.080 | 0.915 ± 0.065 | 0.868 ± 0.190 | 46.096 ± 11.337 |
| lgbm_ensemble5_morgan | 0.388 ± 0.821 | 13.890 ± 6.811 | 10.185 ± 5.042 | 0.916 ± 0.066 | 0.894 ± 0.159 | 50.785 ± 11.269 |
| lgbm_morgan_desc | 0.421 ± 0.789 | 13.406 ± 6.621 | 9.670 ± 4.886 | 0.916 ± 0.064 | 0.871 ± 0.187 | 46.100 ± 11.467 |
| mpnn_ensemble3 | 0.149 ± 0.715 | 19.904 ± 6.631 | 14.753 ± 5.043 | 0.838 ± 0.068 | 0.801 ± 0.147 | 45.035 ± 5.791 |
| mpnn_noproduct | 0.162 ± 0.693 | 19.856 ± 6.436 | 14.726 ± 4.861 | 0.843 ± 0.064 | 0.851 ± 0.138 | 52.476 ± 9.753 |
| id_embed_noproduct | 0.307 ± 0.881 | 15.467 ± 6.538 | 11.202 ± 4.676 | 0.904 ± 0.066 | 0.848 ± 0.180 | 42.373 ± 9.399 |
| mpnn_transformer_noproduct | 0.191 ± 0.824 | 18.402 ± 6.333 | 13.361 ± 4.731 | 0.864 ± 0.060 | 0.887 ± 0.120 | 57.237 ± 9.643 |

## split: loo:role:base+gkf

| experiment | R² | RMSE | MAE | Spearman | 90% coverage | interval width |
|---|---|---|---|---|---|---|
| mean_baseline | -0.106 ± 0.065 | 27.898 ± 3.298 | 23.923 ± 2.279 | nan ± nan | 0.901 ± 0.120 | 77.624 ± 7.942 |
| ridge_morgan | 0.495 ± 0.156 | 18.583 ± 3.270 | 14.897 ± 2.607 | 0.802 ± 0.083 | 0.912 ± 0.065 | 58.415 ± 8.384 |
| rf_morgan | 0.553 ± 0.149 | 17.427 ± 2.927 | 12.142 ± 1.817 | 0.838 ± 0.073 | 0.959 ± 0.035 | 67.498 ± 10.291 |
| lgbm_morgan | 0.528 ± 0.238 | 17.533 ± 3.992 | 12.172 ± 2.289 | 0.844 ± 0.084 | 0.906 ± 0.068 | 54.730 ± 7.888 |
| lgbm_ensemble5_morgan | 0.560 ± 0.213 | 16.965 ± 3.709 | 11.817 ± 2.126 | 0.849 ± 0.085 | 0.942 ± 0.043 | 58.888 ± 6.201 |
| lgbm_morgan_desc | 0.549 ± 0.229 | 17.117 ± 3.777 | 11.962 ± 2.119 | 0.851 ± 0.087 | 0.915 ± 0.067 | 54.955 ± 8.170 |
| mpnn_ensemble3 | 0.555 ± 0.107 | 17.475 ± 2.011 | 12.726 ± 1.241 | 0.809 ± 0.062 | 0.884 ± 0.074 | 47.428 ± 9.455 |
| mpnn_noproduct | 0.545 ± 0.105 | 17.690 ± 1.930 | 12.876 ± 1.255 | 0.803 ± 0.057 | 0.858 ± 0.125 | 48.840 ± 17.219 |
| id_embed_noproduct | 0.490 ± 0.214 | 18.427 ± 3.533 | 12.767 ± 1.910 | 0.835 ± 0.077 | 0.825 ± 0.106 | 43.497 ± 10.727 |
| mpnn_transformer_noproduct | 0.281 ± 0.271 | 22.446 ± 6.356 | 16.209 ± 5.344 | 0.753 ± 0.034 | 0.754 ± 0.112 | 41.866 ± 8.645 |

## leave-one-out folds (R² per held-out group)

| split | experiment | held out | n_test | R² | RMSE | coverage |
|---|---|---|---|---|---|---|
| loo:role:ligand+gkf | mean_baseline | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C(C)(C)C)C( | 986 | -0.147 | 30.13 | 0.85 |
| loo:role:ligand+gkf | mean_baseline | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2 | 990 | -0.963 | 24.55 | 1.00 |
| loo:role:ligand+gkf | mean_baseline | COc1ccc(OC)c(P(C(C)(C)C)C(C)(C)C)c1-c1c( | 989 | -0.045 | 29.35 | 0.86 |
| loo:role:ligand+gkf | mean_baseline | COc1ccc(OC)c(P([C@]23C[C@H]4C[C@H](C[C@H | 990 | -0.000 | 28.38 | 0.91 |
| loo:role:ligand+gkf | ridge_morgan | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C(C)(C)C)C( | 986 | 0.450 | 20.87 | 0.89 |
| loo:role:ligand+gkf | ridge_morgan | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2 | 990 | -1.241 | 26.24 | 0.59 |
| loo:role:ligand+gkf | ridge_morgan | COc1ccc(OC)c(P(C(C)(C)C)C(C)(C)C)c1-c1c( | 989 | 0.678 | 16.28 | 0.98 |
| loo:role:ligand+gkf | ridge_morgan | COc1ccc(OC)c(P([C@]23C[C@H]4C[C@H](C[C@H | 990 | 0.608 | 17.76 | 0.97 |
| loo:role:ligand+gkf | rf_morgan | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C(C)(C)C)C( | 986 | 0.561 | 18.64 | 0.95 |
| loo:role:ligand+gkf | rf_morgan | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2 | 990 | -1.163 | 25.78 | 0.57 |
| loo:role:ligand+gkf | rf_morgan | COc1ccc(OC)c(P(C(C)(C)C)C(C)(C)C)c1-c1c( | 989 | 0.934 | 7.37 | 1.00 |
| loo:role:ligand+gkf | rf_morgan | COc1ccc(OC)c(P([C@]23C[C@H]4C[C@H](C[C@H | 990 | 0.829 | 11.73 | 1.00 |
| loo:role:ligand+gkf | lgbm_morgan | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C(C)(C)C)C( | 986 | 0.764 | 13.66 | 0.93 |
| loo:role:ligand+gkf | lgbm_morgan | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2 | 990 | -1.053 | 25.11 | 0.54 |
| loo:role:ligand+gkf | lgbm_morgan | COc1ccc(OC)c(P(C(C)(C)C)C(C)(C)C)c1-c1c( | 989 | 0.928 | 7.71 | 1.00 |
| loo:role:ligand+gkf | lgbm_morgan | COc1ccc(OC)c(P([C@]23C[C@H]4C[C@H](C[C@H | 990 | 0.898 | 9.05 | 1.00 |
| loo:role:ligand+gkf | lgbm_ensemble5_morgan | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C(C)(C)C)C( | 986 | 0.760 | 13.79 | 0.96 |
| loo:role:ligand+gkf | lgbm_ensemble5_morgan | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2 | 990 | -1.030 | 24.97 | 0.62 |
| loo:role:ligand+gkf | lgbm_ensemble5_morgan | COc1ccc(OC)c(P(C(C)(C)C)C(C)(C)C)c1-c1c( | 989 | 0.935 | 7.30 | 0.99 |
| loo:role:ligand+gkf | lgbm_ensemble5_morgan | COc1ccc(OC)c(P([C@]23C[C@H]4C[C@H](C[C@H | 990 | 0.888 | 9.49 | 1.00 |
| loo:role:ligand+gkf | lgbm_morgan_desc | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C(C)(C)C)C( | 986 | 0.802 | 12.53 | 0.94 |
| loo:role:ligand+gkf | lgbm_morgan_desc | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2 | 990 | -0.943 | 24.43 | 0.55 |
| loo:role:ligand+gkf | lgbm_morgan_desc | COc1ccc(OC)c(P(C(C)(C)C)C(C)(C)C)c1-c1c( | 989 | 0.932 | 7.47 | 1.00 |
| loo:role:ligand+gkf | lgbm_morgan_desc | COc1ccc(OC)c(P([C@]23C[C@H]4C[C@H](C[C@H | 990 | 0.895 | 9.19 | 1.00 |
| loo:role:ligand+gkf | mpnn_ensemble3 | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C(C)(C)C)C( | 986 | 0.733 | 14.55 | 0.92 |
| loo:role:ligand+gkf | mpnn_ensemble3 | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2 | 990 | -0.965 | 24.57 | 0.61 |
| loo:role:ligand+gkf | mpnn_ensemble3 | COc1ccc(OC)c(P(C(C)(C)C)C(C)(C)C)c1-c1c( | 989 | 0.815 | 12.34 | 0.96 |
| loo:role:ligand+gkf | mpnn_ensemble3 | COc1ccc(OC)c(P([C@]23C[C@H]4C[C@H](C[C@H | 990 | 0.015 | 28.16 | 0.72 |
| loo:role:ligand+gkf | mpnn_noproduct | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C(C)(C)C)C( | 986 | 0.706 | 15.24 | 0.96 |
| loo:role:ligand+gkf | mpnn_noproduct | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2 | 990 | -0.921 | 24.29 | 0.64 |
| loo:role:ligand+gkf | mpnn_noproduct | COc1ccc(OC)c(P(C(C)(C)C)C(C)(C)C)c1-c1c( | 989 | 0.824 | 12.05 | 0.99 |
| loo:role:ligand+gkf | mpnn_noproduct | COc1ccc(OC)c(P([C@]23C[C@H]4C[C@H](C[C@H | 990 | 0.037 | 27.84 | 0.81 |
| loo:role:ligand+gkf | id_embed_noproduct | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C(C)(C)C)C( | 986 | 0.690 | 15.65 | 0.87 |
| loo:role:ligand+gkf | id_embed_noproduct | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2 | 990 | -1.214 | 26.08 | 0.55 |
| loo:role:ligand+gkf | id_embed_noproduct | COc1ccc(OC)c(P(C(C)(C)C)C(C)(C)C)c1-c1c( | 989 | 0.870 | 10.35 | 0.99 |
| loo:role:ligand+gkf | id_embed_noproduct | COc1ccc(OC)c(P([C@]23C[C@H]4C[C@H](C[C@H | 990 | 0.881 | 9.79 | 0.99 |
| loo:role:ligand+gkf | mpnn_transformer_noproduct | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C(C)(C)C)C( | 986 | 0.792 | 12.82 | 0.97 |
| loo:role:ligand+gkf | mpnn_transformer_noproduct | CC(C)c1cc(C(C)C)c(-c2ccccc2P(C2CCCCC2)C2 | 990 | -1.193 | 25.96 | 0.68 |
| loo:role:ligand+gkf | mpnn_transformer_noproduct | COc1ccc(OC)c(P(C(C)(C)C)C(C)(C)C)c1-c1c( | 989 | 0.840 | 11.49 | 0.98 |
| loo:role:ligand+gkf | mpnn_transformer_noproduct | COc1ccc(OC)c(P([C@]23C[C@H]4C[C@H](C[C@H | 990 | 0.323 | 23.34 | 0.90 |
| loo:role:base+gkf | mean_baseline | CCN=P(N=P(N(C)C)(N(C)C)N(C)C)(N(C)C)N(C) | 1320 | -0.122 | 26.29 | 0.99 |
| loo:role:base+gkf | mean_baseline | CN(C)C(=NC(C)(C)C)N(C)C | 1317 | -0.020 | 24.91 | 0.98 |
| loo:role:base+gkf | mean_baseline | CN1CCCN2CCCN=C12 | 1318 | -0.177 | 32.49 | 0.73 |
| loo:role:base+gkf | ridge_morgan | CCN=P(N=P(N(C)C)(N(C)C)N(C)C)(N(C)C)N(C) | 1320 | 0.298 | 20.79 | 0.89 |
| loo:role:base+gkf | ridge_morgan | CN(C)C(=NC(C)(C)C)N(C)C | 1317 | 0.680 | 13.96 | 1.00 |
| loo:role:base+gkf | ridge_morgan | CN1CCCN2CCCN=C12 | 1318 | 0.509 | 21.00 | 0.85 |
| loo:role:base+gkf | rf_morgan | CCN=P(N=P(N(C)C)(N(C)C)N(C)C)(N(C)C)N(C) | 1320 | 0.352 | 19.98 | 0.91 |
| loo:role:base+gkf | rf_morgan | CN(C)C(=NC(C)(C)C)N(C)C | 1317 | 0.708 | 13.33 | 1.00 |
| loo:role:base+gkf | rf_morgan | CN1CCCN2CCCN=C12 | 1318 | 0.599 | 18.97 | 0.96 |
| loo:role:base+gkf | lgbm_morgan | CCN=P(N=P(N(C)C)(N(C)C)N(C)C)(N(C)C)N(C) | 1320 | 0.197 | 22.24 | 0.83 |
| loo:role:base+gkf | lgbm_morgan | CN(C)C(=NC(C)(C)C)N(C)C | 1317 | 0.744 | 12.48 | 0.99 |
| loo:role:base+gkf | lgbm_morgan | CN1CCCN2CCCN=C12 | 1318 | 0.644 | 17.88 | 0.90 |
| loo:role:base+gkf | lgbm_ensemble5_morgan | CCN=P(N=P(N(C)C)(N(C)C)N(C)C)(N(C)C)N(C) | 1320 | 0.263 | 21.31 | 0.89 |
| loo:role:base+gkf | lgbm_ensemble5_morgan | CN(C)C(=NC(C)(C)C)N(C)C | 1317 | 0.753 | 12.25 | 0.99 |
| loo:role:base+gkf | lgbm_ensemble5_morgan | CN1CCCN2CCCN=C12 | 1318 | 0.665 | 17.34 | 0.95 |
| loo:role:base+gkf | lgbm_morgan_desc | CCN=P(N=P(N(C)C)(N(C)C)N(C)C)(N(C)C)N(C) | 1320 | 0.227 | 21.83 | 0.83 |
| loo:role:base+gkf | lgbm_morgan_desc | CN(C)C(=NC(C)(C)C)N(C)C | 1317 | 0.740 | 12.58 | 1.00 |
| loo:role:base+gkf | lgbm_morgan_desc | CN1CCCN2CCCN=C12 | 1318 | 0.680 | 16.95 | 0.91 |
| loo:role:base+gkf | mpnn_ensemble3 | CCN=P(N=P(N(C)C)(N(C)C)N(C)C)(N(C)C)N(C) | 1320 | 0.405 | 19.14 | 0.78 |
| loo:role:base+gkf | mpnn_ensemble3 | CN(C)C(=NC(C)(C)C)N(C)C | 1317 | 0.647 | 14.65 | 0.93 |
| loo:role:base+gkf | mpnn_ensemble3 | CN1CCCN2CCCN=C12 | 1318 | 0.613 | 18.64 | 0.94 |
| loo:role:base+gkf | mpnn_noproduct | CCN=P(N=P(N(C)C)(N(C)C)N(C)C)(N(C)C)N(C) | 1320 | 0.397 | 19.27 | 0.68 |
| loo:role:base+gkf | mpnn_noproduct | CN(C)C(=NC(C)(C)C)N(C)C | 1317 | 0.631 | 14.97 | 0.93 |
| loo:role:base+gkf | mpnn_noproduct | CN1CCCN2CCCN=C12 | 1318 | 0.605 | 18.82 | 0.96 |
| loo:role:base+gkf | id_embed_noproduct | CCN=P(N=P(N(C)C)(N(C)C)N(C)C)(N(C)C)N(C) | 1320 | 0.192 | 22.31 | 0.71 |
| loo:role:base+gkf | id_embed_noproduct | CN(C)C(=NC(C)(C)C)N(C)C | 1317 | 0.689 | 13.76 | 0.97 |
| loo:role:base+gkf | id_embed_noproduct | CN1CCCN2CCCN=C12 | 1318 | 0.589 | 19.21 | 0.80 |
| loo:role:base+gkf | mpnn_transformer_noproduct | CCN=P(N=P(N(C)C)(N(C)C)N(C)C)(N(C)C)N(C) | 1320 | 0.383 | 19.49 | 0.73 |
| loo:role:base+gkf | mpnn_transformer_noproduct | CN(C)C(=NC(C)(C)C)N(C)C | 1317 | 0.548 | 16.57 | 0.90 |
| loo:role:base+gkf | mpnn_transformer_noproduct | CN1CCCN2CCCN=C12 | 1318 | -0.090 | 31.28 | 0.63 |
