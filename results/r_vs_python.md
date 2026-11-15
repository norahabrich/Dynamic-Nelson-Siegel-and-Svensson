Python 3.13.18 (NumPy, SciPy) vs R version 4.3.3 (2024-02-29) (base R only), same data and settings. RMSE tables are compared at their saved precision of 0.01 bp.

| Quantity | Python | R | Largest absolute difference |
|---|---:|---:|---:|
| Two-step factors, 348 months × 3 |  |  | 6.4e-14 % |
| In-sample RMSE, Nelson-Siegel (two-step) | 10.44 | 10.44 | 0 bp |
| In-sample RMSE, Svensson (two-step) | 8.40 | 8.40 | 0 bp |
| In-sample RMSE, Nelson-Siegel (Kalman) | 11.15 | 11.16 | 1.0e-02 bp |
| Kalman λ | 0.07777 | 0.07776 | 4.7e-06 |
| Kalman log-likelihood | 3184.552 | 3184.551 | 1.0e-03 |
| Kalman transition matrix A (9 values) |  |  | 7.2e-05 |
| Kalman factor means μ |  |  | 0.04 % |
| Kalman measurement error sd (17 values) |  |  | 4.4e-03 bp |
| Kalman filtered factors, 348 months × 3 |  |  | 7.8e-04 % |
| Out-of-sample RMSE, 5 two-step models and benchmarks (255 values) |  |  | 0 bp |
| Out-of-sample RMSE, Kalman model (51 values) |  |  | 0.21 bp |
| Diebold-Mariano tests with a different 5% verdict |  |  | 0 of 255 |
| Runtime, one Kalman estimation | 102 s | 276 s | |
| Runtime, full pipeline | 8.6 min | 18.8 min | |
