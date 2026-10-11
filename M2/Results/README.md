# Results: quick check, Rayleigh, R = 32, 3 draws

Short run of the matched comparison (seeds 1 to 3) used to check that the full pipeline works end to end before the 50-draw run. With only 3 draws the standard errors are large, so these numbers should not be used for conclusions; the report results are in the 50-draw folder (`matched_rayleigh_R32_50draws`).

## Files

| File | Content |
|---|---|
| `matched_results.mat` | Full result structure `res`: per-draw rates, iteration counts, errors, run times, convergence histories and the configuration used |
| `per_draw.csv` | One row per draw and architecture (SC, GC, FC); column dictionary below |
| `summary.csv` | Mean, standard error (`se`) and number of draws (`n`) for each metric and architecture |
| `config_used.txt` | Every configuration value used for this run |
| `matched_sumrate.png` | Mean sum-rate with standard-error bars: reciprocal vs general, fixed and recomputed precoder |
| `matched_convergence.png` | Sum-rate versus CGA iteration for draw 1: reciprocal (solid) vs general (dashed). The SC general curve lies exactly under the SC reciprocal curve |

## Column dictionary (`per_draw.csv`)

| Column | Meaning |
|---|---|
| `draw`, `seed` | Draw index and the random seed used (`seed = draw`) |
| `arch`, `Rg` | Architecture (SC, GC, FC) and RIS block size (1, 4, 32) |
| `H_TX_fro`, `H_RX_fro` | Frobenius norms of the two channels, a fingerprint showing that all methods used the same channel |
| `rate_init` | Sum-rate of the common starting point (diagonal scattering matrix and MMSE precoder), bps/Hz |
| `rec_fv`, `gen_fv` | Final sum-rate of the reciprocal / general method with the precoder frozen at its initial value |
| `rec_v2`, `gen_v2` | Final sum-rate after recomputing the MMSE precoder at the final scattering matrix |
| `gen_jt` | Final sum-rate of the general method with joint precoder and scattering-matrix design |
| `gap_fv`, `gap_v2` | `gen - rec` for the fixed and recomputed precoder; positive means reciprocity costs rate |
| `gap_fv_pct` | `gap_fv` as a percentage of `gen_fv` |
| `rec_iters`, `gen_iters`, `jt_iters` | Iterations used (CGA iterations; outer iterations for the joint design) |
| `rec_sym_err`, `gen_sym_err` | Frobenius norm of `Theta - Theta^T` for the final scattering matrix |
| `rec_uni_err`, `gen_uni_err` | Frobenius norm of `Theta^H Theta - I` for the final scattering matrix |
| `rec_time_s`, `gen_time_s`, `jt_time_s` | Run time in seconds |

## Configuration

N = 2, K = 2, R = 32, G = 8 (GC block size 4), P_max = 20 dBm, noise = -80 dBm, Rayleigh fading, path loss C0 = -30 dB with exponent 2.2, BS-RIS distance 50 m, RIS-user distance 2.5 m, symmetry-penalty weight nu = 1, CGA tolerance 1e-8, joint design enabled. Full list in `config_used.txt`.

## Results (mean +- standard error over 3 draws, bps/Hz)

| Architecture | Reciprocal (fixed V) | General (fixed V) | Gap (general - reciprocal) | General joint |
|---|---|---|---|---|
| SC (R_g = 1)  | 6.085 +- 1.007 | 6.085 +- 1.007 | 0.000 +- 0.000 | 11.288 +- 0.089 |
| GC (R_g = 4)  | 7.169 +- 1.146 | 7.170 +- 1.146 | 0.001 +- 0.001 | 12.979 +- 0.110 |
| FC (R_g = 32) | 7.322 +- 1.156 | 7.322 +- 1.156 | 0.000 +- 0.000 | 13.613 +- 0.155 |

With the precoder recomputed: reciprocal 7.285 / 8.150 / 8.274 and general 7.285 / 8.162 / 8.274 for SC / GC / FC.

![Mean sum-rate](matched_sumrate.png)

![Convergence, draw 1](matched_convergence.png)

## Notes

- The per-draw values for draws 1 to 3 are identical to draws 1 to 3 of the 50-draw run, which confirms that matched seeds reproduce the same channels and results.
- The large error bars come from the small number of draws and the wide spread of channel realizations, not from a problem in the comparison.
- The general and reciprocal results are almost identical in this sample, and the ordering SC < GC < FC and the large joint-design advantage match the 50-draw results.
