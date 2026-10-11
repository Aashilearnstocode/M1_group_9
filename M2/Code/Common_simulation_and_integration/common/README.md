# Common Simulation and Integration

Shared simulation environment and integration code for the project *Quantifying the Sum-Rate Cost of Reciprocity in Beyond-Diagonal RIS-Aided MU-MISO Systems*.

The purpose of this component is to make the **reciprocal BD-RIS** method (`../reciprocal`) and the **general / non-reciprocal-capable BD-RIS baseline** (`../nonreciprocal`) scientifically comparable. Both methods are fed by the same configuration, the same channel realizations, the same starting point and the same metric functions, so any sum-rate difference between them is caused by the reciprocity constraint and nothing else.

## What was done

1. **Single configuration file.** `configs/default_config.m` holds every setting in one place: antennas (N), users (K), RIS elements (R), number of groups (G), transmit power, noise power, path-loss parameters, fading type, random seed, number of draws, penalty weight and stopping tolerance. Derived quantities (linear power, noise, path loss, block sizes for SC/GC/FC) are recomputed after any override. Field names match `reciprocal/config.m`, so the existing algorithm code accepts it unchanged.
2. **Channel generation.** `common/common_generate_channels.m` builds the BS-to-RIS channel `H_TX` (R x N) and the RIS-to-user channel `H_RX` (K x R) for Rayleigh or Rician fading, with large-scale path loss applied. Random numbers are drawn in the same order as `reciprocal/generate_channels.m`, so a given seed gives identical channels in both code bases.
3. **Noise and path-loss model.** `common/common_pathloss.m` implements `PL = C0 (d/d0)^(-rho)`; noise power is defined in the configuration.
4. **Shared metrics.** `common/common_sinr.m` and `common/common_sumrate.m` are the single SINR and sum-rate implementation; every reported number is recomputed with them. `common/common_mmse_beamformer.m` and `common/common_init_theta.m` provide the common initial precoder and the common starting scattering matrix.
5. **Matched experiment runner.** `common/run_matched_experiment.m` runs, for each draw and each architecture (SC, GC, FC):
   - the reciprocal method (CGA with symmetry penalty, then symmetry/unitarity projection),
   - the general baseline with the precoder frozen (same CGA machinery, no symmetry penalty),
   - the general baseline with joint precoder/scattering design (Li et al., Algorithm 1).
6. **Standard result storage.** `common/save_matched_results.m` writes `.mat`, `.csv` and a configuration record; `common/plot_matched_results.m` produces the summary and convergence figures.
7. **Verification.** `tests/test_common_framework.m` checks the common functions against the reciprocal implementation (same seed gives the same channels, starting matrix, beamformer, SINR and sum-rate), the power constraint, reproducibility, path-loss scaling for Rayleigh and Rician fading, and configuration overrides.
8. **First matched comparison.** A 50-draw Rayleigh run was completed (results below).

## Fairness rules enforced

1. Draw `d` uses `rng(seed_base + d, 'twister')`; with `seed_base = 0` this is the same channel as `rng(d)` in the algorithm folders' own scripts.
2. The order of random calls is fixed: channels, then starting matrix, then initial precoder.
3. Every method starts from the same `Theta_init` and the same initial precoder `V0`.
4. Same dimensions, architecture block sizes, transmit power, noise power and fading for all methods.
5. All final rates are recomputed with the common metric functions.
6. The per-draw CSV stores channel norms, so identical channels can be verified.

## Test configuration

N = 2, K = 2, R = 32, G = 8 (GC block size 4), P_max = 20 dBm, noise = -80 dBm, Rayleigh fading, path loss C0 = -30 dB with exponent 2.2, BS-RIS distance 50 m, RIS-user distance 2.5 m, symmetry-penalty weight nu = 1, CGA tolerance 1e-8, 50 draws (seeds 1 to 50). The full list is in `results/matched_rayleigh_R32_50draws/config_used.txt`.

## Results (50 matched draws, Rayleigh, mean +- standard error, bps/Hz)

Notation: *fixed V* = rate with the frozen initial precoder; *recomputed V* = rate after recomputing the MMSE precoder at the final scattering matrix; *joint* = Li's jointly optimized precoder and scattering matrix.

| Architecture | Reciprocal (fixed V) | General (fixed V) | Gap (general - reciprocal) | Gap in % | General joint |
|---|---|---|---|---|---|
| SC (R_g = 1)  | 8.517 +- 0.267 | 8.517 +- 0.267 | 0.000 +- 0.000 | 0.00 % | 11.290 +- 0.078 |
| GC (R_g = 4)  | 9.909 +- 0.298 | 10.012 +- 0.310 | 0.103 +- 0.026 | 0.88 % | 13.100 +- 0.074 |
| FC (R_g = 32) | 10.310 +- 0.322 | 10.331 +- 0.324 | 0.021 +- 0.012 | 0.17 % | 13.711 +- 0.074 |

With the precoder recomputed at the end: reciprocal 9.145 / 10.449 / 10.857 and general 9.145 / 10.575 / 10.878 for SC / GC / FC, giving gaps of 0.000, 0.126 +- 0.028 (1.07 %) and 0.021 +- 0.012 (0.17 %).

Figures: `results/matched_rayleigh_R32_50draws/matched_sumrate.png` (mean sum-rate with standard error) and `results/matched_rayleigh_R32_50draws/matched_convergence.png` (convergence for draw 1; the SC general curve lies exactly under the SC reciprocal curve).

### Observations

- **SC shows no gap, as expected.** A diagonal scattering matrix is always symmetric, so the reciprocity constraint removes nothing and both methods return identical results.
- **The measured cost of reciprocity is small.** For GC it is about 0.10 bps/Hz (roughly 0.9 %), about four standard errors from zero, with the general design ahead in 48 of 50 draws. For FC it is about 0.02 bps/Hz (roughly 0.2 %), under two standard errors; the median FC gap is essentially zero and the mean is driven by a handful of draws (for example draws 37, 41 and 46, with gaps of 0.27 to 0.38 bps/Hz).
- **The FC gap is not larger than the GC gap.** This is a mean over 50 draws of local-optimum solutions with a frozen precoder, so it should be read as an initial indication, not as a final ranking of architectures. In a few draws the reciprocal solution is slightly better than the general one (2 of 50 for GC, 6 of 50 for FC, at most 0.27 bps/Hz). Both methods stop at local optima, so this reflects optimization differences, not a benefit of reciprocity.
- **Constraints hold.** Reciprocal solutions are symmetric to about 1e-14 (Frobenius norm of Theta - Theta^T) and all solutions are unitary to about 1e-14 or better. The general solutions are genuinely non-symmetric for GC and FC (mean symmetry error about 6 and 2.8), confirming that the baseline uses the extra freedom that the reciprocal design gives up.
- **Joint design is a different comparison.** Li's joint design reaches 11.29 / 13.10 / 13.71 bps/Hz, about 2.8 to 3.4 bps/Hz above the fixed-V results. That difference comes from updating the precoder, not from reciprocity, so it is a reference baseline and should not be read as the reciprocity gap.
- **Iteration limits.** The reciprocal CGA reached its 8000-iteration limit in 4 of 50 GC draws (draws 21, 26, 31 and 34); no SC or FC run did. The joint design reached its 300-outer-iteration limit in almost every draw, so its numbers are slightly conservative.

## Reading the result tables

- `*_fv`: sum-rate with the frozen initial precoder; `*_v2`: sum-rate after recomputing the MMSE precoder; `*_jt`: Li's joint design.
- Reciprocity gap (positive means reciprocity costs rate): `gap_fv = gen_fv - rec_fv` and `gap_v2 = gen_v2 - rec_v2`.
