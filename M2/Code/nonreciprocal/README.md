# General (non-reciprocal-capable) BD-RIS baseline — Member 3

Baseline for the project *Quantifying the Sum-Rate Cost of Reciprocity in BD-RIS-aided MU-MISO*.
Reference: H. Li, S. Shen, B. Clerckx, "Beyond Diagonal RIS: From Transmitting and Reflecting Modes to
Single-, Group-, and Fully-Connected Architectures" (arXiv 2205.02866), **Algorithm 1 + Algorithm 2, reflective mode (K_t = 0)**.

Scattering matrix constraint here: block-diagonal **unitary** (SC: Rg=1, GC: Rg=R/G, FC: Rg=R) and **no** symmetry
constraint (Theta = Theta.' is *not* imposed). The reciprocal method adds exactly that constraint.

## Two modes (both run from the same channel, the same start and the same initial V)
| Mode | Code | What it is | Use it for |
|---|---|---|---|
| **Joint** (paper-faithful) | `li_joint_optimizer.m` | Li Alg. 1: block-coordinate ascent over iota, tau, W (closed form + bisection), Theta (Riemannian CG) | the "general baseline" and paper validation |
| **Matched / V-frozen** | `general_fixedV_optimizer.m` | Aashi's `run_cga_optimizer` with `nu = 0` + plain unitary projection | like-for-like reciprocity gap against Aashi's pipeline |

Aashi's reciprocal pipeline freezes V at the initial MMSE precoder, Li's algorithm updates V every outer iteration.
Compare reciprocal vs general **within the same mode**; the joint-vs-frozen difference is a V-update effect, not a reciprocity effect.

## Files
- `li_joint_optimizer.m` — Algorithm 1 driver. Channels are rescaled internally (noise = 1) so all matrices are O(1); SINR is invariant.
- `li_aux_update.m` (eq. 17-18), `li_precoder_update.m` (eq. 19-20, bisection on lambda),
  `li_theta_terms.m` (eq. 23), `li_theta_objective.m`, `li_theta_gradient.m`, `li_theta_update.m` (Alg. 2: tangent projection, Polak-Ribiere, Armijo, polar retraction).
- `project_unitary_blocks.m` — exact block-unitary projection (counterpart of `enforce_symmetry_unitarity.m`, without symmetrising).
- `init_general_theta.m` — SC diagonal random-phase start (Li Sec. IV-F) or random block-unitary.
- `add_reciprocal_path.m` — puts the sibling `M2/Code/reciprocal` on the path (config, generate_channels, metrics).
- Experiments: `run_general_convergence.m`, `run_matched_monte_carlo.m`, `validate_li_trends.m`. Checks: `test_general_baseline.m`.

## Test configuration (project M2)
N = 2, K = 2, R = 32, Pmax = 20 dBm, Rayleigh, noise -80 dBm, d_BS-RIS = 50 m, d_RIS-user = 2.5 m, zeta_0 = -30 dB, eps = 2.2
(all from Aashi's `config.m`). GC uses Rg = 4 (G = 8 groups, as in Li Fig. 9).
Matched draws: `rng(d,'twister') -> generate_channels -> initialize_theta(R,'SC',R) -> mmse_beamformer`; keep this call order.

## Algorithm settings
Outer loop: up to 300 iterations, stop when |change in sum-rate| < 1e-6 bps/Hz. Inner CG: up to 30 iterations per outer iteration,
Armijo c = 1e-4, halving, polar retraction. These are defaults of `li_joint_optimizer.m` (`opts`).

## Deviations from the paper (simplifications for M2)
1. Reflective mode only (K_t = 0, Phi_t = 0), no direct BS-user link — same as the project model.
2. All G blocks are updated **jointly** on the product manifold, not one group at a time (Alg. 2, steps 4-15). Same stationary points; Aashi's pipeline also updates all blocks jointly.
3. Tangent projection is the standard (block-)unitary one, Gamma - Theta*sym(Theta^H Gamma). Eq. (31) in the arXiv text (`chdiag`) is not the unitary-manifold projection as printed.
4. Retraction: polar factor via SVD (identical to the paper's (Phi + d Xi)(I + d^2 Xi^H Xi)^(-1/2) for tangent Xi).
5. Aashi's code multiplies noise etc. in physical units (hence her tiny Armijo constant 2e-11); this folder normalises internally. Final rates are evaluated with her `calculate_sinr` / `calculate_sumrate` in physical units (checked in the tests).

## Left for M3
Rician fading, sweeps over R, G, Pmax, 50-200 draw Monte-Carlo, runtime/complexity, SC special-case algorithm (Li Alg. 3), hybrid (STAR) mode if the scope grows.

## Validation status
`test_general_baseline.m` — finite-difference gradient check, surrogate vs eq. (21a), monotone ascent, unitarity ~1e-15, power constraint, block structure, metric agreement with Aashi's functions, and FC >= GC >= SC on average.
`validate_li_trends.m` — qualitative trends of Li Fig. 9/10 (rate grows with P; FC > GC > SC in Rayleigh).
