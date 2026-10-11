# M2 Implementation Plan — Group 9 BD-RIS Reciprocity-Gap Project

**Project:** Quantifying the Sum-Rate Cost of Reciprocity in Beyond-Diagonal RIS-Aided MU-MISO Systems  
**Purpose:** Provide an actionable build order, function contracts, validation tests T1–T14, schedule, fair-comparison protocol, and video outline.  
**Status:** Planning artifact. A test listed here is a required test, not evidence that it has already passed. Do not report a solver as implemented or a result as reproduced until the code has been run and its outputs saved.

---

## 1. M2 objective and boundary

M2 should establish the technical machinery and produce an initial reproducible comparison. It does **not** need to complete every M3 sweep or the final Pareto analysis.

The project compares:

- **General/non-reciprocal-capable design:** optimize the BD-RIS scattering matrix over the relevant unitary architecture.
- **Reciprocal design:** optimize over the same architecture while enforcing both unitarity and symmetry, \(\Theta=\Theta^T\).
- **Controlled comparison:** use the same channel realization, system dimensions, architecture, power, noise, and documented fixed beamformer for each pair.

The global optimum over the general unitary feasible set cannot be lower than the global optimum over its reciprocal subset. Local optimizers may still produce a negative *measured* gap because their convergence quality differs. Such cases must be diagnosed and reported, not clipped.

### M2 minimum evidence of progress

1. A common configuration and channel-generation framework.
2. Shared SINR and sum-rate implementation.
3. General and reciprocal scattering-matrix initialization and optimizer pipeline, or clearly documented working subproblems if full solvers are not yet integrated.
4. Numerical unitarity and symmetry checks.
5. At least one small-configuration convergence plot and one matched comparison, with raw outputs and settings.
6. A report that distinguishes completed tests from planned tests.

---

## 2. Build order and ownership

| Step | Owner | Deliverable | Acceptance condition |
|---|---|---|---|
| 1. Freeze conventions | Members 1–4 | `docs/model_decisions.md` | Matrix dimensions, channel orientation, noise units, beamformer rule, and meaning of \(G\) are written down |
| 2. Shared configuration | Member 4 | `configs/default_config.m` | One structure controls all experiment parameters and random seed |
| 3. Channel generation | Member 4 | `src/common/generate_channels.m` | Reproducible Rayleigh and Rician channels; documented path-loss scaling |
| 4. Shared rate evaluator | Member 4 | `src/common/compute_sinr.m`, `compute_sum_rate.m` | Matches hand calculation for a small deterministic test |
| 5. Architecture constraints | Members 2–4 | `src/common/project_architecture.m` or architecture-specific constructors | CWSC/CWGC/CWFC structures are valid and documented |
| 6. General design | Member 3 | `src/nonreciprocal/` | Feasible unitary initialization, objective, update loop, and convergence record |
| 7. Reciprocal design | Member 2 | `src/reciprocal/` | Feasible symmetric-unitary initialization, objective, update loop, and convergence record |
| 8. Matched runner | Member 4 | `experiments/run_matched_comparison.m` | Both methods receive the exact same channel, beamformer, and configuration |
| 9. Initial results | All | `results/raw/`, `results/figures/` | Raw `.mat`/`.csv`, seeds, parameters, and plotted outputs saved |
| 10. Narrative integration | Member 1 | Report sections and figure captions | Literature claims match verified sources; plots are described without overclaiming |

**Important:** the starter code is a scaffold, not a reproduction of either paper's complete optimizer. Algorithm-specific gradients, manifold operations, penalty coefficients, line search, and final enforcement must be transcribed from the original source and independently checked.

---

## 3. Suggested repository structure

```text
bdris-reciprocity-gap/
├── README.md
├── configs/
│   └── default_config.m
├── src/
│   ├── common/
│   │   ├── generate_channels.m
│   │   ├── compute_sinr.m
│   │   ├── compute_sum_rate.m
│   │   ├── make_architecture.m
│   │   └── validate_theta.m
│   ├── reciprocal/
│   │   ├── initialize_theta_reciprocal.m
│   │   ├── objective_reciprocal.m
│   │   └── optimize_reciprocal.m
│   └── nonreciprocal/
│       ├── initialize_theta_unitary.m
│       ├── objective_unitary.m
│       └── optimize_unitary.m
├── experiments/
│   ├── run_single_solver_test.m
│   └── run_matched_comparison.m
├── tests/
│   └── run_all_tests.m
├── results/
│   ├── raw/
│   ├── figures/
│   └── tables/
└── docs/
    ├── model_decisions.md
    ├── sota_comparison.md
    ├── mathematical_formulation.md
    └── validation_log.md
```

Use the team's agreed language and existing repository conventions if they differ; do not create duplicate implementations merely to match this suggested tree.

---

## 4. Function contracts

### 4.1 Configuration

`cfg = default_config()`

Must specify at least:

- \(N\): BS transmit antennas.
- \(K\): single-antenna users.
- \(M\): RIS ports/elements.
- `architecture`: `SC`, `GC`, or `FC`.
- \(G\): number of groups (if GC); document group size and divisibility.
- `Pmax_dBm`, `noise_dBm`.
- `fading`: `Rayleigh` or `Rician`; Rician \(K\)-factor.
- Path-loss reference gain, reference distance, exponent, BS–RIS distance, and RIS–user distances.
- Random seed, optimizer tolerance, maximum iterations, line-search parameters, and restart count.
- Fixed-beamformer method and normalization convention.

Never silently mix dBm and watts. Convert once at the configuration boundary and retain units in variable names.

### 4.2 Channel generation

`[H_TX, H_RX, geom] = generate_channels(cfg)`

Expected dimensions:

\[
H_{\mathrm{TX}}\in\mathbb{C}^{M\times N},\qquad
H_{\mathrm{RX}}\in\mathbb{C}^{K\times M}.
\]

Requirements:

- Same channel arrays are passed to both solvers for a matched trial.
- Path loss is a power gain; channel amplitude is multiplied by its square root.
- Rician LoS array response and geometry are explicitly defined.
- The random seed and all geometry are saved with the trial.

### 4.3 Fixed beamformer

`V = make_fixed_beamformer(H_TX,H_RX,cfg)`

Requirements:

- Compute the beamformer once per channel realization using a documented method.
- Normalize it to meet \(\|V\|_F^2=P_{\max}\) when positive power is specified.
- Do not recompute it separately inside the general and reciprocal solvers.
- Do not allow one solver to receive a beamformer optimized for its own current \(\Theta\) while the other receives a different beamformer.

If the approved project protocol later requires a fixed MMSE beamformer based on a reference effective channel, document that reference channel and regularization explicitly.

### 4.4 Shared rate evaluator

`[R, sinr, user_rates] = compute_sum_rate(H_TX,H_RX,Theta,V,noise_power)`

Use

\[
E=H_{\mathrm{RX}}\Theta H_{\mathrm{TX}},\quad
c_{ki}=e_kv_i,
\]

\[
\gamma_k=\frac{|c_{kk}|^2}{\sum_{i\ne k}|c_{ki}|^2+N_0},\qquad
R=\sum_k\log_2(1+\gamma_k).
\]

The shared evaluator must be used by both optimizers and by final reporting. Do not use different rate implementations for the two methods.

### 4.5 Scattering-matrix validation

`metrics = validate_theta(Theta,cfg)`

At minimum record:

\[
e_U=\|\Theta^H\Theta-I_M\|_F,\qquad
e_S=\|\Theta-\Theta^T\|_F.
\]

For CWGC, also verify off-block entries are zero to numerical tolerance. For CWSC, verify off-diagonal entries are zero and diagonal entries have unit magnitude. Report errors after any final projection/enforcement, not only before it.

### 4.6 Optimizer interfaces

Recommended interfaces:

```matlab
[Theta, history, info] = optimize_unitary(H_TX,H_RX,V,cfg,Theta0)
[Theta, history, info] = optimize_reciprocal(H_TX,H_RX,V,cfg,Theta0)
```

`history` should include objective/rate per iteration, gradient norm if available, step size, and feasibility errors. `info` should include convergence status, stop reason, iteration count, runtime, and restart identifier.

The two solvers may use different constraint-aware algorithms, but they must share the same channel, fixed \(V\), power/noise, architecture, and final rate evaluator. A measured gap can reflect both feasible-set restriction and optimizer quality; disclose this limitation.

---

## 5. Validation tests T1–T14

These are acceptance tests to run. Mark each as `PASS`, `FAIL`, or `NOT RUN` in `docs/validation_log.md`, with a date, commit hash, parameters, and saved output. A test is not passed merely because code exists.

| ID | Test | Procedure | Expected outcome |
|---|---|---|---|
| T1 | Configuration/units | Convert \(20\) dBm and \(-80\) dBm using \(P[\mathrm W]=10^{(P[\mathrm{dBm}]-30)/10}\) | \(20\) dBm \(=0.1\) W; \(-80\) dBm \(=10^{-11}\) W |
| T2 | Channel dimensions | Generate channels for \(N=2,K=2,M=32\) | `H_TX` is \(32\times2\), `H_RX` is \(2\times32\) |
| T3 | Seed reproducibility | Generate twice with same seed and config | Channels match to numerical precision |
| T4 | Fading/path-loss sanity | Run Rayleigh and Rician samples and inspect empirical power | Finite values; average scale is consistent with configured path loss and fading normalization |
| T5 | Beamformer power | Compute fixed \(V\) | \(\|V\|_F^2=P_{\max}\) within tolerance, unless a documented zero-power case |
| T6 | Deterministic rate check | Use a small hand-constructed effective channel and compare manual SINR/rate with shared function | Same SINR and sum-rate to numerical tolerance |
| T7 | General unitary feasibility | Validate general \(\Theta\) | \(\|\Theta^H\Theta-I\|_F\) below stated tolerance |
| T8 | Reciprocal feasibility | Validate reciprocal \(\Theta\) after final enforcement | Both unitarity and \(\|\Theta-\Theta^T\|_F\) below stated tolerance |
| T9 | Architecture structure | Check SC diagonal, GC block diagonal, FC full dimensions | Forbidden entries are zero within tolerance; block sizes match config |
| T10 | CWSC zero-gap control | Run both formulations for same SC instance using equivalent feasible solutions/solver settings | Reciprocity constraint adds no feasible-set restriction; gap should be zero within numerical/optimization tolerance |
| T11 | Gradient finite difference | Compare analytical directional derivative with centered finite difference on a small problem | Relative derivative error is small and decreases as step size is adjusted, subject to conditioning |
| T12 | Monotonicity/convergence | Inspect accepted optimization updates | Objective/rate does not decrease beyond tolerance when algorithm's line search promises ascent; otherwise document exact objective being monitored |
| T13 | Matched trial | Run both methods with saved common config, channels, \(V\), and seed | Metadata confirms identical problem instance and same shared rate evaluator |
| T14 | Reproducible output | Run experiment twice with same seed and deterministic settings | Rates and saved metadata reproduce within numerical tolerance; any nondeterminism is documented |

### Suggested tolerance policy

Choose tolerances based on numerical precision and algorithm implementation; do not invent one after seeing the result. As an initial engineering convention for double precision, feasibility errors around \(10^{-6}\) may be a reasonable target for small examples, but tolerances must be adjusted and justified for larger matrices and the specific algorithm.

### Validation log template

```text
Test ID:
Date / commit:
MATLAB version / toolbox:
Configuration and random seed:
Command or script:
Observed output:
Expected output:
PASS / FAIL / NOT RUN:
If failed, diagnosis and next action:
Saved evidence path:
```

---

## 6. Fair-comparison protocol

1. **Pair by realization.** For trial \(\ell\), generate \(H_{\mathrm{TX}}^{(\ell)}\) and \(H_{\mathrm{RX}}^{(\ell)}\) once and reuse them for both methods.
2. **Freeze all shared settings.** \(N,K,M,G,P_{\max},N_0\), fading, path loss, geometry, and architecture must be identical.
3. **Freeze the beamformer.** Compute \(V^{(\ell)}\) once, normalize once, and pass the exact matrix to both optimizers.
4. **Use the same rate evaluator.** Evaluate final feasible matrices with the same code.
5. **Record algorithmic effort.** Document restart count, stopping criteria, initialization policy, maximum iterations, and runtime.
6. **Enforce before measuring.** Evaluate rate using the final feasible \(\Theta\), after any symmetry/unitarity enforcement.
7. **Do not hide negative gaps.** Report them and inspect convergence, initialization, and feasibility.
8. **Keep reproduction and project comparisons separate.** When reproducing a source paper, use that paper's own parameters. For the new controlled experiment, use one common parameter set.
9. **Save raw data.** Each result row should include seed, architecture, \(N,K,M,G\), power, fading, \(R_U\), \(R_S\), absolute/relative gaps, feasibility errors, runtime, iterations, and stop reason.
10. **Do not overstate novelty.** The direct-comparison paper identified in the SOTA notes must be inspected in full before the gap is finalized.

---

## 7. Initial M2 configuration

The four-member plan specifies this starting point for Members 2 and 3:

- \(N=2\) transmit antennas.
- \(K=2\) single-antenna users.
- \(M=32\) RIS elements/ports.
- \(P_{\max}=20\) dBm.
- Rayleigh fading.
- Begin with one architecture, then extend to SC/GC/FC.
- Use the same channel realization and system parameters for both methods.

The proposal's broader ranges (\(N=K\in\{2,\ldots,6\}\), \(M\in[8,64]\), \(G\in\{1,2,4,M\}\), and \(P_{\max}\in[0,20]\) dBm) are suitable for the later sweep after the small case works. Confirm whether \(G\) denotes group count or group size in the final implementation and use one convention consistently.

---

## 8. Schedule and completion checklist

### Suggested five-work-session schedule

| Session | Main work | Evidence to save |
|---|---|---|
| 1 | Lock model conventions; inspect paper equations/algorithms; create config and validation log | `model_decisions.md`, config, reference notes |
| 2 | Implement and test channels, fixed beamformer, SINR and sum-rate | T1–T6 logs and deterministic rate example |
| 3 | Implement reciprocal and general initializers/updates on the small case | T7–T9 logs, convergence histories |
| 4 | Run matched comparison; debug failed tests; save raw outputs and first plots | T10–T14 logs, `.mat`/`.csv`, figures |
| 5 | Integrate report, code README, contribution commits and video | final M2 PDF, repository, 20-minute walkthrough |

This is a work-session sequence, not a guarantee that all algorithms can be implemented in five sessions. If a solver is not working, report the blocker and preserve verified subcomponents rather than relabeling incomplete work as complete.

### Checklist

- [ ] SOTA table and gap updated against the full competing paper.
- [ ] Model conventions documented.
- [ ] Shared configuration and channel generation run.
- [ ] Shared SINR/sum-rate evaluator passes deterministic check.
- [ ] Fixed beamformer rule is documented and identical for both methods.
- [ ] General solver's feasible initialization and objective run.
- [ ] Reciprocal solver's feasible initialization and objective run.
- [ ] Unitarity and symmetry errors are logged.
- [ ] At least one convergence experiment is saved.
- [ ] At least one matched comparison is saved.
- [ ] Raw data, configuration, random seeds, and figures are committed.
- [ ] Each member has descriptive incremental commits and can explain their own work.
- [ ] Report separates completed results from future M3 work.

---

## 9. Video walkthrough outline (20 minutes)

| Time | Presenter/topic | What to show |
|---|---|---|
| 0:00–2:00 | Member 1: motivation, SOTA and qualified research gap | Project question; comparison matrix; why the direct competing paper matters |
| 2:00–5:00 | Member 1: mathematical formulation | Signal model, SINR, sum-rate, unitary and symmetry constraints, gap KPIs |
| 5:00–9:00 | Member 4: common simulation environment | Configuration, channel dimensions, fading/path-loss, fixed beamformer, shared rate evaluator |
| 9:00–13:00 | Member 2: reciprocal method | Symmetric-unitary initialization, optimizer steps, feasibility checks and actual convergence evidence |
| 13:00–16:00 | Member 3: general baseline | General unitary design, architecture support, convergence evidence and baseline-reproduction status |
| 16:00–18:30 | Member 4: matched result | Same seed/channels/config, initial comparison plot, saved raw results and any negative-gap diagnosis |
| 18:30–20:00 | All: limitations and next steps | What passed, what remains incomplete, M3 sweeps, complexity analysis, and reproducibility |

Do not narrate expected plots as if they were measured results. If a component has not run, say so and show the current blocker or test plan.

---

## 10. M2/M3 boundary and risk handling

**M2:** common environment, small-case optimizer progress, sanity checks, initial convergence and matched comparison, SOTA positioning, report/video/repository.

**M3:** systematic sweeps over architecture, \(M\), \(G\), power and fading; statistical gap characterization; full complexity/performance Pareto analysis; additional robustness and convergence experiments.

**Risk: full solver too time-consuming.** First isolate and validate the scattering-matrix subproblem with a fixed beamformer if this is permitted by the team's agreed scope. Clearly label this as a simplification and compare it to the relevant source-paper setup. Do not claim to reproduce the full joint precoder/RIS algorithm if the precoder was fixed.

**Risk: unstable manifold optimization.** Reduce the test size, verify gradients and retractions, inspect line-search acceptance and constraint errors, and document the stability boundary. Do not conceal failed runs.

**Risk: novelty overlap.** Read the competing paper in full and narrow the research claim to what remains untested under the agreed protocol.

---

## 11. Git commit guidance

Use small, descriptive commits that reflect real work, for example:

- `docs: add SOTA comparison and qualified gap`
- `docs: derive system model and rate equations`
- `common: add reproducible channel generator`
- `common: validate SINR and sum-rate evaluator`
- `reciprocal: add symmetric-unitary initialization`
- `reciprocal: validate symmetry and unitarity`
- `baseline: add unitary initialization`
- `experiments: add matched small-case runner`
- `tests: record T1-T14 validation outcomes`
- `results: add initial convergence plot and raw data`

Only use a commit message after the corresponding work has actually been performed. Individual contribution evidence must reflect actual contributions.
