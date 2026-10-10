# Reciprocal BD-RIS Optimization Pipeline

## Pipeline Explanation
1. **Initialize**: The architecture specifies the grouping $R_g$ (e.g. SC has $R_g=1$, GC2 has $R_g=2$, FC has $R_g=32$). For the Monte Carlo and overlay plots, we use a shared diagonal SC starting point, and compute $V$ from it, so the comparison is perfectly fair across architectures. For the single-draw convergence run and error report, we use a random symmetric unitary FC start.
2. **MMSE Beamformer**: A fixed minimum mean square error (MMSE) precoder `V` is computed based on the initial effective channel to provide an initial transmission baseline.
3. **Auxiliary Variables**: The FP auxiliary variables `tau` (user SINRs) and `y` are computed given the current beamformer and scattering matrix.
4. **Objective & Gradient**: The regularized FP sum-rate surrogate objective is exactly evaluated: $F_{FP}(\Theta) - \nu\|\Theta - \Theta^T\|_F^2$. A Euclidean gradient is computed that includes both the sum-rate derivatives and the penalty term derivative.
5. **Tangent Projection**: The Euclidean gradient is projected onto the tangent space of the Stiefel (unitary) manifold, yielding the Riemannian gradient. Note that the gradient is only compared and evaluated inside the diagonal blocks, because off-block entries are zero by construction for group-connected and single-connected architectures.
6. **Conjugate Gradient Ascent (CGA)**: The search direction is updated using the Polak-Ribière rule: $\beta = \max(0, \frac{\langle r_{new}, r_{new} - r_{transported} \rangle}{\langle r_{transported}, r_{transported} \rangle})$. If the direction ceases to be an ascent direction (i.e., $\langle r, \Xi \rangle \le 0$), it resets to the steepest ascent direction ($r$).
7. **Armijo Line Search**: A backtracking line search guarantees sufficient ascent on the FP surrogate objective with $\tau$, $y$, and $V$ frozen. The Armijo sufficient increase coefficient is set to $c=2e-11$ (small enough to accept valid ascent steps in a flat non-convex landscape), and the first trial step is scaled by $1 / \|\Xi\|_F$ to prevent massive initial overshoots. The 200-step limit is strictly monitored and logged (0 of 50 draws reached the limit in our tests).
8. **Iterate**: After an accepted step, the loop refreshes the auxiliary variables ($\tau$, $y$), computes new gradients, and applies retraction steps until the absolute difference in the surrogate objective between iterations falls below $\epsilon = 10^{-8}$.
9. **Final Projection**: The objective contains a penalty $-\nu\|\Theta - \Theta^T\|_F^2$ which only softly pushes $\Theta$ toward symmetry. At termination, `enforce_symmetry_unitarity.m` applies a strict projection. First, it explicitly symmetrizes the final block components via $\Theta_{sym} = 0.5(\Theta_g + \Theta_g^T)$. Then it applies an SVD ($U \Sigma V^H$) and sets $\Theta_g = U V^H$. This produces a matrix that is precisely unitary. Furthermore, because $U V^H$ is the unitary polar factor of $\Theta_{sym}$, and $\Theta_{sym}$ is symmetric, its polar factor is also symmetric. Thus, the resulting matrix strictly satisfies both unitarity and symmetry.

## System and Algorithm Configuration
- **System Parameters (`config.m`)**: `N=2` (TX antennas), `K=2` (Users), `R=16/32` (Elements), `Pmax=20 dBm`. Path loss uses `C0=-30 dB` at `1m`, `rho=2.2`. Distance to RIS is `50m`, RIS to users is `2.5m`. Noise `N0=-80 dBm`. Fading is Rayleigh (`Kf=-Inf`) or Rician (`Kf=2 dB`).
- **Algorithm Settings (`run_cga_optimizer.m`, `armijo_search.m`)**: 
  - Penalty weight `nu = 1.0`
  - Convergence tolerance `eps = 1e-8` (absolute objective difference)
  - Maximum iterations = `8000`
  - Armijo maximum steps = `200`
  - Armijo contraction factor = `0.75`
  - Armijo sufficient increase coefficient = `2e-11`

### Deviations from Default / Custom Design Choices
- **FP Variables ($\tau$ and $y$)**: The algorithm treats the auxiliary variable $\tau$ exactly as the SINR $\gamma_k$. The variable $y_k$ is computed as $y_k = \frac{\mathbf{e}_k \mathbf{v}_k}{\sum_i |\mathbf{e}_k \mathbf{v}_i|^2 + \sigma^2}$, directly matching the FP transformation conditions for a fixed precoder.
- **First-Step Scaling**: The initial Armijo step `alpha_0` is scaled by `1/||Xi||_F` to prevent extreme overshoot on the first line-search evaluation.
- **QR Retraction Phase Fix**: QR decomposition is only unique up to a diagonal phase matrix. Forcing `diag(R)` to be strictly positive real (`diag(R)./abs(diag(R))`) makes the retraction well-defined. This ensures $R(\Theta, 0) = \Theta$ and makes the retraction a smooth, contiguous map back to the Stiefel manifold.
- **Vector Transport & PR Update**: The conjugate direction `Xi` is transported to the new tangent space simply by re-projecting it onto the updated tangent space after a step. Furthermore, in the Polak-Ribière formula, the denominator specifically uses the transported old gradient (`r_transported`) rather than the un-transported old gradient. This is a deliberate convention choice for manifold stability.
- **Gradient Penalty Convention**: The Wirtinger derivative definition in Eq. (36) is $2 \cdot \frac{\partial F}{\partial \Theta^*}$. Because the penalty is $-\nu\|\Theta - \Theta^T\|_F^2$, its Wirtinger derivative with respect to $\Theta^*$ incorporates a factor of 2, resulting in the penalty gradient component $-4\nu(\Theta - \Theta^T)$.

## Results and Verification
- **Finite-Difference Accuracy**: The analytic Euclidean gradient matches central finite differences to a relative error of $\sim 4.5 \times 10^{-9}$ across block sizes $R_g \in \{1, 4, 32\}$ (reproducible via `test_gradient_fd.m`). Note that since $1 \times 1$ blocks are trivially symmetric, the penalty gradient evaluation is only strictly meaningful for the non-scalar blocks $R_g \in \{4, 32\}$.
- **Manifold Validation**: The skew-Hermitian tangent condition and unitary retraction constraints are both maintained to machine precision ($\sim 10^{-15}$).
- **Explicit Unitarity/Symmetry Report (FC R=32 Single Draw)**:
  - Initial (random): $||\Theta - \Theta^T||_F = 0$, $||\Theta^H \Theta - I||_F = 5.1 \times 10^{-15}$
  - Optimized (before projection): $||\Theta - \Theta^T||_F = 8.2 \times 10^{-3}$, $||\Theta^H \Theta - I||_F = 3.4 \times 10^{-15}$
  - Final (after projection): $||\Theta - \Theta^T||_F = 1.9 \times 10^{-14}$, $||\Theta^H \Theta - I||_F = 6.1 \times 10^{-15}$
- **Projection Cost**: The final symmetrization/SVD step strictly enforces symmetry at a negligible sum-rate penalty of $\sim 1.1 \times 10^{-3} \text{ bps/Hz}$ (computed over the Monte Carlo ensemble), confirming that $\nu = 1.0$ is sufficient.
- **Iteration Spread (FC R=32)**: The 50-draw Monte Carlo convergence yields a mean of 594 iterations (std = 789), with a range spanning 59 to 4395 iterations. 
- **Known Discrepancy (10.32 vs ~14 bps/Hz)**: The FC $R=32$ architecture yields an average converged sum-rate of $10.32 \text{ bps/Hz}$ over 50 draws, with individual draws reaching up to $14.8 \text{ bps/Hz}$ (or an average of $10.87 \text{ bps/Hz}$ if the precoder $V$ is optimally recomputed at termination). (Note: The 20-draw overlay figure reports an average of 9.62 bps/Hz due to a different draw set). This falls short of the $\sim 14 \text{ bps/Hz}$ target in Fig. 3 of the paper. We do not tune parameters to close this gap. We hypothesize this gap is due to differing unstated simulation setup assumptions:
  1. **Fixed vs Joint Beamformer**: We hold $V$ fixed at its initial MMSE value throughout the CGA iterations (a valid fallback in the proposal's Section 11), while the paper might optimize $V$ and $\Theta$ jointly. Recomputing $V$ at the end already recovers $+0.55 \text{ bps/Hz}$.
  2. **Power Definition**: $P_{max}$ might define per-user power rather than total base station power.

## How to Run
- `plot_convergence.m`: Generates the average sum-rate convergence overlay (20 draws) for SC, GC(2), GC(4), and FC architectures.
- `run_monte_carlo.m`: Evaluates 50 independent channel draws for FC to generate statistically reliable performance gaps and statistics.
- `test_stopping_rule.m`: Prints iterations, rate, alpha and gnorm at three tolerances for inspection.
- `explicit_error_report.m`: Validates numerical constraints of the algorithm before and after final exact projection.
- `test_gradient_fd.m`: A reproducible finite-difference check of the analytic gradient for $R_g \in \{1, 4, 32\}$.

## File Organization
- `initialize_theta.m`, `generate_channels.m`, `mmse_beamformer.m`: Initialization components.
- `objective.m`, `gradient_theta.m`, `fp_variables.m`: FP transform evaluation and analytic gradients.
- `tangent_projection.m`, `retraction.m`, `enforce_symmetry_unitarity.m`: Riemannian geometry operators on the Stiefel manifold.
- `run_cga_optimizer.m`, `armijo_search.m`: The primary loop and Armijo backtracking search.
- `test_blocks.m`, `test_milestone3.m`, `run_convergence_experiment.m`: Diagnostics and convergence plotting.
