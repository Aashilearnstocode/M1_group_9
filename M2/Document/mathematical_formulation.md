# Mathematical Formulation — Group 9 BD-RIS Reciprocity-Gap Project 

Status legend used throughout:
**[P1]** = taken from Li–Shen–Clerckx (checked against our PDF copy) ·
**[P2]** = taken from Fidanovski et al., arXiv 2509.20246 v4 (read through a web-fetch summary — re-check against the PDF) ·
**[DERIVED]** = derived by us here (check by hand / numerically) ·
**[TO VERIFY]** = ambiguity that must be settled before experiments.

---
## 0. Notation (follows [P2]; mapping to [P1] and to M1)

| Symbol | Meaning | Size | [P1] name | M1 name |
|---|---|---|---|---|
| N | BS antennas | – | N | N |
| K | single-antenna users (K = N in our sweeps) | – | K | K |
| M | RIS elements/cells (= R in [P2]) | – | M | M |
| G | number of groups | – | G | G |
| $\bar M = M/G$ | group size (= $R_G$ in [P2]) | – | $\bar M$ | – |
| $\mathbf H_{TX}$ | BS→RIS channel | M×N | $\mathbf G$ | $\mathbf G$ |
| $\mathbf H_{RX}$ | RIS→users; row k is $\mathbf h_k^T$ | K×M | $\mathbf h_k^H$ rows | $\mathbf h_{r,k}^H$ |
| $\boldsymbol\Theta$ | RIS scattering matrix | M×M | $\boldsymbol\Phi_r$ (reflective mode) | $\boldsymbol\Theta$ |
| $\mathbf V=[\mathbf v_1..\mathbf v_K]$ | transmit beamformer | N×K | $\mathbf W$ | $\mathbf w_k$ |
| $P_{\max}$ | total transmit power | W | $P$ | $P_{\max}$ |
| $N_0$ | noise power per user | W | $\sigma_k^2$ | $N_0$ |

**Recommendation:** use $\mathbf H_{TX}$ for the BS→RIS channel so that the letter G only ever means "number of groups".

Note: [P1] writes $\mathbf h_k^H$ (Hermitian) and [P2] writes $\mathbf h_k^T$ (transpose); they differ only by a conjugation of the channel, which does not change statistics for i.i.d. fading. Pick one and state it.

---
## 1. System model (MU-MISO downlink, blocked direct link)

**Why blocked direct link.** Both source papers assume it ([P1] Sec. IV-A, footnote 6; [P2] Fig. 1). M1 wrote $\mathbf h_{d,k}$ into the received signal; we recommend dropping it for reproduction (it can be an extension). **[TO VERIFY by team]**

**Transmit signal.** $\mathbf x=\mathbf V\mathbf s$, with $\mathbf s\in\mathbb C^K$, $\mathbb E[\mathbf s\mathbf s^H]=\mathbf I_K$ (one unit-power stream per user), and power constraint
$$\|\mathbf V\|_F^2\le P_{\max}. \tag{1}$$ [P1 (13c), P2 (12b)]

**Received signal** at user k [P2 (1), P1 (11)]:
$$ r_k=\mathbf h_k^T\,\boldsymbol\Theta\,\mathbf H_{TX}\,\mathbf V\mathbf s+n_k,\qquad n_k\sim\mathcal{CN}(0,N_0). \tag{2}$$
Stack the users: $\mathbf r=\mathbf H_{RX}\boldsymbol\Theta\mathbf H_{TX}\mathbf V\mathbf s+\mathbf n$.

Reading (2) left to right: $\mathbf V\mathbf s$ shapes N antenna signals; $\mathbf H_{TX}$ carries them to the M RIS ports; $\boldsymbol\Theta$ mixes the M port signals (this is the "beyond diagonal" part — a diagonal $\boldsymbol\Theta$ would only phase-rotate each port independently); $\mathbf h_k^T$ carries the result to user k.

**Effective (equivalent) channel** [P2 (16)]:
$$\mathbf E=\mathbf H_{RX}\boldsymbol\Theta\mathbf H_{TX}\in\mathbb C^{K\times N},\qquad \mathbf e_k=\mathbf h_k^T\boldsymbol\Theta\mathbf H_{TX}\ (\text{k-th row}).$$
Define the scalar $c_{ki}\triangleq\mathbf e_k\mathbf v_i$ = signal from stream i seen at user k.

**SINR and sum-rate** [P2 (11),(17),(20); P1 (12),(13a)]:
$$\gamma_k=\frac{|c_{kk}|^2}{\sum_{i\ne k}|c_{ki}|^2+N_0},\qquad R=\sum_{k=1}^K\log_2(1+\gamma_k)\ \ \text{[bit/s/Hz]}. \tag{3}$$
Numerator = desired stream; denominator = inter-user interference + noise. $\log_2(1+\mathrm{SINR})$ is the Shannon spectral efficiency of a Gaussian channel with treat-interference-as-noise.

### 1.1 Channel generation (what we must code) 
Large-scale power gain [P1 Sec. V, P2]: $\mathrm{PL}(d)=\zeta_0(d/d_0)^{-\varepsilon}$, $\zeta_0=-30$ dB, $d_0=1$ m, $\varepsilon=2.2$.

Small-scale (our reading, **[TO VERIFY]**): 
$$\mathbf H=\sqrt{\mathrm{PL}}\Big(\sqrt{\tfrac{\kappa}{\kappa+1}}\bar{\mathbf H}_{LoS}+\sqrt{\tfrac{1}{\kappa+1}}\tilde{\mathbf H}_{NLoS}\Big),\quad \tilde H_{ij}\sim\mathcal{CN}(0,1)$$
with $\kappa=0$ for Rayleigh. Rician factor is **5 dB in [P1]** (κ_BI = κ_IU) and **2 dB in [P2]**.

**Ambiguities that the papers do not settle in the text we have** (do *not* assume silently):
1. Whether PL multiplies the channel power (so $\sqrt{\mathrm{PL}}$ on the amplitude) — we assume yes.
2. The LoS array-response model (ULA/UPA, spacing, angles) — not in the text we read.
3. Units: P_max in dBm → watts ($10^{(\cdot)/10-3}$); N₀ = −80 dBm → $10^{-11}$ W. Keep everything in watts.
4. Whether both papers use the same user placement (users "randomly located" at 2.5 m from RIS in [P1]).

**Hand check of the link budget [DERIVED, verified numerically]:**
BS–RIS: $-30-22\log_{10}50=-67.38$ dB; RIS–user: $-30-22\log_{10}2.5=-38.75$ dB; cascade $=-106.13$ dB.
At $P_{\max}=20$ dBm: $20-106.13=-86.13$ dBm vs noise $-80$ dBm ⇒ **−6.1 dB per-link SNR before any array/RIS/beamforming gain**. At 0 dBm it is −26.1 dB. So sum-rates are dominated by array gain from M and N — a useful sanity check on every simulated curve (rates should not be near zero at 20 dBm with M = 32).

---
## 2. Scattering matrix: unitary and symmetric constraints

**Lossless ⇒ unitary [P1 Sec. II, [P2 (6)].** A passive lossless multiport reflects all incident power: $\|\boldsymbol\Theta\mathbf a\|^2=\|\mathbf a\|^2\ \forall\mathbf a\Rightarrow\boldsymbol\Theta^H\boldsymbol\Theta=\mathbf I$ [DERIVED]. (If lossy but passive: only $\boldsymbol\Theta^H\boldsymbol\Theta\preceq\mathbf I$.)

**Reciprocal ⇒ symmetric.** For a network of reciprocal components (no circulators/magnets/active non-reciprocal parts), the scattering matrix satisfies $\boldsymbol\Theta=\boldsymbol\Theta^T$ (Lorentz reciprocity; [P2 (5)]; physical argument is in Shen–Clerckx–Murch, our M1 ref [3], which is [7] inside [P1]). **We have not read [3] ourselves — cite carefully.** [P1] itself never mentions reciprocity (0 occurrences of "reciprocal"/"symmetric" in our PDF copy).

**Feasible sets** (single block, size n):
- $\mathcal U(n)=\{\boldsymbol\Theta:\boldsymbol\Theta^H\boldsymbol\Theta=\mathbf I_n\}$ — unitary group, **real dimension $n^2$**. (Square case of the complex Stiefel manifold.)
- $\mathcal U_S(n)=\{\boldsymbol\Theta\in\mathcal U(n):\boldsymbol\Theta=\boldsymbol\Theta^T\}$ — symmetric unitary, **real dimension $n(n+1)/2$** [DERIVED; consistent with Santamaria et al., arXiv 2604.22469, who state $M(M+1)/2$ circuits for fully-connected].

**Takagi factorisation [DERIVED/standard]:** every $\boldsymbol\Theta\in\mathcal U_S(n)$ can be written $\boldsymbol\Theta=\mathbf Q\mathbf Q^T$ with $\mathbf Q\in\mathcal U(n)$. This gives a way to *generate* feasible random starting points and to test feasibility.

**Tangent spaces [DERIVED by differentiating the constraints]:**
- At $\boldsymbol\Theta\in\mathcal U(n)$: $T=\{\boldsymbol\Xi:\ \boldsymbol\Theta^H\boldsymbol\Xi+\boldsymbol\Xi^H\boldsymbol\Theta=\mathbf 0\}$ (i.e. $\boldsymbol\Theta^H\boldsymbol\Xi$ skew-Hermitian). Orthogonal projection of a Euclidean gradient $\mathbf Z$:
$$\mathrm{Proj}_{\boldsymbol\Theta}(\mathbf Z)=\mathbf Z-\boldsymbol\Theta\,\tfrac12(\boldsymbol\Theta^H\mathbf Z+\mathbf Z^H\boldsymbol\Theta). \tag{4}$$ This is [P2 (42)–(43)] exactly.
- At $\boldsymbol\Theta\in\mathcal U_S(n)$: additionally $\boldsymbol\Xi=\boldsymbol\Xi^T$. The exact projection/retraction/geodesics are worked out in Santamaria et al. (2604.22469). **[P2] does not use this manifold**; it uses $\mathcal U(n)$ plus a penalty (Sec. 4).

> **Flag on [P1]:** as printed, [P1] (28) gives the tangent space as $\Re\{\boldsymbol\Phi_g^H\mathbf T_g\}=\mathbf 0$ and (30) projects with $\mathrm{chdiag}(\cdot)$ (diagonal part only). That does not look like the Stiefel-manifold projection (4); it resembles a column-norm (oblique-type) projection. The retraction (33) is polar-type and only lands exactly on the manifold for proper tangent directions. **We should implement the standard Stiefel/unitary projection (4) (as in Manopt/Boumal) and say so explicitly in the report; do not copy (30) blindly.** [TO VERIFY with Member 3/Prof.]

### 2.1 Architectures as structural masks
| Architecture | Structure of Θ | Groups | Unconstrained set | Reciprocal set | Real dof (U / U_S) |
|---|---|---|---|---|---|
| CWSC | diagonal | $G=M$ | $\{|θ_{mm}|=1\}$ | same set (diagonal ⇒ symmetric) | $M$ / $M$ |
| CWGC | $\mathrm{blkdiag}(\Theta_1..\Theta_G)$, size $\bar M=M/G$ | $1<G<M$ | each $\Theta_g\in\mathcal U(\bar M)$ | each $\Theta_g\in\mathcal U_S(\bar M)$ | $G\bar M^2$ / $G\bar M(\bar M+1)/2$ |
| CWFC | full | $G=1$ | $\mathcal U(M)$ | $\mathcal U_S(M)$ | $M^2$ / $M(M+1)/2$ |

([P2 (7)–(9)]; [P1 (6),(8),(10)] with $\boldsymbol\Phi_t=0$.) **Group-count convention:** [P2] says $G=1$ is FC and $G=R$ is SC; so M1's list $G\in\{1,2,4,M\}$ is a list of *group counts*. [P1] uses $G=8$ groups at $M=32$ (group size 4); [P2] writes GC(2), GC(4) meaning group *size* 2, 4. **Define it once in the report.**

**Consequence [DERIVED] — CWSC cannot show a reciprocity gap**, because a diagonal unitary matrix is already symmetric: $\mathcal U_S=\mathcal U$ for $\bar M=1$. Hence $\Delta_{CWSC}=0$ *exactly* and is a free correctness test. Also for $\bar M=2$: 4 vs 3 real dof (smallest non-trivial case); for large $\bar M$ the constraint removes about half of the dof.

---
## 3. The optimisation problem (our scope: fixed beamformer)

Given $\mathbf H_{TX},\mathbf H_{RX}$ and a **fixed** $\mathbf V$ with $\|\mathbf V\|_F^2=P_{\max}$:
$$\textbf{(P}_{\mathcal A}\textbf{)}\quad R^\star_{\mathcal A}=\max_{\boldsymbol\Theta\in\mathcal A}\ \sum_k\log_2(1+\gamma_k(\boldsymbol\Theta;\mathbf V)),\qquad \mathcal A\in\{\mathcal A^{U}_{arch},\ \mathcal A^{S}_{arch}\}. \tag{5}$$
where $\mathcal A^S_{arch}\subseteq\mathcal A^U_{arch}$ for each architecture (Section 2.1). This is [P2] problem P2/(18) (MMSE-fixed beamformer, Θ only) and the Θ-block of [P1] Algorithm 1.

**Fixed beamformer.** [P2] uses an MMSE beamformer. **[TO VERIFY]:** *for which Θ* is the MMSE beamformer computed (random initial Θ? identity?) — [P1] eq. (42) computes it from the initial diagonal Φ then normalises $\mathbf V\leftarrow\sqrt{P_{\max}}\,\tilde{\mathbf V}/\|\tilde{\mathbf V}\|_F$. For our comparison the rule must be: **one $\mathbf V$ per channel realisation, computed once from the shared initial Θ₀ (or any rule independent of the design), and given identically to both solvers.**

---
## 4. Fractional programming (FP) — how a sum of logs becomes tractable

**Step 1 — Lagrangian dual transform (LDT)** [P1 (14); P2 (22)]. For $\gamma>0$:
$$\log(1+\gamma)=\max_{\iota\ge0}\Big\{\log(1+\iota)-\iota+\frac{(1+\iota)\gamma}{1+\gamma}\Big\}.$$
*Proof [DERIVED]:* set derivative in ι to zero: $\frac1{1+\iota}-1+\frac{\gamma}{1+\gamma}=0\Rightarrow\iota=\gamma$, and substituting gives $\log(1+\gamma)$. Note $\frac{\gamma}{1+\gamma}=\frac{|c_{kk}|^2}{\sum_i|c_{ki}|^2+N_0}$ (sum over *all* i incl. k). The log has left the fraction, and ι is auxiliary (called $\iota_k$ in [P1], $\tau_k$ in [P2]; it is real — [P2] types it complex, a harmless typo).

**Step 2 — Quadratic transform (QT)** [P1 (15); P2 (23)–(24)]. For $B>0$:
$$\frac{|A|^2}{B}=\max_{y\in\mathbb C}\ \big\{2\Re(y^*A)-|y|^2B\big\},\quad y^\star=A/B.$$
*Proof [DERIVED]:* complete the square, $-B|y-A/B|^2+|A|^2/B$. Apply with $A=c_{kk}$, $B=\sum_i|c_{ki}|^2+N_0$:
$$\hat\eta_k=\log_2(1+\iota_k)-\frac{\iota_k}{\ln2}+\frac{1+\iota_k}{\ln2}\Big[2\Re\{y_k^*c_{kk}\}-|y_k|^2\Big(\sum_i|c_{ki}|^2+N_0\Big)\Big]. \tag{6}$$
**Block-coordinate updates (closed form):** $\iota_k=\gamma_k$, $\;y_k=c_{kk}/(\sum_i|c_{ki}|^2+N_0)$. Alternate: (ι,y) ↔ Θ. Each step cannot decrease the objective ⇒ monotone, converges to a stationary point (not a proven global optimum — [P1] states convergence "cannot be strictly proved").

**Why it matters for fairness:** both papers use the *same* LDT+QT; M1 attributes this to Shen–Yu [7]. Using one FP outer loop for both manifolds means any rate difference comes from the constraint set, not the outer loop.

Note M1 says [P2] uses "Dinkelbach"; per the paper it is **LDT + QT**, not Dinkelbach. Fix in M1.

---
## 5. Θ-subproblem, gradients and the manifold algorithms

For fixed $(\iota,y,\mathbf V)$, define $\mathbf g_i=\mathbf H_{TX}\mathbf v_i\in\mathbb C^M$, so $c_{ki}=\mathbf h_k^T\boldsymbol\Theta\mathbf g_i$. Up to constants (ignoring $1/\ln 2$ scaling), maximise
$$f(\boldsymbol\Theta)=\sum_k(1+\iota_k)\Big[2\Re\{y_k^*c_{kk}\}-|y_k|^2\sum_i|c_{ki}|^2\Big]. \tag{7}$$

**Gradient convention.** For real $f$, define $\nabla f$ by $df=\Re\,\mathrm{tr}(\nabla f^H d\boldsymbol\Theta)$.
- $d|c_{ki}|^2=2\Re(c_{ki}^*\,\mathbf h_k^Td\boldsymbol\Theta\,\mathbf g_i)\Rightarrow\nabla|c_{ki}|^2=2c_{ki}\,\mathbf h_k^*\mathbf g_i^H$ [matches [P2] Lemma 2].
- $\nabla\,2\Re\{y^*c_{kk}\}=2y\,\mathbf h_k^*\mathbf g_k^H$.
$$\boxed{\nabla f=\sum_k(1+\iota_k)\Big[2y_k\mathbf h_k^*\mathbf g_k^H-2|y_k|^2\sum_ic_{ki}\mathbf h_k^*\mathbf g_i^H\Big]}\tag{8}$$ [DERIVED; agrees with [P2] (36) after conjugation bookkeeping.] For group g, use only the g-th block rows/cols: $\mathbf h_k^{(g)}$ and $\mathbf g_i^{(g)}$ — **but the cross-group terms matter**: $c_{ki}=\sum_g\mathbf h_k^{(g)T}\boldsymbol\Theta_g\mathbf g_i^{(g)}$, so $c_{ki}$ must be recomputed with all blocks (this is the "groups are coupled through the shared equivalent channel" point). [P2] leaves $\mathbf W^{(g)}$ undefined in the text we saw — it is presumably $\mathbf H_{TX}^{(g)}$ [TO VERIFY].

**Penalty for symmetry** [P2 (21a),(25a),(31)]: subtract $\nu\|\boldsymbol\Theta-\boldsymbol\Theta^T\|_F^2$, gradient contribution $-4\nu(\boldsymbol\Theta-\boldsymbol\Theta^T)$ [DERIVED: with $A=\Theta-\Theta^T$, $d\|A\|^2=4\Re\,\mathrm{tr}(A^Hd\Theta)$ ✓.]

### 5.1 Unconstrained design (baseline from [P1], reflective mode)
Maximise (7) over $\mathcal U(\bar M)$ per group using Riemannian conjugate gradient: Riemannian gradient = projection (4) of (8); direction $\boldsymbol\Xi\leftarrow$ gradient direction + β·(transported previous direction), Polak–Ribière β; retraction [P1 (33)] $(\boldsymbol\Phi+\delta\boldsymbol\Xi)(\mathbf I+\delta^2\boldsymbol\Xi^H\boldsymbol\Xi)^{-1/2}$ (polar-type); backtracking step. (Note [P1] writes a minimisation of $\tilde f=-f$ form, eq. (25a).)

### 5.2 Reciprocal design ([P2] Algorithm 1)
Same FP loop, objective = (7) − ν‖Θ−Θᵀ‖² (ν = 1 in [P2]), manifold $\mathcal U(\bar M)$ **not** $\mathcal U_S$, QR retraction $\mathcal Q(\boldsymbol\Theta+\alpha\boldsymbol\Xi)$ [P2 (40)–(41)], PR-CGA with Armijo (settings: ε=1e−8, ≤8000 iterations, ≤200 Armijo steps, c_init=1, c_dec=0.75, sufficient-increase coefficient 2e−11). Final exact enforcement: $\boldsymbol\Theta_{sym}=(\boldsymbol\Theta+\boldsymbol\Theta^T)/2=\mathbf U\boldsymbol\Sigma\mathbf V^H\Rightarrow\boldsymbol\Theta_{opt}=\mathbf U\mathbf V^H$.

**Things we found to double-check in [P2] (from the web summary; confirm in the PDF):**
1. The PR denominator is printed as $\langle r^{(i)},\Xi^{(i)}\rangle$ rather than the usual $\langle r^{(i)},r^{(i)}\rangle$, and the sign conventions (direction starts at $-r$ while maximising) look inconsistent.
2. The explicit Armijo inequality is not shown.
3. Claim that the reformulated objective "is convex" in Θ is unproven and doubtful (we maximise a quadratic with indefinite structure on a non-convex set).
4. The post-processing $\mathbf U\mathbf V^H$ of a symmetric matrix is *unitary and symmetric*? $\Theta_{sym}$ symmetric ⇒ $\mathbf U\mathbf V^H$ is its polar unitary factor, which is symmetric [DERIVED: polar factor of a symmetric matrix is symmetric]. (Proof sketch via Takagi: $\mathbf A=\mathbf Q\boldsymbol\Sigma\mathbf Q^T\Rightarrow$ polar factor $=\mathbf Q\mathbf Q^T$.) So $\mathbf U\mathbf V^H$ is the unitary matrix nearest to $\boldsymbol\Theta_{sym}$ and lies in $\mathcal U_S$; but it is not guaranteed to be the nearest point of $\mathcal U_S$ to the *un-symmetrised* iterate, and **it changes the sum-rate**; the reported rate must be computed *after* it.

### 5.3 CWSC special case
Diagonal unit-modulus: each $\theta_m=e^{j\varphi_m}$. [P1] Algorithm 3 is closed-form per element (for hybrid mode with amplitude splitting; in reflective mode the amplitude is fixed to 1 and only the phase update remains). **Validation:** unconstrained and reciprocal CWSC must give the same rate.

---
## 6. Reciprocity gap — exact definition (for the report)

Fix architecture $a$, $(N,K,M,G,P_{\max},\text{fading})$. For channel realisation $\omega$ with shared beamformer $\mathbf V(\omega)$:
$$R_{U}(\omega)=\max_{\text{restarts }s}\ R\big(\boldsymbol\Theta^{U}_s(\omega)\big),\qquad R_{S}(\omega)=\max_{s}\ R\big(\boldsymbol\Theta^{S}_s(\omega)\big)$$
with $\boldsymbol\Theta^S$ **evaluated after exact symmetrisation/unitarisation** and all matrices checked: $\|\boldsymbol\Theta^H\boldsymbol\Theta-\mathbf I\|_F\le10^{-10}$ (U) and additionally $\|\boldsymbol\Theta-\boldsymbol\Theta^T\|_F\le10^{-10}$ (S).

$$\Delta_{abs}(\omega)=R_U(\omega)-R_S(\omega)\ \ [\text{bit/s/Hz}],\qquad \Delta_{rel}(\omega)=100\,\frac{R_U(\omega)-R_S(\omega)}{R_U(\omega)}\ \ [\%].$$
Report $\bar\Delta_{abs}=\frac1{L}\sum_\omega\Delta_{abs}(\omega)$ (paired over L realisations, with 95% CI from the paired differences), and $\bar\Delta_{rel}=100\,\bar\Delta_{abs}/\bar R_U$ (ratio of means; also give mean of per-realisation ratios).

**Subtle points to put in the report:**
- Since $\mathcal A^S\subseteq\mathcal A^U$, the *global* optima satisfy $R^\star_U\ge R^\star_S$. The solvers are local, so **a per-realisation $\Delta<0$ can occur**. Therefore (i) use multiple random restarts for both, (ii) report the fraction of realisations with $\Delta<0$, (iii) never clip negative values silently.
- Use *identical* restart budgets, tolerances and iteration caps where the algorithms allow; document where they cannot (ν penalty, QR vs polar retraction).
- Because the algorithms are not the same code, a measured gap mixes *constraint* effect and *algorithm* effect. A cleaner control: also run the unconstrained design with **[P2]'s own machinery (ν = 0)** and the reciprocal design with the same FP + CGA — then the only change is ν and the final projection. We recommend this as the main "controlled" comparison, and [P1]'s own Riemannian-CG as a secondary cross-check. [DESIGN PROPOSAL — not from a paper]
- $\Delta_{CWSC}=0$ exactly (Sec. 2.1).

---
## 7. Validation checklist Member 1 owns (tests Member 4 / 2 / 3 must pass)
1. Path-loss/SNR hand check of Sec. 1.1 (−6.13 dB at 20 dBm) reproduced by code.
2. $\mathbf V$ satisfies $\|\mathbf V\|_F^2=P_{\max}$ to machine precision.
3. Sum-rate evaluator agrees with a brute-force loop implementation of (3).
4. Random $\boldsymbol\Theta=\mathbf Q\mathbf Q^T$ passes unitary+symmetric tests; Gaussian-then-QR passes unitary only.
5. Finite-difference check of gradient (8) on random data (relative error < 1e-6).
6. FP objective (6) equals the sum-rate at $\iota_k=\gamma_k$, $y_k=y_k^\star$ (equality check).
7. Objective non-decreasing across FP outer iterations (to numerical tolerance).
8. CWSC: $\Delta=0$. FC ≥ GC ≥ SC in the unconstrained case on average (ordering of [P1]).
9. [P2] reproduction: Fig. 3 (convergence, K=N=2, R=32, 20 dBm, Rayleigh), Fig. 4/7 trends (GC(4) ≥ SotA etc. — trends only; we do not know numeric values).
10. [P1] reproduction: trends of Fig. 7, 9, 10 (note their hybrid-mode setting vs our reflective-mode setting).
