# SOTA Synthesis, Gap Statement and Positioning — Group 9 (M2, Member 1)

*Search date: 10 Oct 2026. Verification levels:* **FULL** = read the PDF · **ARXIV-FETCH** = read arXiv abstract/HTML through a web-fetch summary (must be re-read before quoting numbers) · **LISTED-ONLY** = taken from our M1 reference list, not re-checked.

## 1. Synthesis (organised by the questions the project cares about)

**(a) Physical model → constraint.** A passive lossless multiport has a unitary scattering matrix, $\Theta^H\Theta=I$. Li–Shen–Clerckx (TWC 2023, **FULL**) build the nine mode/architecture BD-RIS framework on this and impose *nothing else* (their text never mentions reciprocity). Reciprocity — no non-reciprocal components — adds $\Theta=\Theta^T$ (Shen–Clerckx–Murch, TWC 2022, **LISTED-ONLY**). Hence the "unconstrained" design is an upper-bound-style benchmark that needs non-reciprocal hardware (circulators etc.), not necessarily "active" hardware as M1 phrases it.

**(b) Architecture ↔ structure.** SC = diagonal, GC = block-diagonal, FC = full. Under reciprocity SC is *unchanged* (diagonal ⇒ symmetric), so any reciprocity cost can only appear for GC (block size ≥ 2) and FC. Real degrees of freedom: $n^2$ (unitary) vs $n(n+1)/2$ (symmetric unitary) per block.

**(c) Optimisation families.**
1. *FP + manifold CG* — Li et al. (FULL) for $\mathcal U(\cdot)$; Fidanovski et al. 2509.20246 (ARXIV-FETCH) for reciprocal, via quadratic penalty + final symmetrise/SVD. Same LDT+QT outer loop in both.
2. *Exact symmetric-unitary manifold* — Santamaria et al. 2604.22469 (tangent space, retraction, geodesics, low-rank result when $M>N_t+N_r$) and Fidanovski et al. 2607.22170 (blockwise tangent projections + Takagi retractions for GC).
3. *Projection / ADMM* — Wu & Clerckx (arXiv 2412.15965, LISTED-ONLY), Zhou et al. 2509.17804 (SOSUP; reciprocal-only; ARXIV-FETCH).
4. *Closed form* — Yahya et al. (OJ-COMS 2024, LISTED-ONLY) passive MRT / interference nulling.

**(d) What is known about the *cost* of reciprocity.**
- Fidanovski 2509.20246: evaluates only reciprocal designs against SotA [Yahya; Wu–Clerckx]; **no** non-reciprocal U(R) run reported (ARXIV-FETCH). Uses fixed MMSE beamformer in some figures.
- Fidanovski et al. 2607.22170 (ARXIV-FETCH; date not stated in what we read): *does* compare reciprocal vs non-reciprocal (Manopt over U(R)) for SC/GC/FC; reports the curves "nearly coincide"/"no noticeable loss"; setup N_t=N_r=4, P_max=0 dBm, 100 realisations, sweep over R only, isotropic signalling (no beamformer), no numeric gap; authors say other system models are left for a journal extension and cite earlier reports of non-reciprocal gains in interference-coupled settings.
- Santamaria et al. 2604.22469: optimising on $\mathcal U(M)$ then projecting to symmetric causes a "significant rate penalty" — i.e. *how* reciprocity is enforced matters; no direct $\mathcal U$ vs $\mathcal U_S$ gap number (ARXIV-FETCH).
- Zhou et al. 2509.17804: states reciprocity prevents reaching the theoretical upper bound when DoF>1; no gap measurement.

## 2. Comparison table

| # | Work | Venue / ID | Verif. | Objective / setting | Constraint on Θ | Architectures | Method | Evidence on reciprocity cost | Limitation for us |
|---|---|---|---|---|---|---|---|---|---|
| 1 | Li, Shen, Clerckx | IEEE TWC 22(4):2311–2324, Apr 2023 | FULL | max sum-rate, joint W & Φ, MU-MISO, direct link blocked | unitary only | 9 mode×arch (CWSC/GC/FC) | FP (LDT+QT) + BCD; Riemannian PR-CG (Alg. 2); closed form CWSC (Alg. 3) | none (no reciprocity) | printed tangent/projection (28),(30) look non-standard; hybrid-mode results; joint design |
| 2 | Fidanovski et al. | arXiv 2509.20246 (M1 cites IEEE TWC vol 25, pp 20463–20477 — **unverified**) | ARXIV-FETCH | max sum-rate, MU-MISO, N=K=2, R=32, 0–20 dBm | symmetric + unitary | SC, GC(2), GC(4), FC | FP (LDT+QT) + penalty ν‖Θ−Θᵀ‖² on Stiefel, CGA+Armijo, final sym+SVD | none vs unconstrained | benchmarks only reciprocal SotA |
| 3 | Fidanovski, …, Clerckx | arXiv 2511.07683 (Nov 2025) | abstract only | FP + manifold for reciprocal BD-RIS | symmetric unitary | not seen | LDT+QT, manifold | not seen | read full text |
| 4 | Fidanovski, Venkataramanaiah, Jorswieck, de Abreu | arXiv 2607.22170 | ARXIV-FETCH | sum-rate, group-connected reciprocal BD-RIS MIMO; isotropic signalling | symmetric unitary blocks | SC/GC/FC | blockwise tangent proj. + Takagi retraction | **reciprocal ≈ non-reciprocal** (graphical, N=4, 0 dBm) | no beamformer, one power, no numbers |
| 5 | Santamaria et al. | arXiv 2604.22469 | ARXIV-FETCH | Riemannian opt. on $\mathcal U_S$; sum-gain, rate, MSE | symmetric unitary | FC, low-rank stem | geodesics, LS/PO | penalty of "optimise on U then project" | no U vs U_S gap number |
| 6 | Zhou, Fang, Mao, Clerckx | arXiv 2509.17804 | ARXIV-FETCH | stem/cluster-connected, SOSUP | symmetric unitary | many | projection | theoretical upper bound unreachable | reciprocal only |
| 7 | Morales Sandoval et al. | arXiv 2605.16006 | ARXIV-FETCH | cell-free MIMO + beamforming | reciprocal | – | FP+MO | none | different system |
| 8 | Shen, Clerckx, Murch | IEEE TWC 21(2):1229–1243, 2022 | LISTED-ONLY | S-parameter RIS model | unitary / symmetric | SC/GC/FC | theory | physical basis of reciprocity | not read by us |
| 9 | Yahya et al. | IEEE OJ-COMS 5:7613–7627, 2024 | LISTED-ONLY | passive MRT & nulling | – | – | closed form | – | – |
| 10 | Wu & Clerckx | arXiv 2412.15965 | LISTED-ONLY | universal framework (pp-ADMM) | general | arbitrary | ADMM | – | – |
| 11 | Nerini & Clerckx | IEEE CL 27(10):2842–2846, 2023 | LISTED-ONLY | performance–complexity Pareto | – | SC/GC/FC | analysis | – | – |
| 12 | Shen & Yu | IEEE TSP 66(10):2616–2630, 2018 | LISTED-ONLY | FP theory | – | – | LDT/QT | – | – |
| 13 | Boumal | CUP 2023 | LISTED-ONLY | manifold optimisation text | – | – | – | – | – |

**Work still to read before the report:** 2511.07683 (full), 2607.22170 (full, and its date/ journal-extension status), 2604.22469, Shen–Clerckx–Murch. A non-reciprocal BD-RIS full-duplex paper (arXiv 2411.18523) appeared in search but we did not read it.

## 3. Gap statement (revised — replaces M1's "no published work compares…")

**Established by the literature we checked:** the *qualitative* statement "reciprocal and non-reciprocal BD-RIS curves nearly coincide in one setup" exists (#4). So the claims "unanswered" and "first characterisation" in M1 are **not defensible** as written.

**Still open to the best of our knowledge (hypotheses to test, not facts):**
- H1: A *quantified* gap (bit/s/Hz and %, with confidence intervals) across $M$, $G$, $P_{\max}$ for CWGC and CWFC.
- H2: The gap under a **fixed MMSE beamformer** (the existing comparison uses isotropic signalling).
- H3: Matched Rayleigh *and* Rician evaluation; behaviour at $N=K\in\{2,\dots,6\}$ and $M\le 64$ (where #5 predicts low-rank/saturation effects).
- H4: Whether the gap is sensitive to *how* reciprocity is enforced (penalty + projection of #2 vs exact manifold of #4/#5).
- H5: Performance-vs-complexity trade-off for reciprocal vs non-reciprocal at each architecture (circuit complexity counts: $M$, $G\bar M(\bar M+1)/2$, $M(M+1)/2$ vs unitary variants).

**Draft statement for the report:**
> The sum-rate cost of imposing $\Theta=\Theta^T$ on a lossless BD-RIS has so far been reported only graphically, for one operating point and without transmit beamforming (arXiv 2607.22170). We quantify this cost in bit/s/Hz and per cent for CWSC, CWGC and CWFC under a common fixed MMSE beamformer, matched Rayleigh/Rician channels and matched restarts, as a function of $M$, $G$ and $P_{\max}$, and test whether the near-zero gap persists outside the regime already reported. CWSC serves as an exact-zero control. A null result (gap ≈ 0) is a valid finding.

## 4. How our work differs (for the SOTA-positioning slide)
| | Li 2023 | Fidanovski 2025 | Fidanovski 2607 | **Ours** |
|---|---|---|---|---|
| Reciprocity enforced | no | yes (penalty) | yes (exact) | **both, side by side** |
| Non-reciprocal benchmark in same paper | n/a | no | yes (graphical) | **yes, numeric** |
| Beamformer | joint | MMSE/uniform | isotropic | **fixed MMSE, identical** |
| Fading | Rayleigh+Rician | Rayleigh+Rician | not stated | **both** |
| Sweeps | P, M | P, R | R only | **P, M, G, N=K** |
| Output | rates | rates | curves | **Δ_abs, Δ_rel, CI, complexity** |

## 5. Corrections to carry into the final M1→M2 text
1. [2] uses LDT+QT with a symmetry penalty on the *ordinary* Stiefel manifold — not "Dinkelbach" and not a "symmetric Stiefel" manifold.
2. [1] contains no reciprocity discussion; reciprocity ⇒ symmetric comes from Shen–Clerckx–Murch.
3. "Requires active hardware" → "requires non-reciprocal components".
4. Verify the IEEE TWC volume/page numbers for [2] before submission (we only confirmed the arXiv preprint).
5. Define group count vs group size once.
6. Match M1's system model (direct link) to what is actually simulated (blocked direct link).
