"""
M1.py - BD-RIS Reciprocity Gap Study
ECE 310 Wireless Communications | Group 9 | Ahmedabad University

M1 Scope: Channel generation + basic sum-rate computation.
Full optimisers will be implemented in M2.

Run demo : python M1.py --demo
Run tests: python M1.py --test
"""

import numpy as np
import argparse
import sys

# ─────────────────────────────────────────────────────────────────────────────
# PARAMETERS  (from project proposal)
# ─────────────────────────────────────────────────────────────────────────────

N = 4          # BS antennas
K = 4          # users
M = 16         # RIS elements

# Path loss: PL(d) = zeta_0 * (d0/d)^epsilon
ZETA_0    = 10 ** (-30 / 10)   # reference gain at d0 = 1m  (-30 dB)
D0        = 1.0                 # reference distance (m)
EPSILON   = 2.2                 # path loss exponent

D_BS_RIS  = 50.0    # BS  → RIS distance (m)
D_RIS_USR = 2.5     # RIS → user distance (m)
D_BS_USR  = 52.0    # BS  → user distance (m)  direct link

NOISE_W   = 10 ** ((-80 - 30) / 10)   # -80 dBm noise in Watts
PMAX_W    = 10 ** ((10 - 30) / 10)    # 10 dBm transmit power in Watts

SEED      = 42


def path_loss(d):
    """Large-scale path loss: zeta_0 * (d0/d)^epsilon"""
    return ZETA_0 * (D0 / d) ** EPSILON


# ─────────────────────────────────────────────────────────────────────────────
# CHANNEL GENERATION
# ─────────────────────────────────────────────────────────────────────────────

def generate_channel(rows, cols, distance, rng):
    """
    Generate one Rayleigh fading channel matrix of shape (rows, cols).

    Each entry is CN(0, 1) — circularly symmetric complex Gaussian.
    Scaled by sqrt(path_loss(distance)) for large-scale attenuation.
    """
    real = rng.standard_normal((rows, cols))
    imag = rng.standard_normal((rows, cols))
    small_scale = (real + 1j * imag) / np.sqrt(2)   # CN(0,1)
    return np.sqrt(path_loss(distance)) * small_scale


def generate_all_channels(rng):
    """
    Generate the three channel matrices for one realisation.

    Returns:
        G   (M x N) : BS  → RIS channel
        H_r (K x M) : RIS → users channel
        H_d (K x N) : BS  → users direct channel
    """
    G   = generate_channel(M, N, D_BS_RIS,  rng)
    H_r = generate_channel(K, M, D_RIS_USR, rng)
    H_d = generate_channel(K, N, D_BS_USR,  rng)
    return G, H_r, H_d


# ─────────────────────────────────────────────────────────────────────────────
# SUM-RATE COMPUTATION
# ─────────────────────────────────────────────────────────────────────────────

def effective_channel(G, H_r, H_d, Theta):
    """
    Compute the effective channel each user sees.

    H_eff[k] = H_d[k] + H_r[k] @ Theta @ G

    This combines the direct path and the RIS-reflected path.
    Theta (M x M) is the BD-RIS scattering matrix we will optimise in M2.
    """
    return H_d + H_r @ Theta @ G    # shape: (K, N)


def mrt_beamformer(H_eff):
    """
    Maximum Ratio Transmission: point each beam at its user.

    w_k = h_eff_k* / ||h_eff_k||   (matched filter)
    Normalised so total power = PMAX_W.
    """
    W = H_eff.conj().T                                      # (N, K)
    norms = np.linalg.norm(W, axis=0, keepdims=True)
    W = W / np.where(norms < 1e-12, 1.0, norms)            # unit columns
    W = W * np.sqrt(PMAX_W / (np.linalg.norm(W, 'fro') ** 2))
    return W


def compute_sum_rate(G, H_r, H_d, Theta):
    """
    Compute sum-rate (bps/Hz) for a given scattering matrix Theta.

    Steps:
      1. Effective channel H_eff = H_d + H_r @ Theta @ G
      2. MRT beamformer W from H_eff
      3. SINR_k = signal power / (interference + noise)
      4. R = sum_k  log2(1 + SINR_k)
    """
    H_eff = effective_channel(G, H_r, H_d, Theta)
    W     = mrt_beamformer(H_eff)

    sum_rate = 0.0
    for k in range(K):
        gains        = np.abs(H_eff[k].conj() @ W) ** 2   # power from each beam
        signal       = gains[k]
        interference = gains.sum() - signal
        sinr_k       = signal / (interference + NOISE_W)
        sum_rate     += np.log2(1 + sinr_k)

    return sum_rate


# ─────────────────────────────────────────────────────────────────────────────
# PLACEHOLDER SCATTERING MATRICES  (to be replaced by optimisers in M2)
# ─────────────────────────────────────────────────────────────────────────────

def identity_theta():
    """Theta = I  (RIS has no effect — baseline reference)."""
    return np.eye(M, dtype=complex)


def random_unitary_theta(rng):
    """
    Random unitary Theta — placeholder for the unconstrained design (Li et al.).
    In M2 this will be the output of the FP + Stiefel manifold optimiser.
    """
    A = rng.standard_normal((M, M)) + 1j * rng.standard_normal((M, M))
    Q, R = np.linalg.qr(A)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


def random_symmetric_unitary_theta(rng):
    """
    Random symmetric unitary Theta — placeholder for the reciprocal design (Fidanovski et al.).
    Construction: Theta = U @ diag(e^{i*phi}) @ U.T
    Satisfies both Theta^H Theta = I  and  Theta = Theta.T
    In M2 this will be the output of the FP + symmetric manifold optimiser.
    """
    A = rng.standard_normal((M, M)) + 1j * rng.standard_normal((M, M))
    Q, _ = np.linalg.qr(A)
    phi   = rng.uniform(0, 2 * np.pi, M)
    return Q @ np.diag(np.exp(1j * phi)) @ Q.T


# ─────────────────────────────────────────────────────────────────────────────
# DEMO
# ─────────────────────────────────────────────────────────────────────────────

def run_demo():
    rng = np.random.default_rng(SEED)
    G, H_r, H_d = generate_all_channels(rng)

    Theta_I   = identity_theta()
    Theta_unc = random_unitary_theta(rng)
    Theta_rec = random_symmetric_unitary_theta(rng)

    sr_I   = compute_sum_rate(G, H_r, H_d, Theta_I)
    sr_unc = compute_sum_rate(G, H_r, H_d, Theta_unc)
    sr_rec = compute_sum_rate(G, H_r, H_d, Theta_rec)

    print("\n" + "=" * 52)
    print("  BD-RIS Reciprocity Gap — M1 Demo  (Group 9)")
    print("=" * 52)
    print(f"  N={N} antennas  |  K={K} users  |  M={M} RIS elements")
    print(f"  P_max = 10 dBm  |  Noise = -80 dBm  |  Rayleigh fading")
    print(f"  BS-RIS = {D_BS_RIS}m  |  RIS-user = {D_RIS_USR}m\n")
    print(f"  {'Theta type':<30}  {'Sum-rate (bps/Hz)':>18}")
    print(f"  {'-'*30}  {'-'*18}")
    print(f"  {'Identity (no RIS effect)':<30}  {sr_I:>18.4f}")
    print(f"  {'Random unitary (unconstrained)':<30}  {sr_unc:>18.4f}")
    print(f"  {'Random symm-unitary (reciprocal)':<30}  {sr_rec:>18.4f}")
    print()
    print("  NOTE: Theta values are random placeholders.")
    print("  Optimised results will appear from M2 onwards.")
    print("=" * 52 + "\n")


# ─────────────────────────────────────────────────────────────────────────────
# TESTS
# ─────────────────────────────────────────────────────────────────────────────

def run_tests():
    passed = failed = 0

    def check(name, condition, detail=""):
        nonlocal passed, failed
        if condition:
            print(f"  PASS  {name}")
            passed += 1
        else:
            print(f"  FAIL  {name}" + (f" — {detail}" if detail else ""))
            failed += 1

    rng = np.random.default_rng(SEED)
    G, H_r, H_d = generate_all_channels(rng)

    # Channel shapes
    check("G shape is (M, N)",   G.shape   == (M, N))
    check("H_r shape is (K, M)", H_r.shape == (K, M))
    check("H_d shape is (K, N)", H_d.shape == (K, N))

    # Channels are complex
    check("G is complex",   np.iscomplexobj(G))
    check("H_r is complex", np.iscomplexobj(H_r))
    check("H_d is complex", np.iscomplexobj(H_d))

    # Path loss: BS-RIS (50m) should be weaker than RIS-user (2.5m)
    check("BS-RIS channel weaker than RIS-user (path loss)",
          np.mean(np.abs(G)**2) < np.mean(np.abs(H_r)**2))

    # Different seeds give different channels
    G2, _, _ = generate_all_channels(np.random.default_rng(99))
    check("Different seeds give different channels", not np.allclose(G, G2))

    # Same seed gives same channels (reproducibility)
    Ga, _, _ = generate_all_channels(np.random.default_rng(SEED))
    Gb, _, _ = generate_all_channels(np.random.default_rng(SEED))
    check("Same seed gives same channel", np.allclose(Ga, Gb))

    # Effective channel shape
    Theta = identity_theta()
    H_eff = effective_channel(G, H_r, H_d, Theta)
    check("Effective channel shape is (K, N)", H_eff.shape == (K, N))

    # With identity Theta, H_eff = H_d + H_r @ G
    expected = H_d + H_r @ np.eye(M) @ G
    check("Identity Theta: H_eff = H_d + H_r @ G", np.allclose(H_eff, expected))

    # Sum-rate is positive and finite
    sr = compute_sum_rate(G, H_r, H_d, Theta)
    check("Sum-rate is positive",  sr > 0,        f"got {sr:.4f}")
    check("Sum-rate is finite",    np.isfinite(sr), f"got {sr}")

    # Unitary Theta: Theta^H Theta = I
    Theta_unc = random_unitary_theta(np.random.default_rng(1))
    err = np.linalg.norm(Theta_unc.conj().T @ Theta_unc - np.eye(M))
    check("Unitary Theta satisfies Theta^H Theta = I", err < 1e-6, f"err={err:.2e}")

    # Symmetric unitary Theta: both Theta^H Theta = I and Theta = Theta.T
    Theta_rec = random_symmetric_unitary_theta(np.random.default_rng(2))
    err_u = np.linalg.norm(Theta_rec.conj().T @ Theta_rec - np.eye(M))
    err_s = np.linalg.norm(Theta_rec - Theta_rec.T)
    check("Symm-unitary Theta is unitary",   err_u < 1e-6, f"err={err_u:.2e}")
    check("Symm-unitary Theta is symmetric", err_s < 1e-6, f"err={err_s:.2e}")

    # Summary
    print(f"\n  {passed} passed  |  {failed} failed  |  {passed + failed} total\n")
    return failed == 0


# ─────────────────────────────────────────────────────────────────────────────
# ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="BD-RIS M1 — Group 9")
    parser.add_argument("--demo", action="store_true", help="Run demo")
    parser.add_argument("--test", action="store_true", help="Run tests")
    args = parser.parse_args()

    if args.test:
        print("\n  BD-RIS M1 Test Suite\n  " + "─" * 35)
        ok = run_tests()
        sys.exit(0 if ok else 1)

    if args.demo:
        run_demo()

    if not args.test and not args.demo:
        parser.print_help()