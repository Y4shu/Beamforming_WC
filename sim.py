"""
Section VI of lit_review.tex: reproduce the look-direction-mismatch example of
Vorobyov, Gershman & Luo (2003), output SINR versus input SNR.

Run:    .venv/bin/python sim.py
Output: figs/sinr_vs_snr.png
"""

import os
import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt

# Scenario (Section VI-A; follows the setup of Vorobyov et al. 2003)
N = 10                       # ULA elements, half-wavelength spacing
PRESUMED, TRUE = 2.0, 5.0    # desired-signal direction: presumed vs actual (degrees)
INTERFERERS = [30.0, 50.0]   # degrees
INR_DB = 30.0
K = 30                       # snapshots
TRIALS = 100
EPS = 3.0                    # worst-case uncertainty bound ||a - a_bar|| <= EPS
GAMMA = 10.0                 # diagonal loading level (noise power = 1)

rng = np.random.default_rng(0)


def steering(theta):
    """a(theta), eq. (1)."""
    return np.exp(-1j * np.pi * np.arange(N) * np.sin(np.radians(theta)))


def cgauss(*shape):
    return (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / np.sqrt(2)


a_bar, a0 = steering(PRESUMED), steering(TRUE)
R_in = np.eye(N) + sum(10 ** (INR_DB / 10) * np.outer(steering(t), steering(t).conj())
                       for t in INTERFERERS)


def sample_cov(snr_db):
    """K snapshots of eq. (2), desired signal included; returns R_hat."""
    X = cgauss(N, K) + np.sqrt(10 ** (snr_db / 10)) * np.outer(a0, cgauss(K))
    for t in INTERFERERS:
        X += np.sqrt(10 ** (INR_DB / 10)) * np.outer(steering(t), cgauss(K))
    return X @ X.conj().T / K


def sinr(w, snr_db):
    """Output SINR, eq. (4), with the TRUE a0 and R_in."""
    return 10 ** (snr_db / 10) * abs(w.conj() @ a0) ** 2 / np.real(w.conj() @ R_in @ w)


def mvdr(R):
    """SMI, eq. (8)."""
    return np.linalg.solve(R, a_bar)


def worst_case(R):
    """Vorobyov et al. SOCP, eq. (15): min ||R^{1/2} w|| s.t. Re(w^H a_bar) >= EPS ||w|| + 1."""
    U = np.linalg.cholesky(R).conj().T
    w = cp.Variable(N, complex=True)
    cp.Problem(cp.Minimize(cp.norm(U @ w)),
               [cp.real(a_bar.conj() @ w) >= EPS * cp.norm(w) + 1]).solve()
    return w.value


methods = {
    "SMI": mvdr,
    "Diagonal loading": lambda R: mvdr(R + GAMMA * np.eye(N)),
    "Worst-case (SOCP)": worst_case,
}

snrs = np.arange(-20, 31, 5)
curves = {name: [] for name in methods}
for snr in snrs:
    Rs = [sample_cov(snr) for _ in range(TRIALS)]
    for name, f in methods.items():
        curves[name].append(10 * np.log10(np.mean([sinr(f(R), snr) for R in Rs])))
    print(f"SNR {snr:3d} dB: " + ", ".join(f"{n} {c[-1]:.1f}" for n, c in curves.items()))

optimal = [10 * np.log10(10 ** (s / 10) * np.real(a0.conj() @ np.linalg.solve(R_in, a0)))
           for s in snrs]

os.makedirs("figs", exist_ok=True)
plt.plot(snrs, optimal, "k--", label="Optimal")
for name, c in curves.items():
    plt.plot(snrs, c, marker="o", label=name)
plt.xlabel("Input SNR (dB)")
plt.ylabel("Output SINR (dB)")
plt.title(f"Mismatch {TRUE - PRESUMED:.0f}°, K = {K}, {TRIALS} trials")
plt.grid(True)
plt.legend()
plt.savefig("figs/sinr_vs_snr.png", dpi=150)
