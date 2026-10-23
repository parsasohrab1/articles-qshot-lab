"""Shot-allocation estimators for Pauli-sum observables (SRS FR-03..FR-05 and the section-4 adaptive module).

Direct methods measure every Pauli term in its own basis, so the number of +1 outcomes of term i with s_i shots is Binomial(s_i, (1 + <P_i>) / 2); this is
exactly equivalent to sampling shots and is used for speed. Classical shadows sample random single-qubit Pauli bases and snapshots explicitly.
Allocation methods:
  uniform   : equal shots per non-identity term
  weighted  : shots proportional to |c_i|
  oracle_var: shots proportional to |c_i| * sigma_i with the TRUE sigma_i (unrealistic reference, requires knowing the state)
  adaptive  : pilot round (uniform), Beta-posterior variance per term, remaining shots proportional to |c_i| * sigma_i (Neyman allocation), counts pooled
  shadow    : classical shadows with random Pauli bases (Huang-Kueng-Preskill)
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

PAULI = {"I": np.eye(2, dtype=complex), "X": np.array([[0, 1], [1, 0]], dtype=complex), "Y": np.array([[0, -1j], [1j, 0]], dtype=complex), "Z": np.diag([1.0, -1.0]).astype(complex)}
BASIS = {"X": np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2), "Y": np.array([[1, -1j], [1, 1j]], dtype=complex) / np.sqrt(2), "Z": np.eye(2, dtype=complex)}

Hamiltonian = List[Tuple[float, str]]


def pauli_matrix(label: str) -> np.ndarray:
    m = np.array([[1]], dtype=complex)
    for ch in label:
        m = np.kron(m, PAULI[ch])
    return m


def term_expectations(psi: np.ndarray, h: Hamiltonian) -> np.ndarray:
    return np.array([np.real(np.vdot(psi, pauli_matrix(p) @ psi)) for _, p in h])


def exact_expectation(psi: np.ndarray, h: Hamiltonian) -> float:
    return float(np.dot([c for c, _ in h], term_expectations(psi, h)))


def _is_identity(label: str) -> bool:
    return set(label) == {"I"}


def _split(h: Hamiltonian):
    idx = [i for i, (_, p) in enumerate(h) if not _is_identity(p)]
    const = float(sum(c for c, p in h if _is_identity(p)))
    return idx, const


def _alloc(total: int, weights: np.ndarray) -> np.ndarray:
    """Integer shots per term proportional to weights (largest-remainder rounding), at least 1 shot per term."""
    w = np.asarray(weights, dtype=float)
    w = w / w.sum() if w.sum() > 0 else np.full(len(w), 1.0 / len(w))
    s = np.maximum(1, np.floor(total * w).astype(int))
    rem = total - s.sum()
    if rem > 0:
        order = np.argsort(-(total * w - np.floor(total * w)))
        s[order[:rem]] += 1
    return s


def _measure(e: float, shots: np.ndarray, rng) -> np.ndarray:
    """Number of +1 outcomes for each term given its exact expectation e and shot counts."""
    return rng.binomial(shots, np.clip((1 + e) / 2, 0.0, 1.0))


def estimate_uniform(h: Hamiltonian, e: np.ndarray, total: int, rng) -> float:
    idx, const = _split(h)
    s = _alloc(total, np.ones(len(idx)))
    k = _measure(e[idx], s, rng)
    return const + float(sum(h[i][0] * (2 * k_j / s_j - 1) for i, k_j, s_j in zip(idx, k, s)))


def estimate_weighted(h: Hamiltonian, e: np.ndarray, total: int, rng) -> float:
    idx, const = _split(h)
    s = _alloc(total, np.array([abs(h[i][0]) for i in idx]))
    k = _measure(e[idx], s, rng)
    return const + float(sum(h[i][0] * (2 * k_j / s_j - 1) for i, k_j, s_j in zip(idx, k, s)))


def estimate_oracle_var(h: Hamiltonian, e: np.ndarray, total: int, rng) -> float:
    idx, const = _split(h)
    sig = np.sqrt(np.clip(1 - e[idx] ** 2, 1e-6, None))
    s = _alloc(total, np.array([abs(h[i][0]) for i in idx]) * sig)
    k = _measure(e[idx], s, rng)
    return const + float(sum(h[i][0] * (2 * k_j / s_j - 1) for i, k_j, s_j in zip(idx, k, s)))


def estimate_adaptive(h: Hamiltonian, e: np.ndarray, total: int, rng, pilot_fraction: float = 0.2) -> float:
    """Section-4 module: pilot round, posterior variance per term (Beta(1,1) prior), Neyman allocation of the rest, pooled counts."""
    idx, const = _split(h)
    n = len(idx)
    pilot_total = max(n * 4, int(pilot_fraction * total))
    pilot_total = min(pilot_total, total - n) if total > 2 * n else total // 2
    s1 = _alloc(pilot_total, np.ones(n))
    k1 = _measure(e[idx], s1, rng)
    p_hat = (k1 + 1.0) / (s1 + 2.0)
    sig = np.sqrt(np.clip(1 - (2 * p_hat - 1) ** 2, 0.02, None))
    s2 = _alloc(max(total - s1.sum(), n), np.array([abs(h[i][0]) for i in idx]) * sig)
    k2 = _measure(e[idx], s2, rng)
    k, s = k1 + k2, s1 + s2
    return const + float(sum(h[i][0] * (2 * k_j / s_j - 1) for i, k_j, s_j in zip(idx, k, s)))


class _Rot:
    """Caches basis-rotated probability vectors per basis string."""

    def __init__(self, psi: np.ndarray):
        self.psi, self.n, self.cache = psi, int(np.log2(len(psi))), {}

    def probs(self, bases: str) -> np.ndarray:
        if bases not in self.cache:
            u = np.array([[1]], dtype=complex)
            for b in bases:
                u = np.kron(u, BASIS[b])
            p = np.abs(u @ self.psi) ** 2
            self.cache[bases] = p / p.sum()
        return self.cache[bases]


def estimate_shadow(h: Hamiltonian, rot: _Rot, total: int, rng) -> float:
    """Classical shadows: `total` snapshots, each = one random Pauli basis per qubit and one measured bitstring; unbiased estimator 3^k * prod(+-1) when bases match."""
    n = rot.n
    letters = np.array(["X", "Y", "Z"])
    basis_idx = rng.integers(0, 3, size=(total, n))
    bits = np.zeros((total, n), dtype=np.int8)
    keys = ["".join(letters[row]) for row in basis_idx]
    for key in set(keys):
        sel = np.array([i for i, k in enumerate(keys) if k == key])
        outcome = rng.choice(2 ** n, size=len(sel), p=rot.probs(key))
        bits[sel] = ((outcome[:, None] >> np.arange(n - 1, -1, -1)) & 1).astype(np.int8)
    sign = 1 - 2 * bits.astype(float)
    letter_idx = {"X": 0, "Y": 1, "Z": 2}
    est = 0.0
    for c, p in h:
        if _is_identity(p):
            est += c
            continue
        val = np.ones(total)
        for q, ch in enumerate(p):
            if ch != "I":
                val = val * np.where(basis_idx[:, q] == letter_idx[ch], 3.0 * sign[:, q], 0.0)
        est += c * float(val.mean())
    return float(est)


def random_hamiltonian(n: int, n_terms: int, rng, max_weight: int = 2) -> Hamiltonian:
    import itertools
    labels = ["".join(t) for t in itertools.product("IXYZ", repeat=n)]
    labels = [l for l in labels if 0 < sum(ch != "I" for ch in l) <= max_weight]
    pick = rng.choice(len(labels), size=min(n_terms, len(labels)), replace=False)
    return [(float(rng.normal()), labels[i]) for i in pick]


def random_state(n: int, rng) -> np.ndarray:
    v = rng.normal(size=2 ** n) + 1j * rng.normal(size=2 ** n)
    return v / np.linalg.norm(v)
