"""Estimators of <H> = sum_i c_i <P_i> under a total budget of `shots` single-shot measurements.

Each Pauli term is its own measurement setting for `uniform`, `weighted`, `oracle` and `adaptive`: a term measured s times returns
Binomial(s, (1 + e)/2) outcomes of +1, with e the exact expectation, so the simulation needs no state-vector sampling.
`shadows` draws random Pauli bases for all qubits at once (classical shadows) and reuses every snapshot for every term.
All functions return an (reps, n_terms) array of term estimates and the allocation used; the energy is estimates @ coeffs.
"""
import numpy as np

_BASES = {"X": np.array([[1, 1], [1, -1]], complex) / np.sqrt(2),                     # H
          "Y": np.array([[1, 1], [1j, -1j]], complex) / np.sqrt(2),                    # H S^dagger
          "Z": np.eye(2, dtype=complex)}


def allocate(weights, shots, min_each=1):
    """Integer allocation proportional to `weights` with at least `min_each` shots per term; the total is exactly `shots`."""
    w = np.asarray(weights, float)
    m = len(w)
    assert shots >= m * min_each, "budget smaller than one shot per term"
    base = np.full(m, min_each)
    rest = shots - base.sum()
    share = w / w.sum() * rest
    add = np.floor(share).astype(int)
    left = rest - add.sum()
    add[np.argsort(-(share - add))[:left]] += 1
    return base + add


def _measure(e, s, rng, reps):
    k = rng.binomial(s[None, :], ((1 + e) / 2)[None, :].repeat(reps, 0))
    return 2 * k / s[None, :] - 1, k


def uniform(e, c, shots, rng, reps):
    s = allocate(np.ones(len(c)), shots)
    return _measure(e, s, rng, reps)[0], s


def weighted(e, c, shots, rng, reps):
    """Shots proportional to |c_i| (needs no knowledge of the state)."""
    s = allocate(np.abs(c), shots)
    return _measure(e, s, rng, reps)[0], s


def oracle(e, c, shots, rng, reps):
    """Optimal allocation |c_i| * sigma_i with the true variances; a bound that needs knowledge the experimenter does not have."""
    s = allocate(np.abs(c) * np.sqrt(np.maximum(1 - e ** 2, 1e-6)), shots)
    return _measure(e, s, rng, reps)[0], s


def adaptive(e, c, shots, rng, reps, pilot=0.2):
    """Pilot round with weighted allocation, posterior variance per term (Beta(1,1) prior), remaining budget split by |c_i| * sigma_i."""
    m = len(c)
    s1 = allocate(np.abs(c), max(int(shots * pilot), m))
    k1 = rng.binomial(s1[None, :], ((1 + e) / 2)[None, :].repeat(reps, 0))
    p = (k1 + 1) / (s1[None, :] + 2)
    sigma = np.sqrt(4 * p * (1 - p))                                 # per-repetition posterior estimate of the std of a +-1 outcome
    rest = shots - s1.sum()
    est = np.empty((reps, m))
    tot = np.empty(m)
    for r in range(reps):
        s2 = allocate(np.abs(c) * sigma[r] + 1e-9, rest, min_each=0) if rest > 0 else np.zeros(m, int)
        k2 = rng.binomial(s2, (1 + e) / 2)
        est[r] = 2 * (k1[r] + k2) / (s1 + s2) - 1
        tot = s1 + s2
    return est, tot


def shadows(e, c, shots, rng, reps, terms=None, psi=None):
    """Classical shadows with uniformly random Pauli bases; every snapshot estimates every term (estimator 3^w * product of outcomes when the bases match)."""
    n = len(terms[0])
    nb = 3 ** n
    combos = np.array(np.meshgrid(*[np.arange(3)] * n, indexing="ij")).reshape(n, -1).T       # (3^n, n) basis choices, X/Y/Z = 0/1/2
    letters = "XYZ"
    prob = np.empty((nb, 2 ** n))
    for j, cb in enumerate(combos):
        U = np.array([[1.0 + 0j]])
        for q in cb:
            U = np.kron(U, _BASES[letters[q]].conj().T)
        prob[j] = np.abs(U @ psi) ** 2
    q = (prob / nb).ravel()
    q = q / q.sum()
    outcomes = ((np.arange(2 ** n)[:, None] >> (n - 1 - np.arange(n))) & 1)                      # bit of each qubit, qubit 0 = most significant
    sign = 1 - 2 * outcomes                                                                       # (2^n, n) of +-1
    V = np.zeros((nb * 2 ** n, len(terms)))
    for t, p in enumerate(terms):
        sup = [i for i, ch in enumerate(p) if ch != "I"]
        want = np.array([letters.index(p[i]) for i in sup])
        match = (combos[:, sup] == want).all(1)                                                   # (3^n,)
        val = 3.0 ** len(sup) * np.prod(sign[:, sup], axis=1)                                     # (2^n,)
        V[:, t] = (match[:, None] * val[None, :]).ravel()
    counts = rng.multinomial(shots, q, size=reps)
    return counts @ V / shots, np.full(len(terms), shots / len(terms))


STRATEGIES = {"uniform": uniform, "weighted": weighted, "oracle": oracle, "adaptive": adaptive, "shadows": shadows}
