"""Pauli-sum Hamiltonians, random instances, exact expectation values (state-vector, numpy only)."""
import json

import numpy as np

_P = {"I": np.eye(2, dtype=complex), "X": np.array([[0, 1], [1, 0]], complex), "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.diag([1, -1]).astype(complex)}


def matrix(pauli):
    m = np.array([[1.0 + 0j]])
    for c in pauli:
        m = np.kron(m, _P[c])
    return m


def load_hamiltonian(path):
    """JSON {"n_qubits": n, "terms": [{"pauli": "XIZ", "coeff": 0.5}, ...]}"""
    h = json.load(open(path))
    return h["n_qubits"], [t["pauli"] for t in h["terms"]], np.array([t["coeff"] for t in h["terms"]], float)


def random_hamiltonian(n_qubits, n_terms, seed=0, spread=1.0):
    """Distinct non-identity Pauli strings. `spread` > 0 makes coefficients log-normally spread (the case where allocation matters)."""
    rng = np.random.default_rng(seed)
    terms = set()
    while len(terms) < n_terms:
        s = "".join(rng.choice(list("IXYZ"), n_qubits))
        if set(s) != {"I"}:
            terms.add(s)
    terms = sorted(terms)
    c = rng.choice([-1, 1], n_terms) * np.exp(rng.normal(0, spread, n_terms)) * 0.5
    return terms, c


def random_state(n_qubits, seed=7):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=2 ** n_qubits) + 1j * rng.normal(size=2 ** n_qubits)
    return v / np.linalg.norm(v)


def expectation(psi, pauli):
    return float(np.real(np.vdot(psi, matrix(pauli) @ psi)))


def exact_value(psi, terms, coeffs):
    e = np.array([expectation(psi, t) for t in terms])
    return e, float(np.dot(coeffs, e))
