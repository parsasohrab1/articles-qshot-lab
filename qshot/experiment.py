"""Hamiltonian/state I/O (SRS FR-01, FR-02), experiment runner with RMSE and bias per method (FR-06) and CSV output."""

from __future__ import annotations

import csv
import json
import time
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np

from .estimators import (Hamiltonian, _Rot, estimate_adaptive, estimate_oracle_var, estimate_shadow, estimate_uniform, estimate_weighted, exact_expectation,
                         random_hamiltonian, random_state, term_expectations)

METHODS = ["uniform", "weighted", "oracle_var", "adaptive", "shadow"]


def load_hamiltonian_json(path) -> Hamiltonian:
    """FR-01: {"n_qubits": n, "terms": [{"pauli": "XIZ", "coeff": 0.5}, ...]}."""
    d = json.loads(Path(path).read_text(encoding="utf8"))
    n = int(d["n_qubits"])
    h = []
    for t in d["terms"]:
        p = t["pauli"].upper()
        if len(p) != n or not set(p) <= set("IXYZ"):
            raise ValueError(f"invalid Pauli string: {p!r} for n_qubits={n}")
        h.append((float(t["coeff"]), p))
    if not h:
        raise ValueError("no terms")
    return h


def save_hamiltonian_json(h: Hamiltonian, path) -> None:
    n = len(h[0][1])
    Path(path).write_text(json.dumps({"n_qubits": n, "terms": [{"pauli": p, "coeff": c} for c, p in h]}, indent=1), encoding="utf8")


def run_experiment(h: Hamiltonian, psi: np.ndarray, total_shots: int, repeats: int = 200, seed: int = 0, methods: Sequence[str] = METHODS) -> List[Dict]:
    """FR-03..FR-06: repeat each estimator `repeats` times at a fixed total-shot budget; report bias and RMSE against the exact value."""
    exact = exact_expectation(psi, h)
    e = term_expectations(psi, h)
    rot = _Rot(psi)
    rows = []
    for m in methods:
        rng = np.random.default_rng(seed)
        t0 = time.time()
        vals = np.empty(repeats)
        for r in range(repeats):
            if m == "uniform": vals[r] = estimate_uniform(h, e, total_shots, rng)
            elif m == "weighted": vals[r] = estimate_weighted(h, e, total_shots, rng)
            elif m == "oracle_var": vals[r] = estimate_oracle_var(h, e, total_shots, rng)
            elif m == "adaptive": vals[r] = estimate_adaptive(h, e, total_shots, rng)
            elif m == "shadow": vals[r] = estimate_shadow(h, rot, total_shots, rng)
            else: raise ValueError(m)
        rows.append({"method": m, "total_shots": total_shots, "repeats": repeats, "exact": exact, "bias": float(vals.mean() - exact),
                     "rmse": float(np.sqrt(np.mean((vals - exact) ** 2))), "seconds": round(time.time() - t0, 3)})
    return rows


def write_csv(rows: List[Dict], path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def sweep(h: Hamiltonian, psi: np.ndarray, budgets: Sequence[int], repeats: int = 200, seed: int = 0) -> List[Dict]:
    out = []
    for b in budgets:
        out += run_experiment(h, psi, b, repeats, seed)
    return out


def chemistry_hamiltonian(symbols: Sequence[str], r_bohr: float, charge: int = 0):
    """Real molecular Hamiltonian (STO-3G, Jordan-Wigner) and its exact ground state in the right electron sector. Computed, not measured. Needs PennyLane."""
    import pennylane as qml
    H, n = qml.qchem.molecular_hamiltonian(list(symbols), np.array([0, 0, 0, 0, 0, r_bohr], dtype=float), charge=charge, basis="sto-3g")
    coeffs, ops = H.terms()
    h: Hamiltonian = []
    for c, op in zip(coeffs, ops):
        word = ["I"] * n
        for w, o in zip(op.wires, (op.operands if hasattr(op, "operands") else [op])):
            word[int(w)] = o.name[-1] if o.name.startswith("Pauli") else o.name
        h.append((float(np.real(c)), "".join(word)))
    return h, n
