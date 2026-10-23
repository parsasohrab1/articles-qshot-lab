"""Measured experiments: (a) random Pauli Hamiltonians and random states (SIMULATED inputs), (b) real molecular Hamiltonians (computed with PennyLane, STO-3G) with their exact ground states."""
import json
from pathlib import Path

import numpy as np

from qshot import chemistry_hamiltonian, random_hamiltonian, random_state, run_experiment
from qshot.cli import _kron_label


def ground_state(h, n, electrons=2):
    pop = np.array([bin(i).count("1") for i in range(2 ** n)])
    mat = sum(c * _kron_label(p) for c, p in h) + 50.0 * np.diag((pop - electrons) ** 2.0)
    return np.linalg.eigh(mat)[1][:, 0].astype(complex)


def main():
    out = {"random_3q_12terms": {}, "chemistry": {}}
    for shots in (500, 2000, 8000):
        acc = {}
        for seed in range(5):
            rng = np.random.default_rng(seed)
            h, psi = random_hamiltonian(3, 12, rng), random_state(3, rng)
            for r in run_experiment(h, psi, shots, repeats=100, seed=seed):
                acc.setdefault(r["method"], []).append(r["rmse"])
        out["random_3q_12terms"][str(shots)] = {m: round(float(np.mean(v)), 4) for m, v in acc.items()}
    for name, sym, ch, rs in (("H2", ["H", "H"], 0, [1.4, 3.0]), ("HeH+", ["He", "H"], 1, [1.5])):
        for r_b in rs:
            h, n = chemistry_hamiltonian(sym, r_b, ch)
            psi = ground_state(h, n)
            res = {}
            for shots in (500, 2000, 8000):
                res[str(shots)] = {r["method"]: round(r["rmse"], 5) for r in run_experiment(h, psi, shots, repeats=200, seed=0)}
            out["chemistry"][f"{name}_R{r_b}"] = {"n_terms": len(h), "rmse": res}
    Path("results").mkdir(exist_ok=True)
    Path("results/experiments.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
