"""python -m qshot run --hamiltonian h.json [--state s.json | --seed 7] --shots 2000 --out results/run.csv   |   python -m qshot random --qubits 3 --terms 12 --out h.json
python -m qshot chem --molecule H2 --r 1.4 --shots 2000 --out results/h2.csv   (real molecular Hamiltonian, computed with PennyLane)"""

import argparse
import json

import numpy as np

from .estimators import random_hamiltonian, random_state
from .experiment import chemistry_hamiltonian, load_hamiltonian_json, run_experiment, save_hamiltonian_json, write_csv


def _print(rows):
    for r in rows:
        print(f"{r['method']:11s} rmse={r['rmse']:.4f} bias={r['bias']:+.4f} ({r['seconds']}s)")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="qshot")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--hamiltonian", required=True)
    r.add_argument("--state")
    r.add_argument("--seed", type=int, default=7)
    r.add_argument("--shots", type=int, default=2000)
    r.add_argument("--repeats", type=int, default=200)
    r.add_argument("--out", default="results/run.csv")
    g = sub.add_parser("random")
    g.add_argument("--qubits", type=int, default=3)
    g.add_argument("--terms", type=int, default=12)
    g.add_argument("--seed", type=int, default=42)
    g.add_argument("--out", default="configs/hamiltonian.json")
    c = sub.add_parser("chem")
    c.add_argument("--molecule", choices=["H2", "HeH+"], default="H2")
    c.add_argument("--r", type=float, default=1.4)
    c.add_argument("--shots", type=int, default=2000)
    c.add_argument("--repeats", type=int, default=200)
    c.add_argument("--out", default="results/chem.csv")
    a = ap.parse_args(argv)
    if a.cmd == "random":
        h = random_hamiltonian(a.qubits, a.terms, np.random.default_rng(a.seed))
        save_hamiltonian_json(h, a.out)
        print(f"{len(h)} terms -> {a.out}")
    elif a.cmd == "run":
        h = load_hamiltonian_json(a.hamiltonian)
        n = len(h[0][1])
        if a.state:
            psi = np.array(json.loads(open(a.state, encoding="utf8").read()), dtype=complex)
            if len(psi) != 2 ** n:
                raise ValueError("state dimension does not match n_qubits")
        else:
            psi = random_state(n, np.random.default_rng(a.seed))
        rows = run_experiment(h, psi, a.shots, a.repeats, a.seed)
        write_csv(rows, a.out)
        _print(rows)
    else:
        sym, ch = (["H", "H"], 0) if a.molecule == "H2" else (["He", "H"], 1)
        h, n = chemistry_hamiltonian(sym, a.r, ch)
        pop = np.array([bin(i).count("1") for i in range(2 ** n)])
        mat = sum(c * _kron_label(p) for c, p in h) + 50.0 * np.diag((pop - 2.0) ** 2)      # keep the 2-electron sector
        w, v = np.linalg.eigh(mat)
        psi = v[:, 0].astype(complex)
        rows = run_experiment(h, psi, a.shots, a.repeats, 0)
        write_csv(rows, a.out)
        _print(rows)


def _kron_label(label):
    from .estimators import pauli_matrix
    return pauli_matrix(label)


if __name__ == "__main__":
    main()
