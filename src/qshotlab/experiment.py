import numpy as np
import pandas as pd

from . import pauli, strategies


def run_instance(n_qubits, n_terms, shots, reps=200, seed=0, spread=1.0, names=None):
    terms, c = pauli.random_hamiltonian(n_qubits, n_terms, seed, spread)
    psi = pauli.random_state(n_qubits, seed + 1)
    e, exact = pauli.exact_value(psi, terms, c)
    rng = np.random.default_rng(seed + 2)
    rows = []
    for name in names or strategies.STRATEGIES:
        f = strategies.STRATEGIES[name]
        kw = {"terms": terms, "psi": psi} if name == "shadows" else {}
        est, s = f(e, c, shots, rng, reps, **kw)
        err = est @ c - exact
        rows.append({"strategy": name, "shots": shots, "rmse": float(np.sqrt(np.mean(err ** 2))), "bias": float(err.mean()),
                     "bias_se": float(err.std() / np.sqrt(reps)), "exact": exact})
    return rows


def run(n_qubits=3, n_terms=12, shots_list=(500, 2000, 8000), instances=5, reps=200, seed=0, spread=1.0, names=None):
    rows = []
    for i in range(instances):
        for sh in shots_list:
            rows += [{"instance": i, **r} for r in run_instance(n_qubits, n_terms, sh, reps, seed + 10 * i, spread, names)]
    d = pd.DataFrame(rows)
    s = d.groupby(["shots", "strategy"], as_index=False).agg(rmse=("rmse", "mean"), bias=("bias", "mean"))
    base = s[s.strategy == "uniform"].set_index("shots")["rmse"]
    s["rmse_vs_uniform"] = [r.rmse / base[r.shots] for r in s.itertuples()]
    return d, s
