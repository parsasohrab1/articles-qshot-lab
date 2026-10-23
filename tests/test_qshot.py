import time

import numpy as np
import pytest

from qshot import chemistry_hamiltonian, exact_expectation, load_hamiltonian_json, random_hamiltonian, random_state, run_experiment, save_hamiltonian_json, term_expectations
from qshot.estimators import _Rot, estimate_adaptive, estimate_shadow, estimate_uniform, estimate_weighted


def _setup(n=3, terms=12, seed=0):
    rng = np.random.default_rng(seed)
    return random_hamiltonian(n, terms, rng), random_state(n, rng)


def test_json_roundtrip_and_validation(tmp_path):
    h, _ = _setup()
    p = tmp_path / "h.json"
    save_hamiltonian_json(h, p)
    h2 = load_hamiltonian_json(p)
    assert [t[1] for t in h2] == [t[1] for t in h] and np.allclose([t[0] for t in h2], [t[0] for t in h])
    bad = tmp_path / "b.json"
    bad.write_text('{"n_qubits": 2, "terms": [{"pauli": "XQ", "coeff": 1}]}')
    with pytest.raises(ValueError):
        load_hamiltonian_json(bad)


def test_all_estimators_unbiased_on_bell_state():
    bell = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    h = [(1.0, "XX"), (-1.0, "YY"), (0.7, "ZZ"), (0.3, "ZI"), (0.2, "II")]
    exact = exact_expectation(bell, h)
    e = term_expectations(bell, h)
    rot = _Rot(bell)
    for f in (lambda r: estimate_uniform(h, e, 600, r), lambda r: estimate_weighted(h, e, 600, r), lambda r: estimate_adaptive(h, e, 600, r), lambda r: estimate_shadow(h, rot, 600, r)):
        rng = np.random.default_rng(0)
        vals = [f(rng) for _ in range(300)]
        assert abs(np.mean(vals) - exact) < 4 * np.std(vals) / np.sqrt(300) + 2e-3


def test_budget_is_respected_and_rmse_shrinks_with_shots():
    h, psi = _setup()
    a = {r["method"]: r["rmse"] for r in run_experiment(h, psi, 300, repeats=100)}
    b = {r["method"]: r["rmse"] for r in run_experiment(h, psi, 3000, repeats=100)}
    for m in a:
        assert b[m] < a[m], m


def test_performance_3_qubits_12_terms_under_5_seconds():
    h, psi = _setup()
    t = time.time()
    rows = run_experiment(h, psi, 2000, repeats=100)
    assert time.time() - t < 5.0 and {r["method"] for r in rows} == {"uniform", "weighted", "oracle_var", "adaptive", "shadow"}


def test_real_molecular_hamiltonian_h2_is_valid():
    h, n = chemistry_hamiltonian(["H", "H"], 1.4)
    assert n == 4 and len(h) > 5 and all(len(p) == 4 for _, p in h)
