import json

import numpy as np
import pytest

from qshotlab import cli, experiment, pauli, strategies


def test_expectation_of_basis_states():
    zero = np.array([1, 0], complex)
    plus = np.array([1, 1], complex) / np.sqrt(2)
    assert pauli.expectation(zero, "Z") == pytest.approx(1) and pauli.expectation(zero, "X") == pytest.approx(0)
    assert pauli.expectation(plus, "X") == pytest.approx(1)


def test_random_hamiltonian_is_valid_and_json_loads(tmp_path):
    terms, c = pauli.random_hamiltonian(3, 10, seed=1)
    assert len(set(terms)) == 10 and all(len(t) == 3 and set(t) != {"I"} for t in terms)
    f = tmp_path / "h.json"
    f.write_text(json.dumps({"n_qubits": 3, "terms": [{"pauli": t, "coeff": float(x)} for t, x in zip(terms, c)]}))
    n, t2, c2 = pauli.load_hamiltonian(f)
    assert n == 3 and t2 == terms and np.allclose(c2, c)


def test_allocation_exact_budget_and_minimum():
    s = strategies.allocate([5, 1, 0.01, 3], 100)
    assert s.sum() == 100 and s.min() >= 1 and s[0] > s[1] > s[2]


def test_budget_smaller_than_terms_is_rejected():
    with pytest.raises(AssertionError):
        strategies.allocate([1, 1, 1], 2)


@pytest.mark.parametrize("name", ["uniform", "weighted", "oracle", "adaptive", "shadows"])
def test_estimators_are_unbiased(name):
    rows = experiment.run_instance(3, 8, 3000, reps=400, seed=4, names=[name])
    r = rows[0]
    assert abs(r["bias"]) < 4 * r["bias_se"] + 1e-3


def test_budget_is_used_exactly_by_adaptive():
    terms, c = pauli.random_hamiltonian(3, 8, 0)
    psi = pauli.random_state(3)
    e, _ = pauli.exact_value(psi, terms, c)
    est, s = strategies.adaptive(e, c, 1000, np.random.default_rng(0), 5)
    assert est.shape == (5, 8) and s.sum() == 1000


def test_weighting_beats_uniform_when_coefficients_are_spread():
    rows = []
    for i in range(6):
        rows += experiment.run_instance(3, 10, 2000, reps=300, seed=10 * i, spread=1.2, names=["uniform", "weighted", "oracle"])
    import pandas as pd
    m = pd.DataFrame(rows).groupby("strategy").rmse.mean()
    assert m["weighted"] < m["uniform"] and m["oracle"] <= m["weighted"] * 1.05


def test_rmse_falls_with_budget():
    a = experiment.run_instance(3, 8, 500, reps=300, seed=2, names=["uniform"])[0]["rmse"]
    b = experiment.run_instance(3, 8, 8000, reps=300, seed=2, names=["uniform"])[0]["rmse"]
    assert b < a / 2


def test_cli(tmp_path):
    cli.main(["run", "--instances", "2", "--reps", "40", "--shots", "400", "--out", str(tmp_path)])
    assert (tmp_path / "summary.csv").exists()
