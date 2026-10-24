# Shot-Budget Lab – measurement allocation for Pauli-sum observables

Reproducible simulation that compares how a fixed budget of shots is split among the terms of a Hamiltonian H = Σ cᵢ Pᵢ:
**uniform**, **weighted** (∝ |cᵢ|), **oracle** (∝ |cᵢ|·σᵢ, needs the true variances; a bound), **adaptive** (pilot round + posterior variance) and **classical shadows**
(random Pauli bases, every snapshot reused for every term). Output: RMSE and bias against the exact value.

**Maturity: TRL 4** (laboratory validation by simulation for 3–8 qubits; no quantum hardware, as the SRS in [docs/SRS.md](docs/SRS.md) states). See [docs/TRL.md](docs/TRL.md).
Moving to TRL 5 needs the same protocol on a noisy simulator or a real device (readout error, qubit-wise commuting groups).

## Pipeline
| SRS | What it does | Where |
|-----|--------------|-------|
| FR-01 | Hamiltonian from JSON `{n_qubits, terms:[{pauli, coeff}]}` or random generator | `pauli.py` |
| FR-02 | Random states with a seed | `pauli.random_state` |
| FR-03 | Uniform allocation | `strategies.uniform` |
| FR-04 | Weighted allocation (∝ \|cᵢ\|), oracle variance weighting, adaptive pilot/posterior allocation | `strategies.py` |
| FR-05 | Classical shadows, exact simulation of all 3ⁿ bases | `strategies.shadows` |
| FR-06 | RMSE and bias per strategy and budget, CSV | `experiment.py`, `cli.py` |

```bash
pip install -e ".[dev]"
pytest -q
qshotlab run --qubits 3 --terms 12 --instances 6 --reps 300 --out results
```

## Lab result (3 qubits, 12 terms, log-normal coefficient spread 1.0, 6 instances × 300 repetitions; not a hardware claim)
RMSE relative to uniform allocation: weighted 0.72–0.73, adaptive 0.71–0.73, oracle 0.70–0.73 at 500, 2 000 and 8 000 shots; classical shadows 1.28–1.35 (worse for these low-weight terms).
With nearly equal coefficients (4 qubits, 20 terms, spread 0.3) the gain shrinks to 0.90–0.95. All estimators are unbiased within sampling error (tests). The adaptive strategy did **not** beat plain weighting here.
