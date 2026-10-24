# TRL assessment – Shot-Budget Lab

Level reached: **TRL 4 – component/subsystem validated in the laboratory (simulation).**

| Evidence | Status |
|----------|--------|
| Working prototype: Hamiltonians, five strategies, experiment runner, CLI | done |
| Automated tests: exact expectations of basis states, JSON loader, allocation sums to the budget, unbiasedness of every strategy (within 4 standard errors), adaptive uses the whole budget, weighting beats uniform on spread coefficients and the oracle is not worse, RMSE falls with budget, CLI | `pytest` – 13 tests, run in CI |
| Validation against known ground truth | exact state-vector expectation values |

Gaps before TRL 5: noisy simulation and hardware runs, commuting-group measurements, larger systems (>8 qubits), the graph-based adaptive method described in the SRS.
