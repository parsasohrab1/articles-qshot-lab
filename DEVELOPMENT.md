# qshot-lab — development status (branch `develop`)
First development pass from the SRS in README.md (2026-10-23).
| SRS item | Status |
|---|---|
| FR-01 JSON Hamiltonian ({n_qubits, terms:[{pauli, coeff}]}) | done with validation (`load_hamiltonian_json`) |
| FR-02 random input states with seed | done (`random_state`; also a state-JSON input in the CLI) |
| FR-03 uniform allocation | done |
| FR-04 weighted allocation (|c_i|, or variance) | done: `weighted` (|c_i|) and `oracle_var` (|c_i| x TRUE sigma_i — unrealistic reference needing the state) |
| FR-05 classical shadows | done (random Pauli bases, vectorised, grouped sampling) |
| FR-06 RMSE and bias to CSV | done (`python -m qshot run ...`) |
| FR-07 plot | not done |
| SRS section 4: adaptive allocation | implemented: pilot round, Beta-posterior variance per term, Neyman allocation of the remaining shots, pooled counts |
| Performance: 3 qubits, 12 terms < 5 s | met (100 repeats of all five methods in the test) |
Tests: `python -m pytest tests -q` (5 pass). Direct methods sample term outcomes from the exact Binomial (equivalent to measuring each term separately).
## Measured (`results/experiments.json`; RMSE of the energy estimate)
**Simulated inputs** (random Pauli Hamiltonians, 3 qubits, 12 terms, random states; mean of 5 draws x 100 repeats):
| shots | uniform | weighted | oracle-var | adaptive | shadows |
|---|---|---|---|---|---|
| 500 | 0.384 | 0.325 | 0.312 | 0.318 | 0.352 |
| 2000 | 0.181 | 0.165 | 0.164 | **0.158** | 0.172 |
| 8000 | 0.100 | 0.080 | 0.078 | 0.079 | 0.090 |
**Real molecular Hamiltonians** (computed STO-3G Hamiltonians of H2 and HeH+ with their exact ground states; 200 repeats; NOT measured data):
| system | shots | uniform | weighted | oracle-var | adaptive | shadows |
|---|---|---|---|---|---|---|
| H2 R=1.4 (15 terms) | 2000 | 0.0103 | 0.0146 | 0.0080 | 0.0120 | 0.0317 |
| H2 R=3.0 (15 terms) | 2000 | 0.0108 | 0.0124 | 0.0078 | 0.0102 | 0.0275 |
| HeH+ R=1.5 (27 terms) | 2000 | 0.0261 | 0.0333 | 0.0194 | 0.0245 | 0.0486 |
**The ranking changes with the input.** On random Hamiltonians weighting and adaptivity beat uniform and shadows are in between. On the real chemistry Hamiltonians shadows are about 2-3x worse than uniform at every budget, |c_i|-weighting is **worse than uniform**, and adaptive allocation ranges from about 9% better than uniform (8000 shots) to 20% worse (500 shots); only the oracle (true variances) is clearly better. This corrects the earlier draft, which had only random Hamiltonians. Few systems (3), one seed per run, no grouping of commuting terms, no noise — **[to be measured]**.
## Real measured data
None found for this question: searches (Zenodo, IEEE DataPort, arXiv) found IBM/Rigetti circuit-histogram datasets and Pauli-circuit datasets with 20,000-shot distributions, but no raw Pauli-term outcomes suitable for comparing shot allocation on a real device. The molecular Hamiltonians are real physics but computed, not experimental. No hardware data were used.
## Synthetic data policy
Random Hamiltonians/states are labelled SIMULATED inputs; real-chemistry results are reported separately.
## Next steps
Commuting-term grouping, FR-07 plots, hardware or noise-model data (needs an account or a published raw-count dataset), more molecules (BeH2 with sparse methods), confidence intervals.
